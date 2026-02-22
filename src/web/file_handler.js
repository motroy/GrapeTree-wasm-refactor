// file_handler.js - Parse GrapeTree input files
// Handles MLST profiles, FASTA files, and metadata

class FileHandler {
    constructor() {
        this.worker = new Worker('worker.js');
        this.workerResolves = new Map();
        this.worker.onmessage = this._handleWorkerMessage.bind(this);
        this.msgId = 0;
    }

    _handleWorkerMessage(e) {
        const { type, id, result, error } = e.data;
        if (this.workerResolves.has(id)) {
            const { resolve, reject } = this.workerResolves.get(id);
            this.workerResolves.delete(id);
            if (type === 'error') reject(new Error(error));
            else resolve(result);
        }
    }

    async _parseInWorker(content, fileName) {
        const id = this.msgId++;
        return new Promise((resolve, reject) => {
            this.workerResolves.set(id, { resolve, reject });
            this.worker.postMessage({
                type: 'parse_file',
                id,
                data: { content, fileName }
            });
        });
    }

    /**
     * Parse a file and return profile data
     * @param {File} file - Input file
     * @returns {Promise<Object>} Parsed data
     */
    async parse(file) {
        const fileName = file.name.toLowerCase();
        const content = await this._readFile(file);

        if (fileName.endsWith('.fasta') || fileName.endsWith('.fa') ||
            fileName.endsWith('.fna')) {
            return this._parseInWorker(content, fileName);
        } else if (fileName.endsWith('.json')) {
            return this.parseJson(content);
        } else {
            // Assume tab or comma delimited profile
            return this._parseInWorker(content, fileName);
        }
    }

    // parseProfile and parseFasta are now handled in worker.js

    /**
     * Parse JSON format (GrapeTree session or pre-computed lineage)
     */
    parseJson(content) {
        try {
            const data = JSON.parse(content);

            if (data.strains && data.profiles) {
                return data;
            }

            // Pre-computed tree format: {links: [{source, target, distance}, ...]}
            if (data.links && Array.isArray(data.links)) {
                return this._extractFromLinksFormat(data);
            }

            // Try to extract from GrapeTree session format
            if (data.tree && data.metadata) {
                // Convert from session format
                return this._extractFromSession(data);
            }

            throw new Error('JSON file must contain strains and profiles fields');

        } catch (error) {
            throw new Error(`Failed to parse JSON: ${error.message}`);
        }
    }

    /**
     * Extract pre-computed tree from links format:
     * {links: [{source: 0, target: 1, distance: 3.07}, ...]}
     *
     * Optionally supports:
     *   nodes: ["name0", "name1", ...]  — explicit node labels
     *   metadata: {"name0": {col: val}, ...}  — embedded metadata for color-by
     */
    _extractFromLinksFormat(data) {
        if (data.links.length === 0) {
            throw new Error('Links array is empty');
        }

        // Collect all unique node indices
        const nodeSet = new Set();
        for (const link of data.links) {
            nodeSet.add(link.source);
            nodeSet.add(link.target);
        }
        const nodeIndices = Array.from(nodeSet).sort((a, b) => a - b);

        // Map original indices to sequential 0-based indices
        const indexMap = new Map();
        nodeIndices.forEach((idx, i) => { indexMap.set(idx, i); });

        // Use explicit node labels if provided, otherwise fall back to string indices
        const hasLabels = Array.isArray(data.nodes) && data.nodes.length > 0;
        const strains = nodeIndices.map(idx =>
            hasLabels && data.nodes[idx] !== undefined ? String(data.nodes[idx]) : String(idx)
        );

        const edges = data.links.map(link => ({
            from: indexMap.get(link.source),
            to: indexMap.get(link.target),
            distance: link.distance
        }));

        const result = {
            strains,
            type: 'precomputed',
            precomputedTree: {
                edges,
                nNodes: strains.length,
                nEdges: edges.length,
                newick: null
            }
        };

        // Pass through embedded metadata (keyed by node label) so the UI
        // can populate the "Color by" dropdown without a separate metadata file.
        if (data.metadata && typeof data.metadata === 'object') {
            result.embeddedMetadata = data.metadata;
        }

        return result;
    }

    /**
     * Parse metadata file (tab or comma delimited)
     */
    parseMetadata(content) {
        const lines = content.trim().split('\n');

        if (lines.length < 2) {
            throw new Error('Metadata file must have at least a header and one data line');
        }

        const delimiter = this._detectDelimiter(lines[0]);
        const header = lines[0].split(delimiter).map(s => s.trim());

        // Find ID column
        let idColumnIndex = header.findIndex(h =>
            h.toLowerCase() === 'id' || h.toLowerCase() === 'strain'
        );

        if (idColumnIndex === -1) {
            idColumnIndex = 0; // Use first column as ID
        }

        const metadata = {};
        const fields = header.filter((_, i) => i !== idColumnIndex);

        for (let i = 1; i < lines.length; i++) {
            const line = lines[i].trim();
            if (!line) continue;

            const values = line.split(delimiter).map(s => s.trim());
            const id = values[idColumnIndex];

            metadata[id] = {};
            fields.forEach((field, fieldIdx) => {
                const valueIdx = fieldIdx < idColumnIndex ? fieldIdx : fieldIdx + 1;
                metadata[id][field] = values[valueIdx] || '';
            });
        }

        return {
            metadata,
            fields
        };
    }

    /**
     * Parse a pHierCC hierarchical clustering file.
     * Accepts plain-text tab-delimited files or gzip-compressed (.gz) files.
     *
     * Expected format:
     *   #ST_id  HC0  HC2  HC5  ...
     *   1       1    1    1    ...
     *   2       2    1    1    ...
     *
     * Returns:
     *   {
     *     clusterData: { stId: { HC0: '1', HC2: '1', ... }, ... },
     *     hcLevels: ['HC0', 'HC2', 'HC5', ...]
     *   }
     */
    async parseHierCC(file) {
        let content;
        if (file.name.toLowerCase().endsWith('.gz')) {
            content = await this._readGzipFile(file);
        } else {
            content = await this._readFile(file);
        }
        return this._parseHierCCContent(content);
    }

    /**
     * Decompress a gzip file and return its text content.
     * Uses the Web Streams DecompressionStream API (Chrome 80+, Firefox 113+, Safari 16.4+).
     */
    async _readGzipFile(file) {
        const buffer = await file.arrayBuffer();
        const ds = new DecompressionStream('gzip');
        const blob = new Blob([buffer]);
        const stream = blob.stream().pipeThrough(ds);
        const reader = stream.getReader();
        const chunks = [];
        while (true) {
            const { done, value } = await reader.read();
            if (done) break;
            chunks.push(value);
        }
        const totalLength = chunks.reduce((sum, c) => sum + c.length, 0);
        const combined = new Uint8Array(totalLength);
        let offset = 0;
        for (const chunk of chunks) {
            combined.set(chunk, offset);
            offset += chunk.length;
        }
        return new TextDecoder().decode(combined);
    }

    /**
     * Parse tab-delimited pHierCC content into a structured object.
     */
    _parseHierCCContent(content) {
        const lines = content.trim().split('\n');
        if (lines.length < 2) {
            throw new Error('pHierCC file must have a header and at least one data row');
        }

        // Parse header — first column is ST id, rest are HC levels
        const header = lines[0].split('\t').map(h => h.trim());
        // Strip leading '#' from first column; remainder are HC level names
        const hcLevels = header.slice(1).filter(h => h.length > 0);

        if (hcLevels.length === 0) {
            throw new Error('pHierCC file must have at least one HC level column');
        }

        const clusterData = {};
        for (let i = 1; i < lines.length; i++) {
            const line = lines[i].trim();
            if (!line || line.startsWith('#')) continue;

            const fields = line.split('\t');
            const stId = fields[0].trim();
            if (!stId) continue;

            clusterData[stId] = {};
            hcLevels.forEach((level, j) => {
                const val = fields[j + 1] !== undefined ? fields[j + 1].trim() : '';
                clusterData[stId][level] = val;
            });
        }

        if (Object.keys(clusterData).length === 0) {
            throw new Error('No data rows found in pHierCC file');
        }

        return { clusterData, hcLevels };
    }

    // Private helper methods

    async _readFile(file) {
        return new Promise((resolve, reject) => {
            const reader = new FileReader();

            reader.onload = (e) => resolve(e.target.result);
            reader.onerror = (e) => reject(new Error('Failed to read file'));

            reader.readAsText(file);
        });
    }

    _detectDelimiter(line) {
        const tabCount = (line.match(/\t/g) || []).length;
        const commaCount = (line.match(/,/g) || []).length;

        return tabCount > commaCount ? '\t' : ',';
    }

    _extractFromSession(sessionData) {
        // Extract profile data from GrapeTree session JSON
        const strains = Object.keys(sessionData.metadata || {});

        // This would need to reconstruct profiles from tree
        // For now, throw an error
        throw new Error(
            'Converting from GrapeTree session format not yet implemented. ' +
            'Please provide profile or FASTA file.'
        );
    }
}

// Example validation
FileHandler.validateProfileData = function(data) {
    if (!data || !data.strains || !data.profiles) {
        return { valid: false, error: 'Missing strains or profiles' };
    }

    if (data.strains.length === 0) {
        return { valid: false, error: 'No strains provided' };
    }

    if (data.strains.length !== data.profiles.length) {
        return {
            valid: false,
            error: 'Number of strains and profiles must match'
        };
    }

    // Check all profiles have same length
    const profileLength = data.profiles[0].length;
    for (let i = 1; i < data.profiles.length; i++) {
        if (data.profiles[i].length !== profileLength) {
            return {
                valid: false,
                error: `Profile ${i} has different length than profile 0`
            };
        }
    }

    return { valid: true };
};

export default FileHandler;
