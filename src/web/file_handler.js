// file_handler.js - Parse GrapeTree input files
// Handles MLST profiles, FASTA files, and metadata

class FileHandler {
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
            return this.parseFasta(content);
        } else if (fileName.endsWith('.json')) {
            return this.parseJson(content);
        } else {
            // Assume tab or comma delimited profile
            return this.parseProfile(content);
        }
    }

    /**
     * Parse MLST/cgMLST profile file
     * Format:
     * #Strain Gene_1 Gene_2 ...
     * strain1 1 2 3 ...
     * strain2 1 2 4 ...
     */
    parseProfile(content) {
        const lines = content.trim().split('\n');

        if (lines.length < 2) {
            throw new Error('Profile file must have at least a header and one data line');
        }

        // Detect delimiter (tab or comma)
        const delimiter = this._detectDelimiter(lines[0]);

        // Parse header
        const header = lines[0].split(delimiter);
        if (!header[0].startsWith('#')) {
            throw new Error('Profile file must start with # in the first column');
        }

        const strainColumnName = header[0].substring(1).trim();
        const geneNames = header.slice(1).map(s => s.trim());

        // Parse data lines
        const strains = [];
        const profiles = [];

        for (let i = 1; i < lines.length; i++) {
            const line = lines[i].trim();
            if (!line || line.startsWith('#')) continue;

            const fields = line.split(delimiter);

            if (fields.length !== header.length) {
                console.warn(`Line ${i + 1}: Expected ${header.length} fields, got ${fields.length}`);
                continue;
            }

            const strainName = fields[0].trim();
            const profile = fields.slice(1).map(s => {
                const val = s.trim();
                // Convert missing data markers to 0
                if (val === '-' || val === '' || val === 'N/A') {
                    return 0;
                }
                return parseInt(val, 10) || 0;
            });

            strains.push(strainName);
            profiles.push(profile);
        }

        if (strains.length === 0) {
            throw new Error('No valid data lines found in profile file');
        }

        return {
            strains,
            profiles,
            geneNames,
            type: 'profile'
        };
    }

    /**
     * Parse FASTA alignment file
     * Returns sequences converted to numerical profiles using p-distance
     */
    parseFasta(content) {
        const sequences = [];
        const strains = [];

        let currentStrain = null;
        let currentSequence = '';

        const lines = content.split('\n');

        for (const line of lines) {
            const trimmed = line.trim();

            if (trimmed.startsWith('>')) {
                // Save previous sequence
                if (currentStrain !== null) {
                    strains.push(currentStrain);
                    sequences.push(currentSequence);
                }

                // Start new sequence
                currentStrain = trimmed.substring(1).trim().split(/\s+/)[0];
                currentSequence = '';
            } else if (trimmed) {
                currentSequence += trimmed.toUpperCase();
            }
        }

        // Save last sequence
        if (currentStrain !== null) {
            strains.push(currentStrain);
            sequences.push(currentSequence);
        }

        if (sequences.length === 0) {
            throw new Error('No sequences found in FASTA file');
        }

        // Verify all sequences have same length
        const seqLength = sequences[0].length;
        for (let i = 1; i < sequences.length; i++) {
            if (sequences[i].length !== seqLength) {
                throw new Error(
                    `All sequences must have the same length. ` +
                    `Sequence ${strains[i]} has length ${sequences[i].length}, ` +
                    `expected ${seqLength}`
                );
            }
        }

        // Convert sequences to numerical profiles
        // Each position becomes a "gene", each nucleotide an "allele"
        const nucleotideMap = { 'A': 1, 'C': 2, 'G': 3, 'T': 4, '-': 0, 'N': 0 };

        const profiles = sequences.map(seq =>
            Array.from(seq).map(nuc => nucleotideMap[nuc] || 0)
        );

        return {
            strains,
            profiles,
            sequences,
            type: 'fasta'
        };
    }

    /**
     * Parse JSON format (GrapeTree session)
     */
    parseJson(content) {
        try {
            const data = JSON.parse(content);

            if (data.strains && data.profiles) {
                return data;
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
