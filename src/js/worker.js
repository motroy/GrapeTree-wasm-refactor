// worker.js - Web Worker for GrapeTree WASM computations
// Handles CPU-intensive tasks in a background thread

// Load the Emscripten-generated glue code
try {
    importScripts('grapetree.js');
} catch (e) {
    console.warn('Failed to load grapetree.js:', e);
}

let module = null;

// Helper functions for parsing
function detectDelimiter(line) {
    const tabCount = (line.match(/\t/g) || []).length;
    const commaCount = (line.match(/,/g) || []).length;
    return tabCount > commaCount ? '\t' : ',';
}

function parseProfile(content) {
    const lines = content.trim().split('\n');

    if (lines.length < 2) {
        throw new Error('Profile file must have at least a header and one data line');
    }

    // Detect delimiter (tab or comma)
    const delimiter = detectDelimiter(lines[0]);

    // Parse header - handle both '#Strain' and 'Strain' (without # prefix)
    const header = lines[0].split(delimiter);
    let strainColumnName;
    if (header[0].startsWith('#')) {
        strainColumnName = header[0].substring(1).trim();
    } else {
        strainColumnName = header[0].trim();
        // console.warn('Profile file header missing # prefix in first column');
    }
    const geneNames = header.slice(1).map(s => s.trim());

    // Parse data lines
    const strains = [];
    const profiles = [];

    for (let i = 1; i < lines.length; i++) {
        const line = lines[i].trim();
        if (!line || line.startsWith('#')) continue;

        const fields = line.split(delimiter);

        if (fields.length !== header.length) {
            // Skip malformed lines
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

function parseFasta(content) {
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

// Handle messages from main thread
onmessage = async function(e) {
    const { type, id, data } = e.data;

    try {
        if (type === 'init') {
            if (!module) {
                // GrapeTreeWASMModule is defined globally by grapetree.js
                if (typeof GrapeTreeWASMModule === 'undefined') {
                     throw new Error('GrapeTreeWASMModule is not defined. grapetree.js failed to load.');
                }
                module = await GrapeTreeWASMModule();
            }
            postMessage({ type: 'init_done', id });
        }
        else if (type === 'compute_tree') {
            if (!module) throw new Error('WASM module not initialized');

            // The C++ code will call postMessage for progress updates
            // We just need to capture the final result
            const resultJson = module.compute_tree(
                JSON.stringify(data.data),
                data.method,
                data.matrix,
                data.missing,
                data.heuristic
            );

            postMessage({ type: 'result', id, result: resultJson });
        }
        else if (type === 'compute_distance_matrix') {
            if (!module) throw new Error('WASM module not initialized');

            const resultJson = module.compute_distance_matrix(
                JSON.stringify(data.data),
                data.matrixType,
                data.missing
            );

            postMessage({ type: 'result', id, result: resultJson });
        }
        else if (type === 'parse_file') {
            const { content, fileName } = data;
            let result;

            if (fileName.endsWith('.fasta') || fileName.endsWith('.fa') || fileName.endsWith('.fna')) {
                result = parseFasta(content);
            } else {
                // Assume profile
                result = parseProfile(content);
            }

            postMessage({ type: 'result', id, result });
        }
    } catch (error) {
        postMessage({
            type: 'error',
            id,
            error: error.message || String(error)
        });
    }
};
