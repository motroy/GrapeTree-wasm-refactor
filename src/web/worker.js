// worker.js - Web Worker for GrapeTree WASM computations
// Handles CPU-intensive tasks in a background thread

// Load the Emscripten-generated glue code
importScripts('grapetree.js');

let module = null;

// Handle messages from main thread
onmessage = async function(e) {
    const { type, id, data } = e.data;

    try {
        if (type === 'init') {
            if (!module) {
                // GrapeTreeWASMModule is defined globally by grapetree.js
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
    } catch (error) {
        postMessage({
            type: 'error',
            id,
            error: error.message || String(error)
        });
    }
};
