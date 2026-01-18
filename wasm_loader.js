// wasm_loader.js - WebAssembly module loader for GrapeTree
// Provides JavaScript API for WASM tree computation

class GrapeTreeWASM {
    constructor() {
        this.module = null;
        this.isInitialized = false;
    }
    
    /**
     * Initialize the WASM module
     * @returns {Promise<GrapeTreeWASM>}
     */
    async init() {
        if (this.isInitialized) {
            return this;
        }
        
        try {
            // Load the WASM module (Emscripten generated)
            this.module = await GrapeTreeWASMModule();
            this.isInitialized = true;
            console.log('GrapeTree WASM module initialized successfully');
            return this;
        } catch (error) {
            console.error('Failed to initialize WASM module:', error);
            throw new Error(`WASM initialization failed: ${error.message}`);
        }
    }
    
    /**
     * Compute phylogenetic tree from profile data
     * @param {Object} options - Tree computation options
     * @param {Object} options.data - Profile data {strains: [], profiles: []}
     * @param {string} options.method - Tree method: 'MSTree', 'MSTreeV2', 'NJ'
     * @param {string} options.matrix - Matrix type: 'symmetric', 'asymmetric'
     * @param {number} options.missing - Missing data handler: 0-3
     * @param {string} options.heuristic - Tiebreak heuristic: 'eBurst', 'harmonic'
     * @returns {Object} Tree result with newick, edges, nodes
     */
    computeTree(options) {
        this._checkInitialized();
        
        const {
            data,
            method = 'MSTreeV2',
            matrix = 'asymmetric',
            missing = 0,
            heuristic = 'harmonic'
        } = options;
        
        // Validate inputs
        this._validateTreeOptions(data, method, matrix, missing, heuristic);
        
        try {
            // Convert data to JSON string
            const profileJson = JSON.stringify(data);
            
            // Call WASM function
            const resultJson = this.module.compute_tree(
                profileJson,
                method,
                matrix,
                missing,
                heuristic
            );
            
            // Parse result
            const result = JSON.parse(resultJson);
            
            if (!result.success) {
                throw new Error(result.error || 'Tree computation failed');
            }
            
            return {
                newick: result.newick,
                edges: result.edges,
                nNodes: result.n_nodes,
                nEdges: result.n_edges
            };
            
        } catch (error) {
            console.error('Tree computation error:', error);
            throw error;
        }
    }
    
    /**
     * Compute distance matrix only
     * @param {Object} data - Profile data
     * @param {string} matrixType - 'symmetric' or 'asymmetric'
     * @param {number} missing - Missing data handler
     * @returns {Object} Distance matrix result
     */
    computeDistanceMatrix(data, matrixType = 'symmetric', missing = 0) {
        this._checkInitialized();
        
        try {
            const profileJson = JSON.stringify(data);
            
            const resultJson = this.module.compute_distance_matrix(
                profileJson,
                matrixType,
                missing
            );
            
            const result = JSON.parse(resultJson);
            
            if (!result.success) {
                throw new Error(result.error || 'Distance computation failed');
            }
            
            return {
                matrix: result.matrix,
                strainNames: result.strain_names,
                nStrains: result.n_strains
            };
            
        } catch (error) {
            console.error('Distance matrix computation error:', error);
            throw error;
        }
    }
    
    /**
     * Export tree to Newick format string
     * @param {Object} tree - Tree result from computeTree
     * @returns {string} Newick format string
     */
    exportNewick(tree) {
        return tree.newick;
    }
    
    /**
     * Export distance matrix to PHYLIP format
     * @param {Object} matrixResult - Result from computeDistanceMatrix
     * @returns {string} PHYLIP format string
     */
    exportPhylip(matrixResult) {
        const { matrix, strainNames } = matrixResult;
        const n = strainNames.length;
        
        let output = `${n}\n`;
        
        for (let i = 0; i < n; i++) {
            // PHYLIP format: name (10 chars) followed by distances
            const name = strainNames[i].padEnd(10).substring(0, 10);
            const distances = matrix[i].map(d => d.toFixed(6)).join(' ');
            output += `${name} ${distances}\n`;
        }
        
        return output;
    }
    
    // Private helper methods
    
    _checkInitialized() {
        if (!this.isInitialized || !this.module) {
            throw new Error('WASM module not initialized. Call init() first.');
        }
    }
    
    _validateTreeOptions(data, method, matrix, missing, heuristic) {
        if (!data || !data.strains || !data.profiles) {
            throw new Error('Invalid data format. Expected {strains: [], profiles: []}');
        }
        
        if (data.strains.length === 0) {
            throw new Error('No strains provided');
        }
        
        if (data.profiles.length !== data.strains.length) {
            throw new Error('Number of profiles must match number of strains');
        }
        
        const validMethods = ['MSTree', 'MSTreeV2', 'NJ'];
        if (!validMethods.includes(method)) {
            throw new Error(`Invalid method: ${method}. Must be one of: ${validMethods.join(', ')}`);
        }
        
        const validMatrix = ['symmetric', 'asymmetric'];
        if (!validMatrix.includes(matrix)) {
            throw new Error(`Invalid matrix type: ${matrix}. Must be one of: ${validMatrix.join(', ')}`);
        }
        
        if (missing < 0 || missing > 3) {
            throw new Error('Missing data handler must be 0-3');
        }
        
        const validHeuristics = ['eBurst', 'harmonic'];
        if (!validHeuristics.includes(heuristic)) {
            throw new Error(`Invalid heuristic: ${heuristic}. Must be one of: ${validHeuristics.join(', ')}`);
        }
    }
}

// Export for module systems
if (typeof module !== 'undefined' && module.exports) {
    module.exports = GrapeTreeWASM;
}

// Export as ES6 module
export default GrapeTreeWASM;
