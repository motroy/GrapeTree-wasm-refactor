// hiercc.js - Hierarchical Cluster Complexes (pHierCC) for GrapeTree
//
// Implements the core pHierCC algorithm from:
//   Zhou et al. (2021) "A genomic epidemiological framework based on core-genome
//   MLST for Neisseria gonorrhoeae surveillance"
//   https://github.com/zheminzhou/pHierCC
//
// The algorithm builds on a pre-computed minimum spanning tree (MST):
//   For each distance threshold T ("HC level"), nodes connected by MST edges
//   with distance ≤ T belong to the same cluster. Union-Find over the edge list
//   makes each level O(n·α(n)) ≈ O(n).
//
// This JS implementation operates on the edge list returned by the WASM
// tree computation, so no WASM rebuild is required.

class HierCC {
    /**
     * Standard pHierCC distance thresholds.
     * These match the commonly reported "stable" HC levels in the original tool.
     */
    static DEFAULT_LEVELS = [0, 2, 5, 10, 20, 50, 100, 200, 400, 800, 1500];

    /**
     * Compute hierarchical cluster assignments from an MST edge list.
     *
     * @param {Array<{from: number, to: number, distance: number}>} edges
     *   MST edges (0-based strain indices, integer-valued allelic distances).
     * @param {number} nStrains - Total number of strains.
     * @param {number[]} [levels] - Distance thresholds to evaluate.
     *   Defaults to HierCC.DEFAULT_LEVELS.
     * @returns {{
     *   levels: string[],
     *   assignments: Object<string, number[]>,
     *   clusterCounts: Object<string, number>
     * }}
     *   `levels` — ordered HC level names, e.g. ['HC0','HC2','HC5',...]
     *   `assignments` — level name → array of 1-based cluster IDs (one per strain)
     *   `clusterCounts` — level name → number of distinct clusters
     */
    static compute(edges, nStrains, levels = HierCC.DEFAULT_LEVELS) {
        const levelNames = levels.map(t => `HC${t}`);
        const assignments = {};
        const clusterCounts = {};

        for (let li = 0; li < levels.length; li++) {
            const name = levelNames[li];
            const ids = HierCC._computeAtThreshold(edges, nStrains, levels[li]);
            assignments[name] = ids;
            clusterCounts[name] = new Set(ids).size;
        }

        return { levels: levelNames, assignments, clusterCounts };
    }

    /**
     * Build a clusterData dict (matching the pHierCC file-upload format) from
     * computed HierCC results.
     *
     * @param {string[]} strains - Ordered strain names.
     * @param {{levels: string[], assignments: Object}} hierccResult
     * @returns {{ clusterData: Object, hcLevels: string[] }}
     */
    static buildClusterData(strains, hierccResult) {
        const clusterData = {};
        strains.forEach((strain, i) => {
            clusterData[strain] = {};
            for (const level of hierccResult.levels) {
                clusterData[strain][level] = String(hierccResult.assignments[level][i]);
            }
        });
        return { clusterData, hcLevels: hierccResult.levels };
    }

    /**
     * Return only the HC levels that have at least `minClusters` distinct
     * clusters and fewer than `maxClusters`. Useful for trimming meaningless
     * boundary levels (everything-in-one or each-strain-alone).
     *
     * @param {{levels: string[], clusterCounts: Object}} hierccResult
     * @param {number} nStrains
     * @param {number} [minClusters=1]
     * @param {number} [maxClusters] - defaults to nStrains
     * @returns {string[]} Filtered level names
     */
    static meaningfulLevels(hierccResult, nStrains, minClusters = 1, maxClusters = null) {
        const max = maxClusters !== null ? maxClusters : nStrains;
        return hierccResult.levels.filter(level => {
            const c = hierccResult.clusterCounts[level];
            return c >= minClusters && c <= max;
        });
    }

    // -------------------------------------------------------------------------
    // Private helpers
    // -------------------------------------------------------------------------

    /**
     * Union-Find with path compression and union-by-rank.
     * @private
     */
    static _makeUnionFind(n) {
        const parent = Array.from({ length: n }, (_, i) => i);
        const rank   = new Int32Array(n);

        function find(x) {
            // Iterative path compression (avoids call-stack overflow for large n)
            let root = x;
            while (parent[root] !== root) root = parent[root];
            while (parent[x] !== root) {
                const next = parent[x];
                parent[x] = root;
                x = next;
            }
            return root;
        }

        function unite(a, b) {
            const ra = find(a), rb = find(b);
            if (ra === rb) return;
            if (rank[ra] < rank[rb]) {
                parent[ra] = rb;
            } else if (rank[ra] > rank[rb]) {
                parent[rb] = ra;
            } else {
                parent[rb] = ra;
                rank[ra]++;
            }
        }

        return { find, unite };
    }

    /**
     * Compute 1-based cluster IDs at a single distance threshold.
     *
     * Cluster IDs are assigned in ascending order of each component's minimum
     * strain index, giving stable and reproducible IDs across runs.
     *
     * @private
     */
    static _computeAtThreshold(edges, nStrains, threshold) {
        const uf = HierCC._makeUnionFind(nStrains);

        // Merge strains connected by short-enough MST edges
        for (const { from, to, distance } of edges) {
            if (distance <= threshold + 1e-9) {
                uf.unite(from, to);
            }
        }

        // Find the minimum strain index (representative) in each component
        const rootToMin = new Map();
        for (let i = 0; i < nStrains; i++) {
            const root = uf.find(i);
            const cur  = rootToMin.get(root);
            if (cur === undefined || i < cur) rootToMin.set(root, i);
        }

        // Sort components by representative (min index) → stable 1-based IDs
        const sortedRoots = Array.from(rootToMin.entries())
            .sort((a, b) => a[1] - b[1])   // sort by min strain index
            .map(([root]) => root);

        const rootToId = new Map(sortedRoots.map((root, idx) => [root, idx + 1]));

        return Array.from({ length: nStrains }, (_, i) => rootToId.get(uf.find(i)));
    }
}

// Export for module systems
if (typeof module !== 'undefined' && module.exports) {
    module.exports = HierCC;
}

export default HierCC;
