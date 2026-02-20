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
    // pHierCCval: NMI matrix and silhouette scores
    // -------------------------------------------------------------------------

    /**
     * Compute a Normalised Mutual Information (NMI) matrix between every pair
     * of HC levels.  NMI(A,B) = 2·I(A;B) / (H(A)+H(B)) ∈ [0,1].
     *
     * @param {{levels: string[], assignments: Object<string, number[]>}} hierccResult
     * @returns {number[][]} n_levels × n_levels symmetric matrix (diagonal = 1)
     */
    static computeNMIMatrix(hierccResult) {
        const levels = hierccResult.levels;
        const n = levels.length;
        const mat = Array.from({ length: n }, () => new Array(n).fill(0));
        for (let i = 0; i < n; i++) {
            mat[i][i] = 1.0;
            for (let j = i + 1; j < n; j++) {
                const v = HierCC._nmi(
                    hierccResult.assignments[levels[i]],
                    hierccResult.assignments[levels[j]]
                );
                mat[i][j] = v;
                mat[j][i] = v;
            }
        }
        return mat;
    }

    /**
     * Compute average silhouette score for each HC level using a pre-computed
     * pairwise distance matrix.  Levels with ≤1 or ≥n clusters get NaN.
     *
     * Complexity: O(L · n²) where L = number of levels.  For n > 1 000 this
     * can be slow; callers should warn the user or sample before calling.
     *
     * @param {{levels: string[], assignments: Object<string, number[]>}} hierccResult
     * @param {number[][]} distMatrix  - n × n symmetric distance matrix
     * @returns {Object<string, number>}  level → average silhouette ∈ [-1, 1]
     */
    static computeSilhouetteScores(hierccResult, distMatrix) {
        const n = distMatrix.length;
        const scores = {};

        for (const level of hierccResult.levels) {
            const asgn = hierccResult.assignments[level];
            const nClusters = new Set(asgn).size;

            if (nClusters <= 1 || nClusters >= n) {
                scores[level] = NaN;
                continue;
            }

            let total = 0;
            for (let i = 0; i < n; i++) {
                const ci = asgn[i];
                const clusterSums   = new Map();
                const clusterCounts = new Map();

                for (let j = 0; j < n; j++) {
                    if (i === j) continue;
                    const cj = asgn[j];
                    const d  = distMatrix[i][j];
                    clusterSums.set(cj,   (clusterSums.get(cj)   || 0) + d);
                    clusterCounts.set(cj, (clusterCounts.get(cj) || 0) + 1);
                }

                const aCount = clusterCounts.get(ci) || 0;
                const a = aCount > 0 ? clusterSums.get(ci) / aCount : 0;

                let b = Infinity;
                for (const [cj, sum] of clusterSums) {
                    if (cj === ci) continue;
                    const mean = sum / clusterCounts.get(cj);
                    if (mean < b) b = mean;
                }
                if (!isFinite(b)) b = 0;

                const denom = Math.max(a, b);
                total += denom === 0 ? 0 : (b - a) / denom;
            }
            scores[level] = total / n;
        }
        return scores;
    }

    /**
     * Identify "chosen" HC levels — local maxima of the silhouette curve that
     * also show a meaningful NMI drop relative to the next coarser level.
     * Returns an array of level indices (into hierccResult.levels).
     *
     * @param {number[][]} nmiMatrix
     * @param {Object<string, number>|null} silhouetteScores
     * @param {string[]} levels
     * @returns {number[]}
     */
    static detectChosenLevels(nmiMatrix, silhouetteScores, levels) {
        const n = levels.length;
        const chosen = [];

        if (silhouetteScores) {
            const s = levels.map(l => silhouetteScores[l]);
            for (let i = 0; i < n; i++) {
                if (isNaN(s[i])) continue;
                const prev = i > 0     ? s[i - 1] : -Infinity;
                const next = i < n - 1 ? s[i + 1] : -Infinity;
                if (s[i] >= prev && s[i] >= next && s[i] > 0) {
                    chosen.push(i);
                }
            }
        } else {
            // Fall back to NMI drop: mark levels where NMI with the next level < 0.9
            for (let i = 0; i < n - 1; i++) {
                if (nmiMatrix[i][i + 1] < 0.9) chosen.push(i);
            }
        }
        return chosen;
    }

    // -------------------------------------------------------------------------
    // Private helpers
    // -------------------------------------------------------------------------

    /**
     * Normalised Mutual Information between two integer-label arrays.
     * @private
     */
    static _nmi(a, b) {
        const n = a.length;
        const jointCounts = new Map();
        const countsA     = new Map();
        const countsB     = new Map();

        for (let i = 0; i < n; i++) {
            const ai = a[i], bi = b[i];
            const key = `${ai}_${bi}`;
            jointCounts.set(key, (jointCounts.get(key) || 0) + 1);
            countsA.set(ai, (countsA.get(ai) || 0) + 1);
            countsB.set(bi, (countsB.get(bi) || 0) + 1);
        }

        let mi = 0;
        for (const [key, cnt] of jointCounts) {
            const sep = key.indexOf('_');
            const ai  = key.slice(0, sep);
            const bi  = key.slice(sep + 1);
            const pxy = cnt / n;
            const px  = countsA.get(+ai) / n;
            const py  = countsB.get(+bi) / n;
            if (pxy > 0 && px > 0 && py > 0) {
                mi += pxy * Math.log2(pxy / (px * py));
            }
        }

        let hA = 0, hB = 0;
        for (const c of countsA.values()) { const p = c / n; if (p > 0) hA -= p * Math.log2(p); }
        for (const c of countsB.values()) { const p = c / n; if (p > 0) hB -= p * Math.log2(p); }

        if (hA + hB === 0) return 1.0;
        return Math.min(1.0, Math.max(0.0, (2 * mi) / (hA + hB)));
    }

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
