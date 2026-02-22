/**
 * clustering.js - Clustering algorithms for heatmap ordering
 *
 * Provides functions to reorder distance matrices based on clustering
 * (specifically using MST traversal which is O(N^2) and consistent with GrapeTree).
 */

class Clustering {
    /**
     * Compute a permutation of indices that clusters similar items together.
     * Uses Prim's algorithm to build an MST on the symmetrized matrix,
     * then traverses it.
     *
     * @param {number[][]} matrix - N x N distance matrix
     * @returns {number[]} - Array of size N containing the new order of indices
     */
    static computeMSTOrder(matrix) {
        const n = matrix.length;
        if (n === 0) return [];
        if (n === 1) return [0];

        // 1. Prim's Algorithm to build MST
        // adj[i] contains {to: neighbor_index, dist: distance}
        const adj = Array.from({ length: n }, () => []);

        const inMST = new Uint8Array(n); // 0 or 1
        const minDist = new Float64Array(n).fill(Infinity);
        const parent = new Int32Array(n).fill(-1);

        // Start from node 0
        minDist[0] = 0;

        // Helper to get symmetric distance
        const getDist = (i, j) => {
            const d1 = matrix[i][j];
            const d2 = matrix[j][i];
            // Use average distance for stability, or min
            return (d1 + d2) / 2;
        };

        for (let iter = 0; iter < n; iter++) {
            // Find node u not in MST with min minDist
            let u = -1;
            let minVal = Infinity;

            for (let v = 0; v < n; v++) {
                if (!inMST[v] && minDist[v] < minVal) {
                    minVal = minDist[v];
                    u = v;
                }
            }

            // If graph is disconnected, pick the first unvisited node
            if (u === -1) {
                for (let v = 0; v < n; v++) {
                    if (!inMST[v]) {
                        u = v;
                        minVal = 0; // Start new component
                        break;
                    }
                }
            }
            if (u === -1) break; // All visited

            inMST[u] = 1;

            // Add edge to MST
            if (parent[u] !== -1) {
                const p = parent[u];
                const d = getDist(p, u);
                adj[p].push({ to: u, dist: d });
                adj[u].push({ to: p, dist: d });
            }

            // Update neighbors
            for (let v = 0; v < n; v++) {
                if (!inMST[v]) {
                    const d = getDist(u, v);
                    if (d < minDist[v]) {
                        minDist[v] = d;
                        parent[v] = u;
                    }
                }
            }
        }

        // 2. Traverse MST to get linear order
        // We'll use a DFS. To improve the ordering ("seriation"),
        // we sort neighbors by distance from current node.

        const order = [];
        const visited = new Uint8Array(n);

        // Sort adjacency lists by distance
        for(let i=0; i<n; i++) {
            adj[i].sort((a, b) => a.dist - b.dist);
        }

        // Handle disconnected components if any
        for (let i = 0; i < n; i++) {
            if (!visited[i]) {
                Clustering._dfs(i, adj, visited, order);
            }
        }

        return order;
    }

    static _dfs(u, adj, visited, order) {
        visited[u] = 1;
        order.push(u);
        for (const edge of adj[u]) {
            if (!visited[edge.to]) {
                Clustering._dfs(edge.to, adj, visited, order);
            }
        }
    }
}

export default Clustering;
