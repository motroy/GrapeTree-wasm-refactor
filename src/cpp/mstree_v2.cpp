// mstree_v2.cpp - Improved Minimum Spanning Tree
// Uses Edmond's algorithm for directed graphs + branch recrafting

#ifndef GRAPETREE_MSTREE_V2_H
#define GRAPETREE_MSTREE_V2_H

#include <vector>
#include <limits>
#include <algorithm>
#include <numeric>
#include <set>
#include <queue>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <functional>
#include <cmath>
#include <tuple>

namespace grapetree {

// Edge structure defined in mstree.cpp
struct Edge;

// ---------------------------------------------------------------------------
// contemporary() — direct C++ port of the Numba-compiled function from
// achtman-lab/GrapeTree/module/MSTrees.py.
//
// Decides whether Model A (hub-spoke: shared ancestor for src and tgt) is at
// least as likely as Model B (chain: ancestor → src → tgt).
//
//   a0 = d(candidate → current)   directional distances between a candidate
//   a1 = d(current  → candidate)  node and the current source/target
//   b  = d(current  → far-end)    distance from current node to the far end
//   c  = d(candidate→ far-end)    distance from candidate to the far end
//   n_loci                        number of loci (raw-count scale)
//
// All a0/a1/b/c values must be raw allele-difference counts ∈ [0, n_loci].
// Returns true when Model A (hub-spoke) log-likelihood ≥ Model B (chain).
// ---------------------------------------------------------------------------
static bool contemporary(double a0, double a1,
                         double b,  double c,
                         double n_loci)
{
    // Clamp to [0.5, n_loci - 0.5] to keep log() arguments valid.
    auto clamp = [&](double x) {
        return std::max(0.5, std::min(x, n_loci - 0.5));
    };
    a0 = clamp(a0);  a1 = clamp(a1);
    b  = clamp(b);   c  = clamp(c);

    // Early exits (same as Python).
    if (b >= a0 + c && b >= a1 + c) return false;
    if (b == c)                      return true;

    // Model A (hub-spoke) MLE parameters.
    double s11 = std::sqrt(1.0 - a0 / n_loci);
    double s12 = (2.0 * n_loci - b - c) /
                 (2.0 * std::sqrt(n_loci * (n_loci - a0)));

    // Model B (chain) MLE parameters.
    double v    = 1.0 - ((n_loci - a1) * (n_loci - c) / n_loci
                         + (n_loci - b)) / (2.0 * n_loci);
    double denom = b - 2.0 * n_loci * v;
    double s21  = 1.0 + a1 * v / denom;
    double s22  = 1.0 + c  * v / denom;

    // Guard against non-positive arguments to log().
    auto safe_log = [](double x) {
        return std::log(std::max(x, 1e-300));
    };

    double p1 =   a0             * safe_log(1.0 - s11 * s11)
               + (n_loci - a0)   * safe_log(s11 * s11)
               + (b + c)         * safe_log(1.0 - s11 * s12)
               + (2.0 * n_loci - b - c) * safe_log(s11 * s12);

    double p2 =   a1             * safe_log(1.0 - s21)
               + (n_loci - a1)   * safe_log(s21)
               +  b              * safe_log(1.0 - s21 * s22)
               + (n_loci - b)    * safe_log(s21 * s22)
               +  c              * safe_log(1.0 - s22)
               + (n_loci - c)    * safe_log(s22);

    return p1 >= p2;
}

class MSTreeV2 {
private:
    int    n_nodes_;
    int    n_loci_;          // number of loci; used to scale [0,1] ratios to
                             // raw counts when calling contemporary()
    std::vector<double> distance_matrix_;

public:
    // Preferred constructor: accepts an already-flat 1-D distance matrix
    // (row-major, size n*n) and takes ownership via move.  No extra copy.
    MSTreeV2(std::vector<double>&& flat_distances, int n, int n_loci = 0)
        : n_nodes_(n),
          n_loci_(n_loci > 0 ? n_loci : n),
          distance_matrix_(std::move(flat_distances)) {}

    // Convenience constructor from 2-D matrix (kept for compatibility).
    explicit MSTreeV2(
        const std::vector<std::vector<double>>& distances,
        int n_loci = 0
    ) : n_nodes_(static_cast<int>(distances.size())),
        n_loci_(n_loci > 0 ? n_loci
                            : static_cast<int>(distances.size()))
    {
        distance_matrix_.reserve(n_nodes_ * n_nodes_);
        for (const auto& row : distances) {
            distance_matrix_.insert(distance_matrix_.end(),
                                    row.begin(), row.end());
        }
    }

public:
    std::vector<Edge> compute(std::function<void(double)> progress_cb = nullptr) {
        if (progress_cb) progress_cb(0.0);

        // Phase 1: Find minimum incoming edge for each node.
        std::vector<Edge> min_incoming = find_minimum_incoming_edges(progress_cb);

        if (progress_cb) progress_cb(50.0);

        // Phase 2: Detect cycles.
        std::vector<int> cycle_id = detect_cycles(min_incoming);

        if (progress_cb) progress_cb(60.0);

        // Phase 3: Contract cycles if present.
        if (has_cycles(cycle_id)) {
            min_incoming = contract_and_solve(min_incoming, cycle_id);
        }

        if (progress_cb) progress_cb(80.0);

        // Phase 4: Local branch recrafting (Algorithm 1, Zhou et al. 2018).
        recraft_branches(min_incoming);

        if (progress_cb) progress_cb(100.0);

        return min_incoming;
    }

private:
    // -----------------------------------------------------------------------
    // Phase 1 helpers
    // -----------------------------------------------------------------------

    // Harmonic-mean score for a node: lower = more central.
    // Matches Python: weights = N / sum(1/(d+0.1))  re-ranked to [0,1].
    // For tiebreaking we only need the raw harmonic mean (not the percentile
    // ranking), so we return n / sum(1/(d+0.1)).  Lower is more central.
    double harmonic_mean_score(int node) const {
        double sum = 0.0;
        int count  = 0;
        for (int i = 0; i < n_nodes_; ++i) {
            if (i == node) continue;
            double d = distance_matrix_[node * n_nodes_ + i];
            sum += 1.0 / (d + 0.1);
            count++;
        }
        return count > 0 ? static_cast<double>(count) / sum : 0.0;
    }

    // Find minimum incoming edge for each node using harmonic-mean tiebreak.
    std::vector<Edge> find_minimum_incoming_edges(
        std::function<void(double)> progress_cb = nullptr)
    {
        std::vector<Edge> edges;
        int report_frequency = std::max(1, n_nodes_ / 100);

        // Node 0 is the root (no incoming edge).
        for (int to = 1; to < n_nodes_; ++to) {
            if (progress_cb && to % report_frequency == 0)
                progress_cb(static_cast<double>(to) / n_nodes_ * 50.0);

            double min_dist    = std::numeric_limits<double>::max();
            int    best_from   = -1;
            double best_score  = std::numeric_limits<double>::max();

            for (int from = 0; from < n_nodes_; ++from) {
                if (from == to) continue;

                double dist = distance_matrix_[from * n_nodes_ + to];

                if (dist < min_dist) {
                    min_dist   = dist;
                    best_from  = from;
                    best_score = harmonic_mean_score(from);
                } else if (std::abs(dist - min_dist) < 1e-10) {
                    // Tiebreak: prefer source with smaller ht (more central).
                    // Paper sorts edges in ascending order of ht(u); smaller = preferred.
                    double score = harmonic_mean_score(from);
                    if (score < best_score) {
                        best_from  = from;
                        best_score = score;
                    }
                }
            }

            if (best_from != -1)
                edges.emplace_back(best_from, to, min_dist);
        }

        return edges;
    }

    // -----------------------------------------------------------------------
    // Phase 2–3: cycle detection and contraction
    // -----------------------------------------------------------------------

    std::vector<int> detect_cycles(const std::vector<Edge>& edges) {
        std::vector<int> parent(n_nodes_);
        std::iota(parent.begin(), parent.end(), 0);

        std::vector<int> cycle_id(n_nodes_, -1);
        int next_cycle_id = 0;

        for (const Edge& e : edges) {
            int root_from = find_root(parent, e.from);
            int root_to   = find_root(parent, e.to);

            if (root_from == root_to && cycle_id[e.to] == -1) {
                mark_cycle(edges, e.to, cycle_id, next_cycle_id);
                next_cycle_id++;
            }

            parent[root_to] = root_from;
        }

        return cycle_id;
    }

    int find_root(std::vector<int>& parent, int node) {
        int root = node;
        while (parent[root] != root) root = parent[root];
        while (parent[node] != root) {
            int next = parent[node];
            parent[node] = root;
            node = next;
        }
        return root;
    }

    void mark_cycle(const std::vector<Edge>& edges, int start,
                    std::vector<int>& cycle_id, int id)
    {
        int current = start;
        std::set<int> visited;

        while (visited.find(current) == visited.end()) {
            visited.insert(current);
            cycle_id[current] = id;

            bool found = false;
            for (const Edge& e : edges) {
                if (e.to == current) { current = e.from; found = true; break; }
            }
            if (!found) break;
        }
    }

    bool has_cycles(const std::vector<int>& cycle_id) {
        for (int id : cycle_id) if (id != -1) return true;
        return false;
    }

    std::vector<Edge> contract_and_solve(
        const std::vector<Edge>& edges,
        const std::vector<int>&  cycle_id)
    {
        std::map<int, int> node_mapping;
        std::vector<std::set<int>> cycles;
        int next_node = 0;

        std::set<int> unique_cycles;
        for (int id : cycle_id) if (id != -1) unique_cycles.insert(id);
        cycles.resize(unique_cycles.size());

        for (int i = 0; i < n_nodes_; ++i) {
            if (cycle_id[i] == -1) {
                node_mapping[i] = next_node++;
            } else {
                cycles[cycle_id[i]].insert(i);
            }
        }

        for (size_t i = 0; i < cycles.size(); ++i) {
            int contracted_node = next_node++;
            for (int node : cycles[i]) node_mapping[node] = contracted_node;
        }

        int new_size = next_node;
        std::vector<double> new_distances(
            new_size * new_size,
            std::numeric_limits<double>::max());

        std::map<std::pair<int,int>, Edge> edge_mapping;

        std::vector<double> cycle_edge_weights(n_nodes_, 0.0);
        for (const auto& e : edges)
            cycle_edge_weights[e.to] = e.distance;

        for (int i = 0; i < n_nodes_; ++i) {
            for (int j = 0; j < n_nodes_; ++j) {
                if (i == j) continue;

                int ni = node_mapping[i];
                int nj = node_mapping[j];

                if (ni != nj) {
                    double dist         = distance_matrix_[i * n_nodes_ + j];
                    double reduced_dist = dist;
                    if (cycle_id[j] != -1)
                        reduced_dist -= cycle_edge_weights[j];

                    if (reduced_dist < new_distances[ni * new_size + nj]) {
                        new_distances[ni * new_size + nj] = reduced_dist;
                        edge_mapping.erase({ni, nj});
                        edge_mapping.insert({{ni, nj}, Edge(i, j, dist)});
                    }
                }
            }
        }

        MSTreeV2 contracted_solver(std::move(new_distances), new_size, n_loci_);
        std::vector<Edge> contracted_edges = contracted_solver.compute();

        std::vector<Edge> final_edges;
        std::set<int> nodes_with_incoming;

        for (const auto& e : contracted_edges) {
            if (edge_mapping.count({e.from, e.to})) {
                Edge original = edge_mapping.at({e.from, e.to});
                final_edges.push_back(original);
                nodes_with_incoming.insert(original.to);
            }
        }

        for (const auto& e : edges) {
            if (!nodes_with_incoming.count(e.to)) {
                final_edges.push_back(e);
                nodes_with_incoming.insert(e.to);
            }
        }

        return final_edges;
    }

    // -----------------------------------------------------------------------
    // Helper: Reroot the component containing 'old_root' to be rooted at
    // 'new_root', ensuring the directed tree structure is maintained.
    // This reverses edges along the path from old_root to new_root.
    // -----------------------------------------------------------------------
    void reroot(int old_root, int new_root,
                std::vector<Edge>& branches,
                std::map<std::pair<int, int>, int>& edge_map,
                const std::unordered_map<int, std::vector<int>>& adj)
    {
        if (old_root == new_root) return;

        // BFS to find path and parent pointers for backtracking
        std::unordered_map<int, int> parent;
        std::queue<int> q;
        q.push(old_root);
        parent[old_root] = -1;

        bool found = false;
        while (!q.empty()) {
            int u = q.front(); q.pop();
            if (u == new_root) { found = true; break; }
            if (adj.count(u)) {
                for (int v : adj.at(u)) {
                    if (parent.find(v) == parent.end()) {
                        parent[v] = u;
                        q.push(v);
                    }
                }
            }
        }

        if (!found) return; // Should generally not happen if connected

        // Reverse edges along the path from new_root up to old_root
        int curr = new_root;
        while (curr != old_root) {
            int p = parent[curr];
            // The edge currently exists as (p -> curr) in the tree.
            // We want to flip it to (curr -> p).

            // Find the edge index using the map (stored as {p, curr})
            // Since we know p is the parent, and branches store parent->child,
            // we look for {p, curr}.
            auto key = std::make_pair(p, curr);
            if (edge_map.find(key) == edge_map.end()) {
                // Try reverse key just in case, though it shouldn't happen in valid dMST
                key = std::make_pair(curr, p);
            }

            if (edge_map.count(key)) {
                int idx = edge_map[key];
                // Flip the edge in branches
                branches[idx].from = curr;
                branches[idx].to   = p;

                // Update map: remove old key, insert new key {curr, p}
                // (Note: we store directed pairs in edge_map to match branches)
                edge_map.erase(key);
                edge_map[{curr, p}] = idx;
            }
            curr = p;
        }
    }

    // -----------------------------------------------------------------------
    // Phase 4: local branch recrafting
    //
    // Faithful C++ port of _branch_recraft() from
    // achtman-lab/GrapeTree/module/MSTrees.py.
    //
    // Distance matrix stores [0,1] ratios; multiply by n_loci_ to obtain the
    // raw allele-count scale expected by contemporary().
    // -----------------------------------------------------------------------
    void recraft_branches(std::vector<Edge>& branches) {
        if (branches.empty()) return;

        const double nl = static_cast<double>(n_loci_);

        // Scale a [0,1] ratio to a raw allele count for contemporary().
        auto raw = [&](double ratio) { return ratio * nl; };

        // Distance accessor.
        auto D = [&](int i, int j) -> double {
            return distance_matrix_[i * n_nodes_ + j];
        };

        // Precompute harmonic-mean weights (lower = more central).
        std::vector<double> weights(n_nodes_);
        for (int i = 0; i < n_nodes_; ++i) weights[i] = harmonic_mean_score(i);

        // Sort branches: primary key = edge length, secondary = sorted pair
        // of endpoint weights.  Matches Python:
        //   key = [dist[src,tgt]] + sorted([weights[src], weights[tgt]])
        std::sort(branches.begin(), branches.end(),
                  [&](const Edge& a, const Edge& b_e) {
                      if (std::abs(a.distance - b_e.distance) > 1e-12)
                          return a.distance < b_e.distance;
                      double a_min = std::min(weights[a.from], weights[a.to]);
                      double a_max = std::max(weights[a.from], weights[a.to]);
                      double b_min = std::min(weights[b_e.from], weights[b_e.to]);
                      double b_max = std::max(weights[b_e.from], weights[b_e.to]);
                      if (std::abs(a_min - b_min) > 1e-12) return a_min < b_min;
                      return a_max < b_max;
                  });

        // group_id[v] = current group representative for node v.
        std::vector<int> group_id(n_nodes_);
        std::iota(group_id.begin(), group_id.end(), 0);

        // groups[rep] = list of all nodes in that component.
        std::unordered_map<int, std::vector<int>> groups;
        // childrens[v]  = tree-adjacent nodes of v (undirected, built lazily).
        std::unordered_map<int, std::vector<int>> childrens;

        // Map directed edge (u, v) -> index in branches.
        // Only tracks edges that have been processed and committed to a component.
        std::map<std::pair<int,int>, int> edge_map;

        for (const Edge& e : branches) {
            for (int v : {e.from, e.to}) {
                if (!groups.count(v)) {
                    groups[v]    = {v};
                    childrens[v] = {};
                }
            }
        }

        // Helper: sort a list of (weight, dist_to_far, node) triples.
        using Triple = std::tuple<double, double, int>;
        auto sort3 = [](std::vector<Triple>& v) { std::sort(v.begin(), v.end()); };

        size_t i = 0;
        while (i < branches.size()) {
            int src = branches[i].from;
            int tgt = branches[i].to;

            // Snapshot the two component lists before any src/tgt changes.
            std::vector<int> sources = groups[group_id[src]];
            std::vector<int> targets = groups[group_id[tgt]];

            // Nodes visited during this branch's recrafting pass.
            std::unordered_set<int> tried;

            double src_tgt = D(src, tgt);

            // ------------------------------------------------------------------
            // Source-side recrafting: find a better source inside sources.
            // ------------------------------------------------------------------
            if (sources.size() > 1) {
                // Try the 3 lightest (most central) members of the source group.
                std::vector<Triple> cands;
                cands.reserve(sources.size());
                for (int s : sources)
                    cands.emplace_back(weights[s], D(s, tgt), s);
                sort3(cands);

                int limit = std::min((int)cands.size(), 3);
                for (int k = 0; k < limit; ++k) {
                    auto [ws, ds_tgt, s] = cands[k];
                    if (s == src) break;
                    if (ds_tgt < 1.5 * src_tgt) {
                        if (contemporary(raw(D(s, src)), raw(D(src, s)),
                                         raw(ds_tgt), raw(src_tgt), nl)) {
                            tried.insert(src);
                            src = s;
                            break;
                        }
                    }
                }

                // Walk children of current src looking for a closer relay.
                while (!tried.count(src)) {
                    tried.insert(src);

                    std::vector<Triple> mid;
                    for (int s : childrens[src]) {
                        if (!tried.count(s) && D(s, tgt) < 2.0 * src_tgt)
                            mid.emplace_back(weights[s], D(s, tgt), s);
                    }
                    sort3(mid);

                    bool moved = false;
                    for (auto& [w, d, s] : mid) {
                        if (d < src_tgt) {
                            if (!contemporary(raw(D(src, s)), raw(D(s, src)),
                                              raw(src_tgt), raw(d), nl)) {
                                tried.insert(src);
                                src = s;
                                moved = true;
                                break;
                            }
                        } else if (w < weights[src]) {
                            if (contemporary(raw(D(s, src)), raw(D(src, s)),
                                             raw(d), raw(src_tgt), nl)) {
                                tried.insert(src);
                                src = s;
                                moved = true;
                                break;
                            }
                        }
                        tried.insert(s);
                    }
                    if (!moved) break;
                }
            }

            // Refresh src→tgt distance after possible src change.
            src_tgt = D(src, tgt);

            // ------------------------------------------------------------------
            // Target-side recrafting: find a better target inside targets.
            // ------------------------------------------------------------------
            if (targets.size() > 1) {
                std::vector<Triple> cands;
                cands.reserve(targets.size());
                for (int t : targets)
                    cands.emplace_back(weights[t], D(src, t), t);
                sort3(cands);

                int limit = std::min((int)cands.size(), 3);
                for (int k = 0; k < limit; ++k) {
                    auto [wt, d_src_t, t] = cands[k];
                    if (t == tgt) break;
                    if (d_src_t < 1.5 * src_tgt) {
                        if (contemporary(raw(D(t, tgt)), raw(D(tgt, t)),
                                         raw(d_src_t), raw(src_tgt), nl)) {
                            tried.insert(tgt);
                            // Reroot the target component to 't' to allow incoming edge from 'src'
                            reroot(tgt, t, branches, edge_map, childrens);
                            tgt = t;
                            break;
                        }
                    }
                }

                while (!tried.count(tgt)) {
                    tried.insert(tgt);

                    std::vector<Triple> mid;
                    for (int t : childrens[tgt]) {
                        if (!tried.count(t) && D(src, t) < 2.0 * src_tgt)
                            mid.emplace_back(weights[t], D(src, t), t);
                    }
                    sort3(mid);

                    bool moved = false;
                    for (auto& [w, d, t] : mid) {
                        if (d < src_tgt) {
                            if (!contemporary(raw(D(tgt, t)), raw(D(t, tgt)),
                                              raw(src_tgt), raw(d), nl)) {
                                tried.insert(tgt);
                                reroot(tgt, t, branches, edge_map, childrens);
                                tgt = t;
                                moved = true;
                                break;
                            }
                        } else if (w < weights[tgt]) {
                            if (contemporary(raw(D(t, tgt)), raw(D(tgt, t)),
                                             raw(d), raw(src_tgt), nl)) {
                                tried.insert(tgt);
                                reroot(tgt, t, branches, edge_map, childrens);
                                tgt = t;
                                moved = true;
                                break;
                            }
                        }
                        tried.insert(t);
                    }
                    if (!moved) break;
                }
            }

            double brlen = D(src, tgt);
            branches[i] = Edge(src, tgt, brlen);

            // If the updated edge is still the shortest among remaining
            // branches, commit it and merge groups; otherwise re-sort.
            bool needs_resort = (i + 1 < branches.size() &&
                                 branches[i + 1].distance < brlen);

            if (!needs_resort) {
                // Record the finalized edge in our map so we can find/flip it later if needed.
                edge_map[{src, tgt}] = i;

                // Merge tgt's group into src's group.
                int tid = group_id[tgt];
                int sid = group_id[src];
                if (tid != sid) {
                    for (int t : groups[tid])
                        group_id[t] = sid;
                    auto& sg = groups[sid];
                    auto& tg = groups[tid];
                    sg.insert(sg.end(), tg.begin(), tg.end());
                    groups.erase(tid);
                }
                childrens[src].push_back(tgt);
                childrens[tgt].push_back(src);
                ++i;
            } else {
                // Only re-sort from position i onward by edge length.
                std::sort(branches.begin() + i, branches.end(),
                          [](const Edge& a, const Edge& b_e) {
                              return a.distance < b_e.distance;
                          });
            }
        }
    }
};

} // namespace grapetree

#endif // GRAPETREE_MSTREE_V2_H
