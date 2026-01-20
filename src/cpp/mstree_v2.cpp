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
#include <utility>

namespace grapetree {

// Edge structure defined in mstree.cpp
struct Edge;

class MSTreeV2 {
private:
    int n_nodes_;
    std::vector<std::vector<double>> distance_matrix_;
    
public:
    explicit MSTreeV2(
        const std::vector<std::vector<double>>& distances
    ) : distance_matrix_(distances),
        n_nodes_(distances.size()) {}
    
    std::vector<Edge> compute() {
        // Phase 1: Find minimum incoming edge for each node
        std::vector<Edge> min_incoming = find_minimum_incoming_edges();
        
        // Phase 2: Detect cycles
        std::vector<int> cycle_id = detect_cycles(min_incoming);
        
        // Phase 3: Contract cycles if present
        if (has_cycles(cycle_id)) {
            min_incoming = contract_and_solve(min_incoming, cycle_id);
        }
        
        // Phase 4: Local branch recrafting optimization
        recraft_branches(min_incoming);
        
        return min_incoming;
    }
    
private:
    // Find minimum incoming edge for each node using harmonic mean tiebreak
    std::vector<Edge> find_minimum_incoming_edges() {
        std::vector<Edge> edges;
        
        // Node 0 is the root (no incoming edge)
        for (int to = 1; to < n_nodes_; ++to) {
            double min_dist = std::numeric_limits<double>::max();
            int best_from = -1;
            double best_score = -1.0;
            
            for (int from = 0; from < n_nodes_; ++from) {
                if (from == to) continue;
                
                double dist = distance_matrix_[from][to];
                
                if (dist < min_dist) {
                    min_dist = dist;
                    best_from = from;
                    best_score = harmonic_mean_score(from);
                } else if (std::abs(dist - min_dist) < 1e-10) {
                    // Tiebreak using harmonic mean
                    double score = harmonic_mean_score(from);
                    if (score > best_score) {
                        best_from = from;
                        best_score = score;
                    }
                }
            }
            
            if (best_from != -1) {
                edges.emplace_back(best_from, to, min_dist);
            }
        }
        
        return edges;
    }
    
    double harmonic_mean_score(int node) {
        double sum = 0.0;
        int count = 0;
        
        for (int i = 0; i < n_nodes_; ++i) {
            if (i == node) continue;
            
            double dist = distance_matrix_[node][i];
            if (dist > 0.0) {
                sum += 1.0 / dist;
                count++;
            }
        }
        
        return count > 0 ? static_cast<double>(count) / sum : 0.0;
    }
    
    // Detect cycles using Union-Find
    std::vector<int> detect_cycles(const std::vector<Edge>& edges) {
        std::vector<int> parent(n_nodes_);
        std::iota(parent.begin(), parent.end(), 0);
        
        std::vector<int> cycle_id(n_nodes_, -1);
        int next_cycle_id = 0;
        
        for (const Edge& e : edges) {
            int root_from = find_root(parent, e.from);
            int root_to = find_root(parent, e.to);
            
            if (root_from == root_to && cycle_id[e.to] == -1) {
                // Found a cycle - mark all nodes in it
                mark_cycle(edges, e.to, cycle_id, next_cycle_id);
                next_cycle_id++;
            }
            
            parent[root_to] = root_from;
        }
        
        return cycle_id;
    }
    
    int find_root(std::vector<int>& parent, int node) {
        if (parent[node] != node) {
            parent[node] = find_root(parent, parent[node]);
        }
        return parent[node];
    }
    
    void mark_cycle(
        const std::vector<Edge>& edges,
        int start,
        std::vector<int>& cycle_id,
        int id
    ) {
        int current = start;
        std::set<int> visited;
        
        while (visited.find(current) == visited.end()) {
            visited.insert(current);
            cycle_id[current] = id;
            
            // Find edge pointing to current
            bool found = false;
            for (const Edge& e : edges) {
                if (e.to == current) {
                    current = e.from;
                    found = true;
                    break;
                }
            }
            
            if (!found) break;
        }
    }
    
    bool has_cycles(const std::vector<int>& cycle_id) {
        for (int id : cycle_id) {
            if (id != -1) return true;
        }
        return false;
    }
    
    // Contract cycles and recursively solve
    std::vector<Edge> contract_and_solve(
        const std::vector<Edge>& edges,
        const std::vector<int>& cycle_id
    ) {
        // Map old nodes to contracted nodes
        std::map<int, int> node_mapping;
        std::vector<std::set<int>> cycles;
        int next_node = 0;
        
        // Identify unique cycles
        std::set<int> unique_cycles;
        for (int id : cycle_id) {
            if (id != -1) unique_cycles.insert(id);
        }
        
        cycles.resize(unique_cycles.size());
        
        // Build mapping
        for (int i = 0; i < n_nodes_; ++i) {
            if (cycle_id[i] == -1) {
                node_mapping[i] = next_node++;
            } else {
                cycles[cycle_id[i]].insert(i);
            }
        }
        
        // Assign contracted node for each cycle
        for (size_t i = 0; i < cycles.size(); ++i) {
            int contracted_node = next_node++;
            for (int node : cycles[i]) {
                node_mapping[node] = contracted_node;
            }
        }
        
        // Build contracted distance matrix and edge mapping
        int new_size = next_node;
        std::vector<std::vector<double>> new_distances(
            new_size,
            std::vector<double>(new_size, 
                               std::numeric_limits<double>::max())
        );
        
        // Map: (contracted_from, contracted_to) -> Original Edge
        // We use int pair for key as Edge might not be comparable
        std::map<std::pair<int, int>, Edge> edge_mapping;

        for (int i = 0; i < n_nodes_; ++i) {
            for (int j = 0; j < n_nodes_; ++j) {
                if (i == j) continue;
                
                int ni = node_mapping[i];
                int nj = node_mapping[j];
                
                if (ni != nj) {
                    double dist = distance_matrix_[i][j];
                    double reduced_dist = dist;

                    // If target is in a cycle, reduce weight
                    if (cycle_id[j] != -1) {
                        double cycle_edge_weight = 0;
                        for(const auto& e : edges) {
                            if (e.to == j) {
                                cycle_edge_weight = e.distance;
                                break;
                            }
                        }
                        reduced_dist -= cycle_edge_weight;
                    }

                    if (reduced_dist < new_distances[ni][nj]) {
                        new_distances[ni][nj] = reduced_dist;
                        // We construct Edge manually to avoid needing default constructor if it doesn't have one
                        // But Edge structure is simple struct in mstree.cpp
                        edge_mapping.erase({ni, nj});
                        edge_mapping.insert({{ni, nj}, Edge(i, j, dist)});
                    }
                }
            }
        }
        
        // Recursively solve on contracted graph
        MSTreeV2 contracted_solver(new_distances);
        std::vector<Edge> contracted_edges = contracted_solver.compute();
        
        // Expand solution back to original graph
        std::vector<Edge> final_edges;
        std::set<int> nodes_with_incoming_edges;

        // 1. Add inter-component edges from contracted solution
        for (const auto& e : contracted_edges) {
            if (edge_mapping.count({e.from, e.to})) {
                Edge original = edge_mapping.at({e.from, e.to});
                final_edges.push_back(original);
                nodes_with_incoming_edges.insert(original.to);
            }
        }

        // 2. Add remaining edges for all nodes that don't have incoming edges yet
        // This ensures every node (except root) has exactly one incoming edge
        for (const auto& e : edges) {
            // Skip if this node already has an incoming edge from step 1
            if (nodes_with_incoming_edges.find(e.to) == nodes_with_incoming_edges.end()) {
                // Add this edge - it's the minimum incoming edge for this node
                final_edges.push_back(e);
                nodes_with_incoming_edges.insert(e.to);
            }
        }

        return final_edges;
    }
    
    // Local branch recrafting to improve tree quality
    void recraft_branches(std::vector<Edge>& tree) {
        bool improved = true;
        int max_iterations = 10;
        int iteration = 0;
        
        while (improved && iteration < max_iterations) {
            improved = false;
            iteration++;
            
            // Try swapping adjacent edges
            for (size_t i = 0; i < tree.size(); ++i) {
                for (size_t j = i + 1; j < tree.size(); ++j) {
                    if (can_swap_edges(tree, i, j)) {
                        double current_cost = tree[i].distance + 
                                            tree[j].distance;
                        double swap_cost = calculate_swap_cost(
                            tree, i, j
                        );
                        
                        if (swap_cost < current_cost - 1e-10) {
                            perform_edge_swap(tree, i, j);
                            improved = true;
                        }
                    }
                }
            }
        }
    }
    
    bool can_swap_edges(
        const std::vector<Edge>& tree,
        size_t idx1,
        size_t idx2
    ) {
        const Edge& e1 = tree[idx1];
        const Edge& e2 = tree[idx2];
        
        // Check if edges share a node
        return (e1.from == e2.from || e1.from == e2.to ||
                e1.to == e2.from || e1.to == e2.to);
    }
    
    double calculate_swap_cost(
        const std::vector<Edge>& tree,
        size_t idx1,
        size_t idx2
    ) {
        const Edge& e1 = tree[idx1];
        const Edge& e2 = tree[idx2];
        
        // Try alternative connections
        double cost1 = distance_matrix_[e1.from][e2.to] +
                      distance_matrix_[e2.from][e1.to];
        double cost2 = distance_matrix_[e1.to][e2.from] +
                      distance_matrix_[e2.to][e1.from];
        
        return std::min(cost1, cost2);
    }
    
    void perform_edge_swap(
        std::vector<Edge>& tree,
        size_t idx1,
        size_t idx2
    ) {
        // Swap edge endpoints to minimize cost
        Edge& e1 = tree[idx1];
        Edge& e2 = tree[idx2];
        
        std::swap(e1.to, e2.to);
        e1.distance = distance_matrix_[e1.from][e1.to];
        e2.distance = distance_matrix_[e2.from][e2.to];
    }
};

} // namespace grapetree

#endif // GRAPETREE_MSTREE_V2_H
