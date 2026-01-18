// mstree.cpp - Classical Minimum Spanning Tree (Prim's algorithm)
// Implements MSTree with eBurst tiebreaking heuristic

#ifndef GRAPETREE_MSTREE_H
#define GRAPETREE_MSTREE_H

#include <vector>
#include <string>
#include <limits>
#include <algorithm>
#include <map>

namespace grapetree {

struct Edge {
    int from;
    int to;
    double distance;
    
    Edge(int f, int t, double d) : from(f), to(t), distance(d) {}
};

class MSTree {
public:
    enum Heuristic {
        EBURST,
        HARMONIC
    };
    
private:
    int n_nodes_;
    std::vector<std::vector<double>> distance_matrix_;
    Heuristic heuristic_;
    
public:
    MSTree(
        const std::vector<std::vector<double>>& distances,
        Heuristic heuristic = EBURST
    ) : distance_matrix_(distances),
        heuristic_(heuristic),
        n_nodes_(distances.size()) {}
    
    std::vector<Edge> compute() {
        std::vector<Edge> tree_edges;
        tree_edges.reserve(n_nodes_ - 1);
        
        std::vector<bool> in_tree(n_nodes_, false);
        std::vector<double> min_distance(n_nodes_, 
                                         std::numeric_limits<double>::max());
        std::vector<int> parent(n_nodes_, -1);
        
        // Start with node 0 (arbitrary choice)
        int start_node = 0;
        in_tree[start_node] = true;
        min_distance[start_node] = 0.0;
        
        // Initialize distances from start node
        for (int i = 0; i < n_nodes_; ++i) {
            if (i != start_node) {
                min_distance[i] = distance_matrix_[start_node][i];
                parent[i] = start_node;
            }
        }
        
        // Build tree: add n-1 edges
        for (int count = 1; count < n_nodes_; ++count) {
            // Find minimum distance node not yet in tree
            double min_dist = std::numeric_limits<double>::max();
            
            for (int i = 0; i < n_nodes_; ++i) {
                if (!in_tree[i] && min_distance[i] < min_dist) {
                    min_dist = min_distance[i];
                }
            }
            
            // Apply tiebreaking heuristic
            int min_node = select_node_with_tiebreak(
                min_distance,
                in_tree,
                min_dist
            );
            
            // Add edge to tree
            in_tree[min_node] = true;
            tree_edges.emplace_back(
                parent[min_node],
                min_node,
                min_dist
            );
            
            // Update distances to remaining nodes
            for (int i = 0; i < n_nodes_; ++i) {
                if (!in_tree[i]) {
                    double new_dist = distance_matrix_[min_node][i];
                    if (new_dist < min_distance[i]) {
                        min_distance[i] = new_dist;
                        parent[i] = min_node;
                    }
                }
            }
        }
        
        return tree_edges;
    }
    
private:
    int select_node_with_tiebreak(
        const std::vector<double>& distances,
        const std::vector<bool>& in_tree,
        double min_dist
    ) {
        // Collect all nodes at minimum distance
        std::vector<int> candidates;
        for (int i = 0; i < n_nodes_; ++i) {
            if (!in_tree[i] && 
                std::abs(distances[i] - min_dist) < 1e-10) {
                candidates.push_back(i);
            }
        }
        
        if (candidates.size() == 1) {
            return candidates[0];
        }
        
        // Apply heuristic for tiebreaking
        if (heuristic_ == EBURST) {
            return apply_eburst_tiebreak(candidates, in_tree, min_dist);
        } else {
            return apply_harmonic_tiebreak(candidates);
        }
    }
    
    // eBurst: select node with most connections at min_dist
    int apply_eburst_tiebreak(
        const std::vector<int>& candidates,
        const std::vector<bool>& in_tree,
        double min_dist
    ) {
        int best_node = candidates[0];
        int max_connections = 0;
        
        for (int node : candidates) {
            int connections = 0;
            
            // Count connections to nodes already in tree
            for (int j = 0; j < n_nodes_; ++j) {
                if (in_tree[j] && 
                    std::abs(distance_matrix_[node][j] - min_dist) < 1e-10) {
                    connections++;
                }
            }
            
            if (connections > max_connections) {
                max_connections = connections;
                best_node = node;
            } else if (connections == max_connections) {
                // Secondary tiebreak: prefer lower node index
                if (node < best_node) {
                    best_node = node;
                }
            }
        }
        
        return best_node;
    }
    
    // Harmonic mean: prefer nodes with smaller average distance
    int apply_harmonic_tiebreak(
        const std::vector<int>& candidates
    ) {
        int best_node = candidates[0];
        double best_score = -1.0;
        
        for (int node : candidates) {
            double score = compute_harmonic_mean_score(node);
            
            if (score > best_score) {
                best_score = score;
                best_node = node;
            } else if (std::abs(score - best_score) < 1e-10) {
                // Tiebreak by node index
                if (node < best_node) {
                    best_node = node;
                }
            }
        }
        
        return best_node;
    }
    
    double compute_harmonic_mean_score(int node) {
        double sum_reciprocals = 0.0;
        int count = 0;
        
        for (int i = 0; i < n_nodes_; ++i) {
            if (i == node) continue;
            
            double dist = distance_matrix_[node][i];
            if (dist > 0.0) {
                sum_reciprocals += 1.0 / dist;
                count++;
            }
        }
        
        if (count == 0) {
            return 0.0;
        }
        
        // Harmonic mean = n / sum(1/x_i)
        return static_cast<double>(count) / sum_reciprocals;
    }
};

} // namespace grapetree

#endif // GRAPETREE_MSTREE_H
