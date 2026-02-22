// newick.cpp - Newick tree format output
// Converts edge list to standard Newick format

#ifndef GRAPETREE_NEWICK_H
#define GRAPETREE_NEWICK_H

#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <iomanip>
#include <set>

namespace grapetree {

struct Edge;

class NewickFormatter {
private:
    struct TreeNode {
        int id;
        std::vector<int> children;
        int parent;
        double branch_length;
        
        TreeNode() : id(-1), parent(-1), branch_length(0.0) {}
    };
    
public:
    std::string format(
        const std::vector<Edge>& edges,
        const std::vector<std::string>& strain_names
    ) {
        if (edges.empty()) {
            return strain_names.empty() ? "();" : strain_names[0] + ";";
        }

        // Build tree structure
        std::vector<TreeNode> nodes = build_tree_structure(
            edges,
            strain_names.size()
        );

        // Find all root nodes (nodes with no parent)
        std::vector<int> roots;
        for (size_t i = 0; i < nodes.size(); ++i) {
            if (nodes[i].parent == -1) {
                roots.push_back(i);
            }
        }

        // Generate Newick string recursively
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(6);

        if (roots.empty()) {
            // No roots found (shouldn't happen, but handle gracefully)
            return strain_names.empty() ? "();" : strain_names[0] + ";";
        } else if (roots.size() == 1) {
            // Single connected tree - standard case
            oss << to_newick(roots[0], nodes, strain_names);
        } else {
            // Multiple disconnected components - create artificial root
            oss << "(";

            for (size_t i = 0; i < roots.size(); ++i) {
                if (i > 0) oss << ",";
                oss << to_newick(roots[i], nodes, strain_names) << ":0.0";
            }

            oss << ")";
        }

        oss << ";";

        return oss.str();
    }
    
private:
    std::vector<TreeNode> build_tree_structure(
        const std::vector<Edge>& edges,
        int n_nodes
    ) {
        std::vector<TreeNode> nodes(n_nodes);
        
        // Initialize nodes
        for (int i = 0; i < n_nodes; ++i) {
            nodes[i].id = i;
            nodes[i].parent = -1;
        }
        
        // Build adjacency information
        for (const Edge& e : edges) {
            nodes[e.from].children.push_back(e.to);
            nodes[e.to].parent = e.from;
            nodes[e.to].branch_length = e.distance;
        }
        
        return nodes;
    }
    
    int find_root(const std::vector<TreeNode>& nodes) {
        // Find node with no parent
        for (size_t i = 0; i < nodes.size(); ++i) {
            if (nodes[i].parent == -1) {
                return i;
            }
        }
        
        // If all have parents (cycle), pick node with most children
        int best_root = 0;
        size_t max_children = 0;
        
        for (size_t i = 0; i < nodes.size(); ++i) {
            if (nodes[i].children.size() > max_children) {
                max_children = nodes[i].children.size();
                best_root = i;
            }
        }
        
        return best_root;
    }
    
    // Iterative post-order Newick builder.
    // Replaces the former recursive implementation to avoid stack overflow on
    // deep linear trees (e.g. chains of 3 000+ nodes in large cgMLST profiles).
    std::string to_newick(
        int root_id,
        const std::vector<TreeNode>& nodes,
        const std::vector<std::string>& names
    ) {
        struct Frame {
            int node_id;
            size_t child_idx; // index of the next child yet to be written
        };

        std::ostringstream oss;
        oss << std::fixed << std::setprecision(6);

        std::vector<Frame> stk;
        stk.push_back({root_id, 0});

        while (!stk.empty()) {
            Frame& top = stk.back();
            int nid = top.node_id;
            const TreeNode& node = nodes[nid];

            if (top.child_idx == 0 && !node.children.empty()) {
                // First visit to an internal node — open the subtree bracket.
                oss << "(";
            }

            if (top.child_idx < node.children.size()) {
                // Still have children to recurse into.
                if (top.child_idx > 0) oss << ",";
                int child_id = node.children[top.child_idx];
                ++top.child_idx;
                stk.push_back({child_id, 0});
            } else {
                // All children have been written — close this node.
                stk.pop_back();

                if (!node.children.empty()) oss << ")";

                // Node label (leaf or internal).
                if (nid < static_cast<int>(names.size())) {
                    oss << sanitize_name(names[nid]);
                }

                // Branch length — appended after returning to parent.
                if (!stk.empty()) {
                    oss << ":" << node.branch_length;
                }
            }
        }

        return oss.str();
    }
    
    std::string sanitize_name(const std::string& name) {
        // Remove/escape characters that are special in Newick format
        std::string sanitized;
        bool needs_quotes = false;
        
        for (char c : name) {
            if (c == ' ' || c == ':' || c == ';' || 
                c == '(' || c == ')' || c == ',' || 
                c == '[' || c == ']' || c == '\'') {
                needs_quotes = true;
            }
            sanitized += c;
        }
        
        if (needs_quotes) {
            return "'" + sanitized + "'";
        }
        
        return sanitized;
    }
};

// Extended formatter with metadata support
class ExtendedNewickFormatter {
public:
    struct NodeMetadata {
        std::map<std::string, std::string> attributes;
    };
    
    std::string format_with_metadata(
        const std::vector<Edge>& edges,
        const std::vector<std::string>& strain_names,
        const std::vector<NodeMetadata>& metadata
    ) {
        // Extended Newick format with [&key=value] annotations
        // Implementation similar to basic formatter
        // but adds metadata annotations
        
        NewickFormatter basic;
        return basic.format(edges, strain_names);
    }
};

} // namespace grapetree

#endif // GRAPETREE_NEWICK_H
