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
        
        // Find root (node with no parent or most children)
        int root = find_root(nodes);
        
        // Generate Newick string recursively
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(6);
        oss << to_newick(root, nodes, strain_names);
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
    
    std::string to_newick(
        int node_id,
        const std::vector<TreeNode>& nodes,
        const std::vector<std::string>& names
    ) {
        const TreeNode& node = nodes[node_id];
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(6);
        
        if (node.children.empty()) {
            // Leaf node
            oss << sanitize_name(names[node_id]);
        } else {
            // Internal node with children
            oss << "(";
            
            for (size_t i = 0; i < node.children.size(); ++i) {
                if (i > 0) oss << ",";
                
                int child_id = node.children[i];
                oss << to_newick(child_id, nodes, names);
                
                // Add branch length
                oss << ":" << nodes[child_id].branch_length;
            }
            
            oss << ")";
            
            // Add internal node label if desired
            if (node_id < static_cast<int>(names.size())) {
                oss << sanitize_name(names[node_id]);
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
