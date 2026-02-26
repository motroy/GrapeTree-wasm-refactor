// wasm_interface.cpp - Emscripten bindings for GrapeTree
// Exposes C++ functionality to JavaScript

#include <emscripten/bind.h>
#include <emscripten/val.h>
#include <emscripten/emscripten.h>
#include <nlohmann/json.hpp>
#include <string>
#include <vector>
#include <sstream>
#include <functional>

// Include our GrapeTree modules
#include "distance.cpp"
#include "mstree.cpp"
#include "mstree_v2.cpp"
#include "newick.cpp"

using namespace emscripten;
using json = nlohmann::json;
using namespace grapetree;

// Helper to report progress back to JavaScript worker
EM_JS(void, report_progress_js, (const char* task, double percent), {
    if (typeof postMessage === 'function') {
        postMessage({
            type: 'progress',
            task: UTF8ToString(task),
            value: percent
        });
    }
});

void report_progress(const char* task, double percent) {
    report_progress_js(task, percent);
}

// Helper function to parse JSON profile data
DistanceMatrix::ProfileData parse_profile_json(const std::string& json_str) {
    json data = json::parse(json_str);
    
    DistanceMatrix::ProfileData profile;
    profile.strain_names = data["strains"].get<std::vector<std::string>>();
    profile.profiles = data["profiles"].get<std::vector<std::vector<int>>>();
    profile.n_strains = profile.strain_names.size();
    profile.n_genes = profile.profiles.empty() ? 0 : profile.profiles[0].size();
    
    return profile;
}

// Convert edges to JSON
json edges_to_json(
    const std::vector<Edge>& edges,
    const std::vector<std::string>& strain_names
) {
    json result = json::array();
    
    for (const Edge& e : edges) {
        json edge_obj;
        edge_obj["from"] = e.from;
        edge_obj["to"] = e.to;
        edge_obj["from_name"] = strain_names[e.from];
        edge_obj["to_name"] = strain_names[e.to];
        edge_obj["distance"] = e.distance;
        result.push_back(edge_obj);
    }
    
    return result;
}

// Main tree computation function
std::string compute_tree(
    const std::string& profile_json,
    const std::string& method,
    const std::string& matrix_type,
    int missing_handler,
    const std::string& heuristic
) {
    try {
        report_progress("Parsing input", 0.0);

        // Parse input
        auto profile_data = parse_profile_json(profile_json);
        
        // Compute distance matrix
        DistanceMatrix dm(profile_data);

        // Distance calculation takes roughly 90% of time for large N
        auto dist_progress = [](double p) {
            report_progress("Computing distance matrix", p * 0.9);
        };

        // Tree construction takes remaining 10%
        auto tree_progress = [](double p) {
            report_progress("Building tree", 90.0 + p * 0.1);
        };

        // Compute tree
        std::vector<Edge> tree_edges;

        if (method == "MSTree") {
            // MSTree holds a const-ref to the 2-D matrix, so we keep it alive
            // for the duration of the call.  No extra copy is created.
            auto distances = dm.compute_symmetric(
                static_cast<DistanceMatrix::MissingHandler>(missing_handler),
                dist_progress
            );
            MSTree::Heuristic h = (heuristic == "harmonic") ?
                MSTree::HARMONIC : MSTree::EBURST;
            MSTree mst(distances, h);
            tree_edges = mst.compute();
            tree_progress(100.0);
        } else if (method == "MSTreeV2") {
            // Use the flat-matrix path: compute one contiguous 1-D vector and
            // move it directly into MSTreeV2.  This avoids having both the 2-D
            // matrix AND a flat copy alive simultaneously, cutting peak distance
            // memory roughly in half (e.g. ~85 MB instead of ~170 MB for
            // N=3 256 strains).
            int n = profile_data.n_strains;
            std::vector<double> flat_distances;
            if (matrix_type == "symmetric") {
                flat_distances = dm.compute_symmetric_flat(
                    static_cast<DistanceMatrix::MissingHandler>(missing_handler),
                    dist_progress
                );
                MSTreeV2 mst2(std::move(flat_distances), n, profile_data.n_genes);
                tree_edges = mst2.compute(tree_progress);
            } else {
                // Compute both the normalized [0,1] ratios (for the algorithm)
                // and the raw allele-difference counts (for branch lengths).
                flat_distances = dm.compute_asymmetric_flat(dist_progress);
                std::vector<double> flat_raw = dm.compute_asymmetric_raw_flat();
                MSTreeV2 mst2(std::move(flat_distances), std::move(flat_raw),
                              n, profile_data.n_genes);
                tree_edges = mst2.compute(tree_progress);
            }
        } else {
            throw std::runtime_error("Unknown method: " + method);
        }
        
        report_progress("Formatting output", 100.0);

        // Format output as Newick
        NewickFormatter formatter;
        std::string newick = formatter.format(
            tree_edges,
            profile_data.strain_names
        );
        
        // Build JSON response
        json response;
        response["success"] = true;
        response["newick"] = newick;
        response["edges"] = edges_to_json(tree_edges, profile_data.strain_names);
        response["n_nodes"] = profile_data.n_strains;
        response["n_edges"] = tree_edges.size();
        
        return response.dump();
        
    } catch (const std::exception& e) {
        json error_response;
        error_response["success"] = false;
        error_response["error"] = e.what();
        return error_response.dump();
    }
}

// Distance matrix computation
std::string compute_distance_matrix(
    const std::string& profile_json,
    const std::string& matrix_type,
    int missing_handler
) {
    try {
        report_progress("Parsing input", 0.0);

        auto profile_data = parse_profile_json(profile_json);
        
        DistanceMatrix dm(profile_data);
        std::vector<std::vector<double>> distances;
        
        auto dist_progress = [](double p) {
            report_progress("Computing distance matrix", p);
        };

        if (matrix_type == "symmetric") {
            distances = dm.compute_symmetric(
                static_cast<DistanceMatrix::MissingHandler>(missing_handler),
                dist_progress
            );
        } else {
            distances = dm.compute_asymmetric(dist_progress);
        }
        
        report_progress("Formatting output", 100.0);

        // Convert to JSON
        json response;
        response["success"] = true;
        response["matrix"] = distances;
        response["strain_names"] = profile_data.strain_names;
        response["n_strains"] = profile_data.n_strains;
        
        return response.dump();
        
    } catch (const std::exception& e) {
        json error_response;
        error_response["success"] = false;
        error_response["error"] = e.what();
        return error_response.dump();
    }
}

// Emscripten bindings
EMSCRIPTEN_BINDINGS(grapetree_module) {
    function("compute_tree", &compute_tree);
    function("compute_distance_matrix", &compute_distance_matrix);
    
    // Also expose individual components if needed
    enum_<DistanceMatrix::MissingHandler>("MissingHandler")
        .value("IGNORE", DistanceMatrix::IGNORE)
        .value("REMOVE_COLUMN", DistanceMatrix::REMOVE_COLUMN)
        .value("TREAT_AS_ALLELE", DistanceMatrix::TREAT_AS_ALLELE)
        .value("ABSOLUTE_DIFF", DistanceMatrix::ABSOLUTE_DIFF);
    
    enum_<MSTree::Heuristic>("Heuristic")
        .value("EBURST", MSTree::EBURST)
        .value("HARMONIC", MSTree::HARMONIC);
}

// For standalone compilation (non-WASM)
#ifndef __EMSCRIPTEN__
int main() {
    // Test harness
    std::string test_data = R"({
        "strains": ["A", "B", "C"],
        "profiles": [
            [1, 2, 3],
            [1, 2, 4],
            [1, 3, 3]
        ]
    })";
    
    std::string result = compute_tree(
        test_data,
        "MSTreeV2",
        "asymmetric",
        0,
        "harmonic"
    );
    
    std::cout << "Result: " << result << std::endl;
    
    return 0;
}
#endif
