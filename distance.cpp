// distance.cpp - Distance matrix computation for GrapeTree
// Supports both symmetric and asymmetric distance calculations

#ifndef GRAPETREE_DISTANCE_H
#define GRAPETREE_DISTANCE_H

#include <vector>
#include <string>
#include <cmath>
#include <limits>

namespace grapetree {

class DistanceMatrix {
public:
    struct ProfileData {
        std::vector<std::string> strain_names;
        std::vector<std::vector<int>> profiles;
        int n_strains;
        int n_genes;
    };
    
    enum MissingHandler {
        IGNORE = 0,        // Ignore missing in pairwise comparison
        REMOVE_COLUMN = 1, // Remove columns with missing
        TREAT_AS_ALLELE = 2, // Treat missing as unique allele
        ABSOLUTE_DIFF = 3  // Count as absolute difference
    };
    
private:
    ProfileData data_;
    
public:
    explicit DistanceMatrix(const ProfileData& data) : data_(data) {}
    
    // Compute symmetric distance matrix (for MSTree and NJ)
    std::vector<std::vector<double>> compute_symmetric(
        MissingHandler handler = IGNORE
    ) {
        std::vector<std::vector<double>> matrix(
            data_.n_strains,
            std::vector<double>(data_.n_strains, 0.0)
        );
        
        for (int i = 0; i < data_.n_strains; ++i) {
            for (int j = i + 1; j < data_.n_strains; ++j) {
                double dist = compute_pairwise_distance(
                    data_.profiles[i],
                    data_.profiles[j],
                    handler
                );
                matrix[i][j] = dist;
                matrix[j][i] = dist;
            }
        }
        
        return matrix;
    }
    
    // Compute asymmetric distance matrix (for MSTreeV2)
    std::vector<std::vector<double>> compute_asymmetric() {
        std::vector<std::vector<double>> matrix(
            data_.n_strains,
            std::vector<double>(data_.n_strains, 0.0)
        );
        
        for (int i = 0; i < data_.n_strains; ++i) {
            for (int j = 0; j < data_.n_strains; ++j) {
                if (i == j) {
                    matrix[i][j] = 0.0;
                } else {
                    matrix[i][j] = compute_directional_distance(
                        data_.profiles[i],
                        data_.profiles[j]
                    );
                }
            }
        }
        
        return matrix;
    }
    
    // Compute p-distance for aligned sequences
    std::vector<std::vector<double>> compute_p_distance(
        const std::vector<std::string>& sequences
    ) {
        int n = sequences.size();
        std::vector<std::vector<double>> matrix(
            n, 
            std::vector<double>(n, 0.0)
        );
        
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                double dist = p_distance(sequences[i], sequences[j]);
                matrix[i][j] = dist;
                matrix[j][i] = dist;
            }
        }
        
        return matrix;
    }
    
private:
    // Pairwise allelic distance
    double compute_pairwise_distance(
        const std::vector<int>& profile1,
        const std::vector<int>& profile2,
        MissingHandler handler
    ) {
        int differences = 0;
        int valid_positions = 0;
        
        for (int k = 0; k < data_.n_genes; ++k) {
            int allele1 = profile1[k];
            int allele2 = profile2[k];
            
            // Missing data represented as 0 or negative values
            bool missing1 = (allele1 <= 0);
            bool missing2 = (allele2 <= 0);
            
            if (missing1 || missing2) {
                switch (handler) {
                    case IGNORE:
                        // Skip this position
                        continue;
                        
                    case REMOVE_COLUMN:
                        // Skip this position (same as IGNORE)
                        continue;
                        
                    case TREAT_AS_ALLELE:
                        // Missing is treated as a unique allele
                        if (missing1 != missing2 || 
                            (!missing1 && allele1 != allele2)) {
                            differences++;
                        }
                        valid_positions++;
                        break;
                        
                    case ABSOLUTE_DIFF:
                        // Count missing as difference
                        differences++;
                        valid_positions++;
                        break;
                }
            } else {
                // Both alleles present
                if (allele1 != allele2) {
                    differences++;
                }
                valid_positions++;
            }
        }
        
        return static_cast<double>(differences);
    }
    
    // Directional distance for MSTreeV2 (asymmetric)
    double compute_directional_distance(
        const std::vector<int>& from_profile,
        const std::vector<int>& to_profile
    ) {
        int differences = 0;
        int missing_in_from = 0;
        
        for (int k = 0; k < data_.n_genes; ++k) {
            int from_allele = from_profile[k];
            int to_allele = to_profile[k];
            
            if (from_allele <= 0) {
                // Missing data in 'from' profile
                missing_in_from++;
            } else if (to_allele > 0 && from_allele != to_allele) {
                // Both present but different
                differences++;
            }
        }
        
        // Asymmetric: penalize missing data in source node
        // This encourages the tree to grow from complete profiles
        return static_cast<double>(differences) + 
               0.5 * static_cast<double>(missing_in_from);
    }
    
    // p-distance for DNA sequences
    double p_distance(
        const std::string& seq1,
        const std::string& seq2
    ) {
        if (seq1.length() != seq2.length()) {
            return std::numeric_limits<double>::max();
        }
        
        int differences = 0;
        int valid_positions = 0;
        
        for (size_t i = 0; i < seq1.length(); ++i) {
            char c1 = std::toupper(seq1[i]);
            char c2 = std::toupper(seq2[i]);
            
            // Skip gaps and ambiguous bases
            if (c1 == '-' || c1 == 'N' || c2 == '-' || c2 == 'N') {
                continue;
            }
            
            if (c1 != c2) {
                differences++;
            }
            valid_positions++;
        }
        
        if (valid_positions == 0) {
            return 0.0;
        }
        
        return static_cast<double>(differences) / 
               static_cast<double>(valid_positions);
    }
};

} // namespace grapetree

#endif // GRAPETREE_DISTANCE_H
