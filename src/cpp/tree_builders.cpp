// tree_builders.cpp - tree-building pipelines shared by the WebAssembly
// interface and the native tests.  Included after distance.cpp, mstree.cpp
// and mstree_v2.cpp.

#ifndef GRAPETREE_TREE_BUILDERS_H
#define GRAPETREE_TREE_BUILDERS_H

#include <functional>
#include <map>
#include <utility>
#include <vector>

namespace grapetree {

// MSTreeV2 with GrapeTree's asymmetric allelic distance.
//
// Follows the original pipeline:
//   * strains with identical profiles are merged into one node; the tree is
//     built on the unique profiles and each duplicate is then attached to its
//     representative by a zero-length edge;
//   * the tree is built from the directed distance
//         d(u->v) = #{l : pi_l(v) != 0 and pi_l(u) != pi_l(v)} * n_loci / N_v
//     on the allele-count scale (n_loci / N_v rescales for missing data);
//   * nodes within a couple of differences of a more central node are attached
//     to it directly (shortcuts);
//   * the best root is chosen automatically;
//   * each branch is then labelled with the symmetric count of loci that
//     differ where both strains have data.
//
// Returned edges refer to the strain indices of `data`.
inline std::vector<Edge> build_mstree_v2_asymmetric(
    const DistanceMatrix::ProfileData& data,
    std::function<void(double)> dist_progress,
    std::function<void(double)> tree_progress)
{
    const int n_all = data.n_strains;
    const int n_loci = data.n_genes;

    // Group identical profiles; `reps` holds the first strain of each group.
    std::map<std::vector<int>, int> group_of_profile;
    std::vector<int> reps;
    std::vector<int> group(n_all);
    std::vector<int> group_sizes;
    for (int i = 0; i < n_all; ++i) {
        auto it = group_of_profile.find(data.profiles[i]);
        if (it == group_of_profile.end()) {
            it = group_of_profile.emplace(data.profiles[i], static_cast<int>(reps.size())).first;
            reps.push_back(i);
            group_sizes.push_back(0);
        }
        group[i] = it->second;
        ++group_sizes[group[i]];
    }

    std::vector<Edge> edges;
    const int n_unique = static_cast<int>(reps.size());

    if (n_unique >= 2) {
        DistanceMatrix::ProfileData unique;
        unique.n_strains = n_unique;
        unique.n_genes = n_loci;
        unique.strain_names.reserve(n_unique);
        unique.profiles.reserve(n_unique);
        for (int r : reps) {
            unique.strain_names.push_back(data.strain_names[r]);
            unique.profiles.push_back(data.profiles[r]);
        }

        DistanceMatrix dm(unique);

        std::vector<double> distances = dm.compute_asymmetric_flat(dist_progress);
        for (double& d : distances) d *= n_loci;

        std::vector<double> branch_lengths = dm.compute_symmetric_link_flat();

        MSTreeV2 mst(std::move(distances), std::move(branch_lengths), n_unique, n_loci);
        mst.set_distances_in_counts(true);
        mst.set_auto_root(true);

        for (Edge e : mst.compute_grapetree(tree_progress, group_sizes)) {
            edges.emplace_back(reps[e.from], reps[e.to], e.distance);
        }
    }

    // Duplicates hang off their representative at distance zero.
    for (int i = 0; i < n_all; ++i) {
        if (reps[group[i]] != i) edges.emplace_back(reps[group[i]], i, 0.0);
    }
    return edges;
}

} // namespace grapetree

#endif // GRAPETREE_TREE_BUILDERS_H
