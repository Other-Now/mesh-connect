#pragma once
// Connectivity analysis: node->element adjacency, connected components
// (before and after connectors), free / non-manifold edges, and coincident
// but unmerged nodes (the classic "looks connected, isn't" mesh defect).
#include <cstdint>
#include <span>
#include <utility>
#include <vector>

#include "mc/connect.hpp"
#include "mc/mesh.hpp"

namespace mc {

struct Csr {
    std::vector<int64_t> offset;  // size nodes+1
    std::vector<int32_t> item;    // element indices
};
Csr node_to_elems(const Mesh& m);

class UnionFind {
public:
    explicit UnionFind(size_t n) : p_(n) {
        for (size_t i = 0; i < n; ++i) p_[i] = static_cast<int32_t>(i);
    }
    int32_t find(int32_t a) {
        while (p_[a] != a) a = p_[a] = p_[p_[a]];  // path halving
        return a;
    }
    bool unite(int32_t a, int32_t b) {
        a = find(a), b = find(b);
        if (a == b) return false;
        if (a > b) std::swap(a, b);  // smaller index wins: deterministic labels
        p_[b] = a;
        return true;
    }

private:
    std::vector<int32_t> p_;
};

struct Components {
    int32_t count = 0;
    std::vector<int32_t> of_elem;                // component label per element, 0..count-1
    std::vector<std::vector<int32_t>> pids;      // sorted distinct part ids per component
    std::vector<int64_t> elems;                  // element count per component
};
// Components of the element graph (elements sharing a node are connected),
// optionally plus connector links and extra node pairs (e.g. merges).
Components components(const Mesh& m, std::span<const WeldResult> welds = {},
                      std::span<const std::pair<int32_t, int32_t>> node_pairs = {});

struct EdgeStats {
    int64_t edges = 0;        // distinct edges
    int64_t free = 0;         // used by exactly one element (part boundary, holes, cracks)
    int64_t nonmanifold = 0;  // used by 3+ elements (T-joints)
};
EdgeStats edge_stats(const Mesh& m);

// Pairs of distinct nodes closer than tol (spatial hash, 27-cell search).
std::vector<std::pair<int32_t, int32_t>> coincident_nodes(const Mesh& m, double tol);

}  // namespace mc
