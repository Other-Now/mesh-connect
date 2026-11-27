#pragma once
// Bounding-volume hierarchy over element boxes, for "which elements are
// within r of this point" queries. Built once (median split on the longest
// axis), then read-only, so many threads can query it at once.
#include <cstdint>
#include <vector>

#include "mc/mesh.hpp"

namespace mc {

class Bvh {
public:
    explicit Bvh(const Mesh& m, int leaf_size = 4);

    // Calls f(element_index, element_box) for every element whose box is
    // within r of p. Boxes are stored in leaf order, next to each other, so a
    // caller can reject an element without touching the mesh's node arrays.
    template <class F>
    void query(Vec3 p, double r, F&& f) const {
        if (nodes_.empty()) return;
        const double r2 = r * r;
        int32_t stack[64];
        int sp = 0;
        stack[sp++] = 0;
        while (sp) {
            const Node& n = nodes_[stack[--sp]];
            if (n.box.dist2(p) > r2) continue;
            if (n.count) {
                for (int32_t i = 0; i < n.count; ++i)
                    if (leaf_box_[n.first + i].dist2(p) <= r2) f(order_[n.first + i], leaf_box_[n.first + i]);
            } else {
                stack[sp++] = n.first;      // left child
                stack[sp++] = n.first + 1;  // right child
            }
        }
    }

    size_t node_count() const { return nodes_.size(); }
    int depth() const { return depth_; }

private:
    struct Node {
        Aabb box;
        int32_t first = 0;  // leaf: start in order_; inner: index of left child
        int32_t count = 0;  // > 0 for leaves
    };
    std::vector<Node> nodes_;
    std::vector<int32_t> order_;
    std::vector<Aabb> leaf_box_;  // element boxes, same order as order_
    int depth_ = 0;
};

}  // namespace mc
