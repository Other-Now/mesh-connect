#include "mc/bvh.hpp"

#include <algorithm>
#include <stdexcept>

namespace mc {

Bvh::Bvh(const Mesh& m, int leaf_size) {
    const int32_t n = static_cast<int32_t>(m.elems());
    if (n == 0) return;
    std::vector<Aabb> box(n);
    std::vector<Vec3> ctr(n);
    for (int32_t e = 0; e < n; ++e) {
        box[e] = m.box(e);
        ctr[e] = box[e].center();
    }
    order_.resize(n);
    for (int32_t e = 0; e < n; ++e) order_[e] = e;
    nodes_.reserve(2 * static_cast<size_t>(n) / leaf_size + 1);
    nodes_.push_back({});

    // Explicit work list: (node, begin, end, depth). Children are allocated
    // as adjacent pairs so an inner node only stores the left index.
    struct Work {
        int32_t node, b, e, d;
    };
    std::vector<Work> work{{0, 0, n, 1}};
    while (!work.empty()) {
        Work w = work.back();
        work.pop_back();
        depth_ = std::max(depth_, w.d);
        Aabb bb, cb;
        for (int32_t i = w.b; i < w.e; ++i) {
            bb.grow(box[order_[i]]);
            cb.grow(ctr[order_[i]]);
        }
        nodes_[w.node].box = bb;
        if (w.e - w.b <= leaf_size) {
            nodes_[w.node].first = w.b;
            nodes_[w.node].count = w.e - w.b;
            continue;
        }
        int ax = cb.longest_axis();
        int32_t mid = w.b + (w.e - w.b) / 2;
        std::nth_element(order_.begin() + w.b, order_.begin() + mid, order_.begin() + w.e,
                         [&](int32_t a, int32_t b) {
                             double ca = axis(ctr[a], ax), cb2 = axis(ctr[b], ax);
                             return ca < cb2 || (ca == cb2 && a < b);
                         });
        int32_t left = static_cast<int32_t>(nodes_.size());
        nodes_.push_back({});
        nodes_.push_back({});
        nodes_[w.node].first = left;
        nodes_[w.node].count = 0;
        if (w.d + 1 >= 64) throw std::runtime_error("BVH too deep");
        work.push_back({left, w.b, mid, w.d + 1});
        work.push_back({left + 1, mid, w.e, w.d + 1});
    }
    leaf_box_.resize(n);
    for (int32_t i = 0; i < n; ++i) leaf_box_[i] = box[order_[i]];
}

}  // namespace mc
