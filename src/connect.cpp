#include "mc/connect.hpp"

#include <algorithm>
#include <cmath>

namespace mc {
namespace {

struct Best {
    int32_t pid;
    int32_t elem;
    Projection proj;
};

// Strict total order: distance, then element index. BVH and brute force visit
// elements in different orders; this makes them pick the same element on ties.
inline bool better(const Projection& a, int32_t ea, const Projection& b, int32_t eb) {
    return a.dist2 < b.dist2 || (a.dist2 == b.dist2 && ea < eb);
}

}  // namespace

WeldResult connect_one(const Mesh& m, const Bvh* bvh, Vec3 p, const ConnectOptions& o) {
    const double tol2 = o.tol * o.tol;
    Best best[16];  // closest element per part; a weld rarely sees more than 3-4 parts
    int nb = 0;
    auto visit = [&](int32_t e, const Aabb& box) {
        int32_t pid = m.pid[e];
        if (o.prune) {
            // The element's box is a lower bound on its distance. Skip the
            // (much dearer) projection if the box is already farther than the
            // tolerance or than this part's best hit so far. Strict '>', so a
            // tie that the element-index rule might break is still projected.
            double lb = box.dist2(p);
            if (lb > tol2) return;
            for (int i = 0; i < nb; ++i)
                if (best[i].pid == pid && lb > best[i].proj.dist2) return;
        }
        Projection pr = project(m, e, p);
        if (pr.dist2 > tol2) return;
        for (int i = 0; i < nb; ++i) {
            if (best[i].pid == pid) {
                if (better(pr, e, best[i].proj, best[i].elem)) best[i] = {pid, e, pr};
                return;
            }
        }
        if (nb < 16) best[nb++] = {pid, e, pr};
    };
    if (bvh) {
        bvh->query(p, o.tol, visit);
    } else {
        for (int32_t e = 0; e < static_cast<int32_t>(m.elems()); ++e) visit(e, m.box(e));
    }
    std::sort(best, best + nb, [](const Best& a, const Best& b) {
        return a.proj.dist2 < b.proj.dist2 || (a.proj.dist2 == b.proj.dist2 && a.pid < b.pid);
    });

    WeldResult r;
    r.nsheets = std::min(nb, std::min(o.max_sheets, 3));
    for (int i = 0; i < r.nsheets; ++i) {
        SheetHit& h = r.sheet[i];
        h.elem = best[i].elem;
        h.pid = best[i].pid;
        h.dist = std::sqrt(best[i].proj.dist2);
        h.point = best[i].proj.point;
        h.w = best[i].proj.w;
        h.normal = best[i].proj.normal;
    }
    // Order the stack through the thickness (along the closest sheet's normal),
    // so adjacent entries are physically adjacent sheets: a 3-sheet weld then
    // joins bottom-middle and middle-top, never bottom-top.
    if (r.nsheets == 3) {
        auto x = m.corners(r.sheet[0].elem);
        Vec3 n = m.nn(r.sheet[0].elem) == 3 ? cross(x[1] - x[0], x[2] - x[0])
                                            : cross(x[2] - x[0], x[3] - x[1]);
        auto s = [&](const SheetHit& h) { return dot(h.point - p, n); };
        std::sort(r.sheet.begin(), r.sheet.begin() + 3, [&](const SheetHit& a, const SheetHit& b) {
            double sa = s(a), sb = s(b);
            return sa < sb || (sa == sb && a.pid < b.pid);
        });
    }
    r.status = r.nsheets >= 2 ? WeldStatus::Ok : (r.nsheets == 1 ? WeldStatus::OneSheet : WeldStatus::NoSheet);
    return r;
}

namespace {

// 21 bits per axis interleaved: nearby points get nearby codes.
uint64_t morton(Vec3 p, const Aabb& b) {
    auto q = [](double v, double lo, double hi) {
        double t = hi > lo ? (v - lo) / (hi - lo) : 0;
        return static_cast<uint64_t>(std::clamp(t, 0.0, 1.0) * 2097151.0);
    };
    auto spread = [](uint64_t x) {
        x &= 0x1fffff;
        x = (x | x << 32) & 0x1f00000000ffffull;
        x = (x | x << 16) & 0x1f0000ff0000ffull;
        x = (x | x << 8) & 0x100f00f00f00f00full;
        x = (x | x << 4) & 0x10c30c30c30c30c3ull;
        x = (x | x << 2) & 0x1249249249249249ull;
        return x;
    };
    return spread(q(p.x, b.lo.x, b.hi.x)) | spread(q(p.y, b.lo.y, b.hi.y)) << 1 | spread(q(p.z, b.lo.z, b.hi.z)) << 2;
}

}  // namespace

std::vector<WeldResult> connect_all(const Mesh& m, const Bvh* bvh, std::span<const Vec3> pts,
                                    const ConnectOptions& o) {
    std::vector<WeldResult> out(pts.size());
    const int64_t n = static_cast<int64_t>(pts.size());
    // Visit welds in Z-curve order so consecutive queries touch the same BVH
    // nodes and elements. Results still go to each weld's own slot, so the
    // output order (and content) does not change.
    std::vector<int64_t> order(pts.size());
    for (int64_t i = 0; i < n; ++i) order[i] = i;
    if (o.spatial_order && n > 1) {
        Aabb b;
        for (Vec3 p : pts) b.grow(p);
        std::vector<std::pair<uint64_t, int64_t>> key(pts.size());
        for (int64_t i = 0; i < n; ++i) key[i] = {morton(pts[i], b), i};
        std::sort(key.begin(), key.end());
        for (int64_t i = 0; i < n; ++i) order[i] = key[i].second;
    }
#pragma omp parallel for schedule(dynamic, 256)
    for (int64_t k = 0; k < n; ++k) {
        int64_t i = order[k];
        out[i] = connect_one(m, bvh, pts[i], o);
    }
    return out;
}

}  // namespace mc
