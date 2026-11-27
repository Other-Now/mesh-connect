#include "mc/gen.hpp"

#include <algorithm>
#include <cmath>

namespace mc {

void add_plate(Mesh& m, const Plate& p, uint64_t seed, int64_t& next_node, int64_t& next_elem) {
    const int nx = std::max(1, static_cast<int>(std::lround((p.x1 - p.x0) / p.h)));
    const int ny = std::max(1, static_cast<int>(std::lround((p.y1 - p.y0) / p.h)));
    const double dx = (p.x1 - p.x0) / nx, dy = (p.y1 - p.y0) / ny;
    Rng rng(seed);
    std::vector<int32_t> idx(static_cast<size_t>(nx + 1) * (ny + 1));
    for (int j = 0; j <= ny; ++j)
        for (int i = 0; i <= nx; ++i) {
            double x = p.x0 + i * dx, y = p.y0 + j * dy;
            if (i > 0 && i < nx && j > 0 && j < ny && p.jitter > 0) {
                x += rng.uni(-p.jitter, p.jitter) * dx;
                y += rng.uni(-p.jitter, p.jitter) * dy;
            }
            idx[static_cast<size_t>(j) * (nx + 1) + i] = m.add_node(next_node++, {x, y, p.z});
        }
    auto at = [&](int i, int j) { return idx[static_cast<size_t>(j) * (nx + 1) + i]; };
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
            int32_t a = at(i, j), b = at(i + 1, j), c = at(i + 1, j + 1), d = at(i, j + 1);
            if (!p.tri) {
                m.add_elem(next_elem++, p.pid, a, b, c, d);
            } else if ((i + j) % 2 == 0) {
                m.add_elem(next_elem++, p.pid, a, b, c);
                m.add_elem(next_elem++, p.pid, a, c, d);
            } else {
                m.add_elem(next_elem++, p.pid, a, b, d);
                m.add_elem(next_elem++, p.pid, b, c, d);
            }
        }
    m.thickness[p.pid] = p.t;
}

// Layout (side view, x to the right):
//
//   part 1 (z=0)  ==========
//   part 2 (z=g)          ==========          neighbours overlap by `flange`
//   part 3 (z=0)                  ==========
//
// Even flanges also get a reinforcement patch on top (z=2g): 3-sheet stacks.
// Odd flanges get a stiffener 2.5 mm below the weld plane: close, but outside
// the search tolerance, so it must NOT be welded.
Assembly gen_assembly(const AssemblySpec& s) {
    Assembly a;
    const double L = 1000, W = 500, flange = 40, g = s.gap;
    const double hs[3] = {1.0, 1.15, 1.3};
    int64_t nn = 1, ne = 1;
    for (int k = 0; k < s.parts; ++k) {
        double x0 = k * (L - flange);
        a.plates.push_back({k + 1, x0, x0 + L, 0, W, (k % 2) * g, s.h * hs[k % 3], k % 2 == 1, 0.25});
    }
    for (int k = 0; k + 1 < s.parts; ++k) {
        double x0 = (k + 1) * (L - flange), x1 = x0 + flange;
        if (k % 2 == 0)
            a.plates.push_back({100 + k, x0, x1, 100, 200, 2 * g, s.h, false, 0.25});  // patch
        else
            a.plates.push_back({200 + k, x0, x1, 300, 400, g / 2 - 2.5, s.h, true, 0.25});  // stiffener
    }
    for (size_t i = 0; i < a.plates.size(); ++i) add_plate(a.mesh, a.plates[i], s.seed * 1000 + i, nn, ne);

    Rng rng(s.seed);
    const double margin = 3;
    for (int w = 0; w < s.welds; ++w) {
        double u = rng.uni();
        Label l;
        if (u < 0.1) {  // single sheet: weld point off the flanges
            int k = static_cast<int>(rng.next() % s.parts);
            double x0 = k * (L - flange);
            l.kind = "1sheet";
            double z = (k % 2) * g;
            l.p = {rng.uni(x0 + flange + margin, x0 + L - flange - margin), rng.uni(margin, W - margin), z + g / 2};
            l.pids = {k + 1};
            l.zs = {z};
        } else {
            int k = static_cast<int>(rng.next() % (s.parts - 1));
            double x0 = (k + 1) * (L - flange);
            double x = rng.uni(x0 + margin, x0 + flange - margin);
            bool three = (k % 2 == 0) && u < 0.2;
            if (three) {
                l.kind = "3sheet";
                l.p = {x, rng.uni(100 + margin, 200 - margin), g};
                l.pids = {k + 1, k + 2, 100 + k};
                l.zs = {(k % 2) * g, ((k + 1) % 2) * g, 2 * g};
            } else {
                double y;
                do y = rng.uni(margin, W - margin);
                while (k % 2 == 0 && y > 100 - margin && y < 200 + margin);  // patch area is 3-sheet
                l.kind = "2sheet";
                l.p = {x, y, g / 2};
                l.pids = {k + 1, k + 2};
                l.zs = {(k % 2) * g, ((k + 1) % 2) * g};
            }
        }
        a.labels.push_back(std::move(l));
    }
    return a;
}

}  // namespace mc
