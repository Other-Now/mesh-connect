#include "mc/topology.hpp"

#include <algorithm>
#include <cmath>
#include <unordered_map>

namespace mc {

Csr node_to_elems(const Mesh& m) {
    Csr c;
    c.offset.assign(m.nodes() + 1, 0);
    for (size_t e = 0; e < m.elems(); ++e)
        for (int i = 0; i < m.nn(e); ++i) ++c.offset[m.conn[e][i] + 1];
    for (size_t n = 0; n < m.nodes(); ++n) c.offset[n + 1] += c.offset[n];
    c.item.resize(c.offset.back());
    std::vector<int64_t> pos(c.offset.begin(), c.offset.end() - 1);
    for (size_t e = 0; e < m.elems(); ++e)
        for (int i = 0; i < m.nn(e); ++i) c.item[pos[m.conn[e][i]]++] = static_cast<int32_t>(e);
    return c;
}

Components components(const Mesh& m, std::span<const WeldResult> welds,
                      std::span<const std::pair<int32_t, int32_t>> node_pairs) {
    UnionFind uf(m.nodes());
    for (size_t e = 0; e < m.elems(); ++e)
        for (int i = 1; i < m.nn(e); ++i) uf.unite(m.conn[e][0], m.conn[e][i]);
    for (const auto& w : welds)
        for (int s = 1; s < w.nsheets; ++s)
            uf.unite(m.conn[w.sheet[0].elem][0], m.conn[w.sheet[s].elem][0]);
    for (auto [a, b] : node_pairs) uf.unite(a, b);

    Components c;
    std::vector<int32_t> label(m.nodes(), -1);
    c.of_elem.resize(m.elems());
    for (size_t e = 0; e < m.elems(); ++e) {
        int32_t root = uf.find(m.conn[e][0]);
        if (label[root] < 0) {
            label[root] = c.count++;
            c.pids.emplace_back();
            c.elems.push_back(0);
        }
        int32_t k = label[root];
        c.of_elem[e] = k;
        c.elems[k]++;
        c.pids[k].push_back(m.pid[e]);
    }
    for (auto& p : c.pids) {
        std::sort(p.begin(), p.end());
        p.erase(std::unique(p.begin(), p.end()), p.end());
    }
    return c;
}

EdgeStats edge_stats(const Mesh& m) {
    // Sort packed (min,max) node pairs and count runs: no hash map needed.
    std::vector<uint64_t> key;
    key.reserve(m.elems() * 4);
    for (size_t e = 0; e < m.elems(); ++e) {
        int k = m.nn(e);
        for (int i = 0; i < k; ++i) {
            uint32_t a = m.conn[e][i], b = m.conn[e][(i + 1) % k];
            if (a > b) std::swap(a, b);
            key.push_back(uint64_t(a) << 32 | b);
        }
    }
    std::sort(key.begin(), key.end());
    EdgeStats s;
    for (size_t i = 0; i < key.size();) {
        size_t j = i;
        while (j < key.size() && key[j] == key[i]) ++j;
        s.edges++;
        if (j - i == 1) s.free++;
        if (j - i > 2) s.nonmanifold++;
        i = j;
    }
    return s;
}

std::vector<std::pair<int32_t, int32_t>> coincident_nodes(const Mesh& m, double tol) {
    std::vector<std::pair<int32_t, int32_t>> out;
    if (m.nodes() == 0 || !(tol > 0)) return out;
    // Cells 16x the tolerance: a node only needs a neighbour cell on an axis
    // where it lies within tol of the cell wall, so most nodes probe 1 cell
    // instead of 27.
    const double cs = 16 * tol;
    auto cell = [&](double v) { return static_cast<int64_t>(std::floor(v / cs)); };
    auto hash = [](int64_t x, int64_t y, int64_t z) {
        uint64_t h = static_cast<uint64_t>(x) * 0x9E3779B97F4A7C15ull;
        h ^= static_cast<uint64_t>(y) * 0xC2B2AE3D27D4EB4Full + (h << 6) + (h >> 2);
        h ^= static_cast<uint64_t>(z) * 0x165667B19E3779F9ull + (h << 6) + (h >> 2);
        return h;
    };
    std::vector<std::pair<uint64_t, int32_t>> grid(m.nodes());
    for (size_t i = 0; i < m.nodes(); ++i) {
        Vec3 p = m.xyz[i];
        grid[i] = {hash(cell(p.x), cell(p.y), cell(p.z)), static_cast<int32_t>(i)};
    }
    std::sort(grid.begin(), grid.end());
    // cell key -> [begin, end) in grid.
    std::unordered_map<uint64_t, std::pair<uint32_t, uint32_t>> range;
    range.reserve(grid.size());
    for (uint32_t b = 0; b < grid.size();) {
        uint32_t e = b;
        while (e < grid.size() && grid[e].first == grid[b].first) ++e;
        range.emplace(grid[b].first, std::make_pair(b, e));
        b = e;
    }
    const double tol2 = tol * tol;
    uint64_t keys[8];
    for (size_t i = 0; i < m.nodes(); ++i) {
        Vec3 p = m.xyz[i];
        int64_t x0 = cell(p.x - tol), x1 = cell(p.x + tol), y0 = cell(p.y - tol), y1 = cell(p.y + tol),
                z0 = cell(p.z - tol), z1 = cell(p.z + tol);
        int nk = 0;
        for (int64_t cx = x0; cx <= x1; ++cx)
            for (int64_t cy = y0; cy <= y1; ++cy)
                for (int64_t cz = z0; cz <= z1; ++cz) keys[nk++] = hash(cx, cy, cz);
        int nu = 0;  // drop duplicate keys (hash collisions between the <= 8 cells)
        for (int k = 0; k < nk; ++k) {
            bool dup = false;
            for (int u = 0; u < nu; ++u) dup = dup || keys[u] == keys[k];
            if (!dup) keys[nu++] = keys[k];
        }
        nk = nu;
        for (int k = 0; k < nk; ++k) {
            auto it = range.find(keys[k]);
            if (it == range.end()) continue;
            for (uint32_t g = it->second.first; g < it->second.second; ++g) {
                int32_t j = grid[g].second;
                if (j > static_cast<int32_t>(i) && norm2(m.xyz[j] - p) <= tol2)
                    out.emplace_back(static_cast<int32_t>(i), j);
            }
        }
    }
    std::sort(out.begin(), out.end());
    return out;
}

}  // namespace mc
