#pragma once
// Shell mesh in structure-of-arrays form. Nodes and elements keep their
// solver ids (Nastran GRID / element ids) but everything internal uses
// dense 0-based indices.
#include <array>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "mc/geom.hpp"

namespace mc {

struct Mesh {
    // nodes
    std::vector<int64_t> node_id;
    std::vector<Vec3> xyz;
    std::unordered_map<int64_t, int32_t> node_index;
    // elements: 3 (TRIA3) or 4 (QUAD4) nodes; conn[e][3] = -1 for triangles
    std::vector<int64_t> elem_id;
    std::vector<int32_t> pid;
    std::vector<std::array<int32_t, 4>> conn;
    // shell thickness per property id (PSHELL), if known
    std::unordered_map<int32_t, double> thickness;

    size_t nodes() const { return xyz.size(); }
    size_t elems() const { return conn.size(); }
    int nn(size_t e) const { return conn[e][3] < 0 ? 3 : 4; }

    int32_t add_node(int64_t id, Vec3 p) {
        auto [it, fresh] = node_index.emplace(id, static_cast<int32_t>(xyz.size()));
        if (!fresh) throw std::runtime_error("duplicate GRID id " + std::to_string(id));
        node_id.push_back(id);
        xyz.push_back(p);
        return it->second;
    }
    // node indices (not ids); n3 = -1 for a triangle
    void add_elem(int64_t id, int32_t p, int32_t n0, int32_t n1, int32_t n2, int32_t n3 = -1) {
        elem_id.push_back(id);
        pid.push_back(p);
        conn.push_back({n0, n1, n2, n3});
    }
    std::array<Vec3, 4> corners(size_t e) const {
        const auto& c = conn[e];
        return {xyz[c[0]], xyz[c[1]], xyz[c[2]], c[3] < 0 ? xyz[c[2]] : xyz[c[3]]};
    }
    Aabb box(size_t e) const {
        Aabb b;
        for (int i = 0; i < nn(e); ++i) b.grow(xyz[conn[e][i]]);
        return b;
    }
};

// Project p onto element e (dispatches on TRIA3 / QUAD4).
inline Projection project(const Mesh& m, size_t e, Vec3 p) {
    auto x = m.corners(e);
    return m.nn(e) == 3 ? project_tri(p, x[0], x[1], x[2]) : project_quad(p, x);
}

}  // namespace mc
