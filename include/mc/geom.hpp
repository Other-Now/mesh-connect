#pragma once
// Small 3D geometry kit: vectors, boxes, and point-to-shell-element projection.
#include <array>
#include <cmath>
#include <limits>

namespace mc {

struct Vec3 {
    double x = 0, y = 0, z = 0;
};
inline Vec3 operator+(Vec3 a, Vec3 b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
inline Vec3 operator-(Vec3 a, Vec3 b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
inline Vec3 operator*(double s, Vec3 a) { return {s * a.x, s * a.y, s * a.z}; }
inline double dot(Vec3 a, Vec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline Vec3 cross(Vec3 a, Vec3 b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline double norm2(Vec3 a) { return dot(a, a); }
inline double norm(Vec3 a) { return std::sqrt(dot(a, a)); }

struct Aabb {
    static constexpr double inf = std::numeric_limits<double>::infinity();
    Vec3 lo{inf, inf, inf}, hi{-inf, -inf, -inf};
    void grow(Vec3 p) {
        lo = {std::fmin(lo.x, p.x), std::fmin(lo.y, p.y), std::fmin(lo.z, p.z)};
        hi = {std::fmax(hi.x, p.x), std::fmax(hi.y, p.y), std::fmax(hi.z, p.z)};
    }
    void grow(const Aabb& b) { grow(b.lo), grow(b.hi); }
    Vec3 center() const { return 0.5 * (lo + hi); }
    // Squared distance from p to the box (0 inside).
    double dist2(Vec3 p) const {
        double dx = std::fmax(std::fmax(lo.x - p.x, 0.0), p.x - hi.x);
        double dy = std::fmax(std::fmax(lo.y - p.y, 0.0), p.y - hi.y);
        double dz = std::fmax(std::fmax(lo.z - p.z, 0.0), p.z - hi.z);
        return dx * dx + dy * dy + dz * dz;
    }
    int longest_axis() const {
        Vec3 d = hi - lo;
        return d.x >= d.y && d.x >= d.z ? 0 : (d.y >= d.z ? 1 : 2);
    }
};
inline double axis(Vec3 v, int a) { return a == 0 ? v.x : (a == 1 ? v.y : v.z); }

// Closest point on one shell element to a query point.
//  w      : interpolation weights of the element's nodes at that point
//           (barycentric for TRIA3, bilinear shape functions for QUAD4).
//           A connector uses them to tie the weld to the element's nodes.
//  normal : true if the point projects along the surface normal into the
//           element interior; false if the closest point was clamped to an edge.
struct Projection {
    Vec3 point;
    double dist2 = std::numeric_limits<double>::infinity();
    std::array<double, 4> w{0, 0, 0, 0};
    bool normal = false;
};

Projection project_tri(Vec3 p, Vec3 a, Vec3 b, Vec3 c);
// Bilinear QUAD4 (nodes in order (-1,-1),(1,-1),(1,1),(-1,1)), possibly warped.
// Newton iteration on (xi, eta), falling back to the straight edges.
Projection project_quad(Vec3 p, const std::array<Vec3, 4>& x);
// Shortcut that splits the quad into two triangles (what a naive tool does).
// Exact for flat parallelograms only; used to measure the error on warped quads.
Projection project_quad_split(Vec3 p, const std::array<Vec3, 4>& x);

std::array<double, 4> quad_shape(double xi, double eta);
Vec3 quad_point(const std::array<Vec3, 4>& x, double xi, double eta);

}  // namespace mc
