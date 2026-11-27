#include "mc/geom.hpp"

#include <algorithm>

namespace mc {
namespace {

constexpr double kInsideEps = 1e-10;  // parametric tolerance for "inside"

// Closest point on segment ab; t in [0,1].
inline double seg_t(Vec3 p, Vec3 a, Vec3 b) {
    Vec3 ab = b - a;
    double l2 = norm2(ab);
    if (l2 <= 0) return 0;
    return std::clamp(dot(p - a, ab) / l2, 0.0, 1.0);
}

inline void take(Projection& best, const Projection& c) {
    if (c.dist2 < best.dist2) best = c;
}

}  // namespace

std::array<double, 4> quad_shape(double xi, double eta) {
    return {0.25 * (1 - xi) * (1 - eta), 0.25 * (1 + xi) * (1 - eta),
            0.25 * (1 + xi) * (1 + eta), 0.25 * (1 - xi) * (1 + eta)};
}

Vec3 quad_point(const std::array<Vec3, 4>& x, double xi, double eta) {
    auto n = quad_shape(xi, eta);
    return n[0] * x[0] + n[1] * x[1] + n[2] * x[2] + n[3] * x[3];
}

Projection project_tri(Vec3 p, Vec3 a, Vec3 b, Vec3 c) {
    Projection best;
    // Unconstrained projection onto the plane, in barycentric coordinates.
    Vec3 e0 = b - a, e1 = c - a, d = p - a;
    double a00 = dot(e0, e0), a01 = dot(e0, e1), a11 = dot(e1, e1);
    double det = a00 * a11 - a01 * a01;
    if (det > 1e-300 * (a00 * a11)) {
        double b0 = dot(d, e0), b1 = dot(d, e1);
        double v = (a11 * b0 - a01 * b1) / det;
        double w = (a00 * b1 - a01 * b0) / det;
        double u = 1 - v - w;
        if (u >= -kInsideEps && v >= -kInsideEps && w >= -kInsideEps) {
            Projection r;
            r.w = {u, v, w, 0};
            r.point = u * a + v * b + w * c;
            r.dist2 = norm2(p - r.point);
            r.normal = true;
            return r;
        }
    }
    // Outside the face region (or degenerate): closest point is on an edge.
    const Vec3 v[3] = {a, b, c};
    for (int i = 0; i < 3; ++i) {
        int j = (i + 1) % 3;
        double t = seg_t(p, v[i], v[j]);
        Projection r;
        r.w[i] = 1 - t;
        r.w[j] = t;
        r.point = (1 - t) * v[i] + t * v[j];
        r.dist2 = norm2(p - r.point);
        take(best, r);
    }
    return best;
}

Projection project_quad(Vec3 p_world, const std::array<Vec3, 4>& xw) {
    // Work relative to the first corner. In model coordinates (~1000 mm) a
    // double resolves ~1e-13 mm, so the iteration must not ask for more than
    // that: it once used an absolute 1e-14 step tolerance, never "converged"
    // away from the origin, and fell back to the edges for 38% of welds.
    const Vec3 o = xw[0];
    const std::array<Vec3, 4> x = {Vec3{}, xw[1] - o, xw[2] - o, xw[3] - o};
    const Vec3 p = p_world - o;
    // Minimise f(xi,eta) = |x(xi,eta) - p|^2. Far from the answer, use
    // Gauss-Newton with a backtracking line search (f must decrease); close to
    // it, add the x_xieta curvature term (full Newton) for fast convergence on
    // warped quads. Plain full Newton from the centre overshot on distorted
    // quads and cycled (2,2) -> (2,-2) -> (0.09,-0.11) -> ... forever: 2 of
    // 100,000 welds landed on an edge instead of inside their element.
    const Vec3 xxe = 0.25 * (x[0] - x[1] + x[2] - x[3]);
    auto f = [&](double a, double b) { return norm2(quad_point(x, a, b) - p); };
    double xi = 0, eta = 0, fc = f(0, 0), last_step = 1e300;
    bool converged = false;
    for (int it = 0; it < 60; ++it) {
        Vec3 xx = 0.25 * ((1 - eta) * (x[1] - x[0]) + (1 + eta) * (x[2] - x[3]));
        Vec3 xe = 0.25 * ((1 - xi) * (x[3] - x[0]) + (1 + xi) * (x[2] - x[1]));
        Vec3 r = quad_point(x, xi, eta) - p;
        double g0 = dot(r, xx), g1 = dot(r, xe);
        double h00 = dot(xx, xx), h11 = dot(xe, xe), h01 = dot(xx, xe);
        if (last_step < 1e-2) {  // near the solution: full Newton if it is a descent step
            double hn = h01 + dot(r, xxe);
            if (h00 * h11 - hn * hn > 0) h01 = hn;
        }
        double det = h00 * h11 - h01 * h01;
        if (!(det > 0)) break;  // degenerate element
        double dxi = -(h11 * g0 - h01 * g1) / det;
        double deta = -(h00 * g1 - h01 * g0) / det;
        double t = 1;
        double nx = xi + dxi, ne = eta + deta, fn = f(nx, ne);
        // Backtrack only far from the answer: near it, f changes by less than
        // its rounding error and the search would stall the iteration early.
        for (int k = 0; k < 40 && last_step >= 1e-2 && fn > fc; ++k) {
            t *= 0.5;
            nx = xi + t * dxi, ne = eta + t * deta, fn = f(nx, ne);
        }
        last_step = std::max(std::fabs(t * dxi), std::fabs(t * deta));
        xi = nx, eta = ne;
        fc = std::min(fc, fn);
        if (last_step < 1e-12) {  // parametric, relative to element
            converged = true;
            break;
        }
        if (std::fabs(xi) > 4 || std::fabs(eta) > 4) break;  // heading far outside: use edges
    }
    if (converged && std::fabs(xi) <= 1 + kInsideEps && std::fabs(eta) <= 1 + kInsideEps) {
        xi = std::clamp(xi, -1.0, 1.0);
        eta = std::clamp(eta, -1.0, 1.0);
        Projection r;
        r.w = quad_shape(xi, eta);
        r.point = quad_point(xw, xi, eta);
        r.dist2 = norm2(p_world - r.point);
        r.normal = true;
        return r;
    }
    // Edges of a bilinear quad are straight, so edge projection is exact.
    Projection best;
    for (int i = 0; i < 4; ++i) {
        int j = (i + 1) % 4;
        double t = seg_t(p_world, xw[i], xw[j]);
        Projection r;
        r.w[i] = 1 - t;
        r.w[j] = t;
        r.point = (1 - t) * xw[i] + t * xw[j];
        r.dist2 = norm2(p_world - r.point);
        take(best, r);
    }
    return best;
}

Projection project_quad_split(Vec3 p, const std::array<Vec3, 4>& x) {
    Projection t0 = project_tri(p, x[0], x[1], x[2]);
    Projection t1 = project_tri(p, x[0], x[2], x[3]);
    Projection r;
    if (t0.dist2 <= t1.dist2) {
        r = t0;
        r.w = {t0.w[0], t0.w[1], t0.w[2], 0};
    } else {
        r = t1;
        r.w = {t1.w[0], 0, t1.w[1], t1.w[2]};
    }
    return r;
}

}  // namespace mc
