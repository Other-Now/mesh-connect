#include "mc/ccx.hpp"

#include <algorithm>
#include <cmath>
#include <ostream>
#include <set>
#include <stdexcept>
#include <vector>

#include "mc/connect.hpp"
#include "mc/gen.hpp"
#include "mc/nastran.hpp"

namespace mc {
namespace {

struct Term {
    int64_t node;
    double c;
};

void write_equation(std::ostream& out, const std::vector<Term>& t, int dof) {
    out << "*EQUATION\n" << t.size() << "\n";
    for (size_t i = 0; i < t.size(); ++i) {
        out << t[i].node << ',' << dof << ',' << fmt_real(t[i].c);
        out << ((i % 4 == 3 || i + 1 == t.size()) ? "\n" : ",");
    }
}

}  // namespace

LapInfo write_lap_inp(std::ostream& out, const LapSpec& s) {
    Mesh m;
    int64_t nn = 1, ne = 1;
    add_plate(m, {1, 0, 100, 0, 40, 0.0, s.hA}, 1, nn, ne);
    const size_t nA = m.nodes();
    add_plate(m, {2, 80, 180, 0, 40, 1.0, s.hB}, 2, nn, ne);
    // Two rows of three welds in the 20 mm overlap, deliberately off the nodes.
    std::vector<Vec3> welds;
    for (double x : {85.3, 94.7})
        for (double y : {7.1, 20.3, 32.9}) welds.push_back({x, y, 0.5});

    LapInfo info;
    info.nodes = static_cast<int>(m.nodes());
    info.elems = static_cast<int>(m.elems());
    info.welds = static_cast<int>(welds.size());

    out << "*HEADING\nmesh-connect lap joint, mode=" << s.mode << "\n*NODE\n";
    for (size_t i = 0; i < m.nodes(); ++i)
        out << m.node_id[i] << ',' << fmt_real(m.xyz[i].x) << ',' << fmt_real(m.xyz[i].y) << ','
            << fmt_real(m.xyz[i].z) << '\n';
    for (int pid : {1, 2}) {
        out << "*ELEMENT, TYPE=S4, ELSET=P" << pid << '\n';
        for (size_t e = 0; e < m.elems(); ++e)
            if (m.pid[e] == pid) {
                out << m.elem_id[e];
                for (int k = 0; k < 4; ++k) out << ',' << m.node_id[m.conn[e][k]];
                out << '\n';
            }
    }
    std::vector<size_t> fix, load;
    for (size_t i = 0; i < m.nodes(); ++i) {
        bool inB = i >= nA;
        if (inB && std::fabs(m.xyz[i].x - 180) < 1e-9) fix.push_back(i);
        if (!inB && std::fabs(m.xyz[i].x) < 1e-9) load.push_back(i);
    }
    auto nset = [&](const char* name, const std::vector<size_t>& v) {
        out << "*NSET, NSET=" << name << '\n';
        for (size_t i = 0; i < v.size(); ++i) out << m.node_id[v[i]] << ((i % 8 == 7 || i + 1 == v.size()) ? "\n" : ",");
    };
    nset("FIX", fix);
    nset("LOAD", load);
    out << "*MATERIAL, NAME=STEEL\n*ELASTIC\n210000.,0.3\n";
    for (int pid : {1, 2}) out << "*SHELL SECTION, ELSET=P" << pid << ", MATERIAL=STEEL\n1.\n";

    // Connectors. A node may be the dependent (first) term of only one
    // equation per dof, and must not carry a boundary condition.
    std::set<int32_t> dependent(fix.begin(), fix.end());
    auto pick_dependent = [&](std::vector<std::pair<int32_t, double>> cand) {
        std::sort(cand.begin(), cand.end(), [](auto& a, auto& b) { return std::fabs(a.second) > std::fabs(b.second); });
        for (auto& c : cand)
            if (!dependent.count(c.first)) {
                dependent.insert(c.first);
                return c.first;
            }
        throw std::runtime_error("no free node for a dependent term");
    };
    ConnectOptions opt;
    opt.tol = 1.0;
    opt.max_sheets = 2;
    for (Vec3 w : welds) {
        std::vector<Term> terms;
        if (s.mode == "shape") {
            WeldResult r = connect_one(m, nullptr, w, opt);
            if (r.status != WeldStatus::Ok) throw std::runtime_error("lap weld not realised");
            const SheetHit& A = r.sheet[0].pid == 1 ? r.sheet[0] : r.sheet[1];
            const SheetHit& B = r.sheet[0].pid == 1 ? r.sheet[1] : r.sheet[0];
            std::vector<std::pair<int32_t, double>> bn;
            for (int k = 0; k < m.nn(B.elem); ++k)
                if (std::fabs(B.w[k]) > 1e-12) bn.push_back({m.conn[B.elem][k], B.w[k]});
            int32_t dep = pick_dependent(bn);
            for (auto& [n, c] : bn)
                if (n == dep) terms.insert(terms.begin(), {m.node_id[n], c});
                else terms.push_back({m.node_id[n], c});
            for (int k = 0; k < m.nn(A.elem); ++k)
                if (std::fabs(A.w[k]) > 1e-12) terms.push_back({m.node_id[m.conn[A.elem][k]], -A.w[k]});
        } else if (s.mode == "nearest") {
            int32_t best[3] = {-1, -1, -1};
            double bd[3] = {1e300, 1e300, 1e300};
            for (size_t e = 0; e < m.elems(); ++e)
                for (int k = 0; k < 4; ++k) {
                    int32_t n = m.conn[e][k];
                    Vec3 d = m.xyz[n] - w;
                    double dxy = d.x * d.x + d.y * d.y;
                    int p = m.pid[e];
                    if (dxy < bd[p] && !(p == 2 && dependent.count(n))) bd[p] = dxy, best[p] = n;
                }
            info.max_snap = std::max({info.max_snap, std::sqrt(bd[1]), std::sqrt(bd[2])});
            dependent.insert(best[2]);
            terms = {{m.node_id[best[2]], 1.0}, {m.node_id[best[1]], -1.0}};
        } else if (s.mode == "patch") {
            // Spread over the weld nugget: mean displacement of every node of
            // each sheet within d/2 of the weld (an RBE3-style patch, as
            // ACM-type spot welds do) instead of a single point.
            const double r2 = 0.25 * s.nugget * s.nugget;
            std::vector<std::pair<int32_t, double>> an, bn;
            for (size_t i = 0; i < m.nodes(); ++i) {
                Vec3 d = m.xyz[i] - w;
                if (d.x * d.x + d.y * d.y <= r2) (i < nA ? an : bn).push_back({static_cast<int32_t>(i), 1.0});
            }
            if (an.empty() || bn.empty()) throw std::runtime_error("patch: no node within the nugget radius");
            for (auto& v : an) v.second = 1.0 / static_cast<double>(an.size());
            for (auto& v : bn) v.second = 1.0 / static_cast<double>(bn.size());
            std::sort(bn.begin(), bn.end(), [&](auto& x, auto& y) {
                return norm2(m.xyz[x.first] - w) < norm2(m.xyz[y.first] - w);
            });
            int32_t dep = -1;
            for (auto& v : bn)
                if (!dependent.count(v.first)) {
                    dep = v.first;
                    break;
                }
            if (dep < 0) throw std::runtime_error("patch: no free dependent node");
            dependent.insert(dep);
            terms.push_back({m.node_id[dep], 1.0 / static_cast<double>(bn.size())});
            for (auto& [n, c] : bn)
                if (n != dep) terms.push_back({m.node_id[n], c});
            for (auto& [n, c] : an) terms.push_back({m.node_id[n], -c});
            info.patch_nodes = std::max(info.patch_nodes, static_cast<int>(an.size() + bn.size()));
        } else if (s.mode != "none") {
            throw std::runtime_error("mode must be shape|nearest|patch|none");
        }
        if (!terms.empty()) {
            for (int dof = 1; dof <= 3; ++dof) write_equation(out, terms, dof);
            info.equations += 3;
        }
    }

    out << "*BOUNDARY\nFIX,1,6\n*STEP\n*STATIC\n*CLOAD\n";
    std::sort(load.begin(), load.end(), [&](size_t a, size_t b) { return m.xyz[a].y < m.xyz[b].y; });
    const double f = -s.load / static_cast<double>(load.size() - 1);
    for (size_t i = 0; i < load.size(); ++i)
        out << m.node_id[load[i]] << ",1," << fmt_real((i == 0 || i + 1 == load.size()) ? f / 2 : f) << '\n';
    out << "*NODE PRINT, NSET=FIX, TOTALS=ONLY\nRF\n*NODE PRINT, NSET=LOAD\nU\n*END STEP\n";
    return info;
}

}  // namespace mc
