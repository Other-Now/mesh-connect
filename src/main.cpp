// mesh-connect CLI. Every subcommand prints one JSON object on stdout so the
// scripts in scripts/ can collect results without parsing prose.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "mc/bvh.hpp"
#include "mc/ccx.hpp"
#include "mc/connect.hpp"
#include "mc/gen.hpp"
#include "mc/nastran.hpp"
#include "mc/topology.hpp"

using namespace mc;
namespace fs = std::filesystem;

namespace {

using Clock = std::chrono::steady_clock;
double secs(Clock::time_point t0) { return std::chrono::duration<double>(Clock::now() - t0).count(); }

struct Args {
    std::vector<std::string> pos;
    std::map<std::string, std::string> opt;
    double num(const std::string& k, double d) const { return opt.count(k) ? std::stod(opt.at(k)) : d; }
    std::string str(const std::string& k, const std::string& d) const { return opt.count(k) ? opt.at(k) : d; }
};
Args parse(int argc, char** argv) {
    Args a;
    for (int i = 2; i < argc; ++i) {
        std::string s = argv[i];
        if (s.rfind("--", 0) == 0 && i + 1 < argc) a.opt[s.substr(2)] = argv[++i];
        else a.pos.push_back(s);
    }
    return a;
}

void set_threads(int n) {
#ifdef _OPENMP
    if (n > 0) omp_set_num_threads(n);
#else
    (void)n;
#endif
}
int hw_threads() {
#ifdef _OPENMP
    return omp_get_num_procs();
#else
    return 1;
#endif
}
int max_threads() {
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

std::vector<Vec3> read_points(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("cannot open " + path);
    std::vector<Vec3> v;
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || !(std::isdigit(static_cast<unsigned char>(line[0])) || line[0] == '-' || line[0] == '.')) continue;
        Vec3 p;
        if (std::sscanf(line.c_str(), "%lf,%lf,%lf", &p.x, &p.y, &p.z) == 3) v.push_back(p);
    }
    return v;
}

std::vector<int> split_ints(const std::string& s, char sep) {
    std::vector<int> v;
    std::stringstream ss(s);
    std::string t;
    while (std::getline(ss, t, sep))
        if (!t.empty()) v.push_back(std::stoi(t));
    return v;
}

std::vector<Label> read_labels(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("cannot open " + path);
    std::vector<Label> v;
    std::string line;
    std::getline(f, line);  // header
    while (std::getline(f, line)) {
        std::stringstream ss(line);
        std::string id, kind, x, y, z, pids, zs;
        std::getline(ss, id, ','), std::getline(ss, kind, ','), std::getline(ss, x, ','), std::getline(ss, y, ','),
            std::getline(ss, z, ','), std::getline(ss, pids, ','), std::getline(ss, zs, ',');
        Label l;
        l.kind = kind;
        l.p = {std::stod(x), std::stod(y), std::stod(z)};
        l.pids = split_ints(pids, ';');
        std::stringstream zz(zs);
        std::string t;
        while (std::getline(zz, t, ';')) l.zs.push_back(std::stod(t));
        v.push_back(std::move(l));
    }
    return v;
}

// FNV-1a over everything a connector depends on: thread-count independence
// is checked by comparing these.
uint64_t hash_results(const std::vector<WeldResult>& r) {
    uint64_t h = 1469598103934665603ull;
    auto mix = [&](const void* p, size_t n) {
        auto b = static_cast<const unsigned char*>(p);
        for (size_t i = 0; i < n; ++i) h = (h ^ b[i]) * 1099511628211ull;
    };
    for (const auto& w : r) {
        mix(&w.status, sizeof w.status);
        mix(&w.nsheets, sizeof w.nsheets);
        for (int s = 0; s < w.nsheets; ++s) {
            mix(&w.sheet[s].elem, sizeof(int32_t));
            mix(&w.sheet[s].point, sizeof(Vec3));
            mix(w.sheet[s].w.data(), sizeof(double) * 4);
        }
    }
    return h;
}

bool same(const WeldResult& a, const WeldResult& b) {
    if (a.status != b.status || a.nsheets != b.nsheets) return false;
    for (int s = 0; s < a.nsheets; ++s)
        if (a.sheet[s].elem != b.sheet[s].elem || std::memcmp(&a.sheet[s].point, &b.sheet[s].point, sizeof(Vec3)) ||
            a.sheet[s].w != b.sheet[s].w)
            return false;
    return true;
}

std::string json_pids(const std::vector<int32_t>& p) {
    std::string s = "[";
    for (size_t i = 0; i < p.size(); ++i) s += (i ? "," : "") + std::to_string(p[i]);
    return s + "]";
}

double median(std::vector<double> v) {
    std::sort(v.begin(), v.end());
    return v.empty() ? 0 : v[v.size() / 2];
}

// ---------------------------------------------------------------- commands

int cmd_gen(const Args& a) {
    if (a.pos.empty()) throw std::runtime_error("gen <out_dir>");
    AssemblySpec s;
    s.parts = static_cast<int>(a.num("parts", 8));
    s.h = a.num("h", 2.0);
    s.welds = static_cast<int>(a.num("welds", 20000));
    s.seed = static_cast<uint64_t>(a.num("seed", 1));
    auto t0 = Clock::now();
    Assembly as = gen_assembly(s);
    fs::create_directories(a.pos[0]);
    fs::path d = a.pos[0];
    {
        std::ofstream f(d / "mesh.bdf");
        write_nastran(f, as.mesh);
    }
    std::ofstream w(d / "welds.csv"), l(d / "labels.csv");
    w << "x,y,z\n";
    l << "id,kind,x,y,z,pids,zs\n";
    for (size_t i = 0; i < as.labels.size(); ++i) {
        const Label& lb = as.labels[i];
        std::string x = fmt_real(lb.p.x), y = fmt_real(lb.p.y), z = fmt_real(lb.p.z);
        w << x << ',' << y << ',' << z << '\n';
        l << i << ',' << lb.kind << ',' << x << ',' << y << ',' << z << ',';
        for (size_t k = 0; k < lb.pids.size(); ++k) l << (k ? ";" : "") << lb.pids[k];
        l << ',';
        for (size_t k = 0; k < lb.zs.size(); ++k) l << (k ? ";" : "") << fmt_real(lb.zs[k]);
        l << '\n';
    }
    std::printf("{\"cmd\":\"gen\",\"nodes\":%zu,\"elems\":%zu,\"parts\":%zu,\"welds\":%zu,\"secs\":%.3f}\n",
                as.mesh.nodes(), as.mesh.elems(), as.plates.size(), as.labels.size(), secs(t0));
    return 0;
}

int cmd_check(const Args& a) {
    if (a.pos.empty()) throw std::runtime_error("check <mesh.bdf>");
    auto t0 = Clock::now();
    NastranStats st;
    Mesh m = read_nastran_file(a.pos[0], &st);
    double t_read = secs(t0);
    t0 = Clock::now();
    Csr csr = node_to_elems(m);
    double t_csr = secs(t0);
    t0 = Clock::now();
    Components c = components(m);
    double t_comp = secs(t0);
    t0 = Clock::now();
    EdgeStats es = edge_stats(m);
    double t_edge = secs(t0);
    double tol = a.num("merge-tol", 1e-3);
    t0 = Clock::now();
    auto pairs = coincident_nodes(m, tol);
    double t_coin = secs(t0);
    Components merged = components(m, {}, pairs);
    int64_t orphans = 0;
    for (size_t n = 0; n < m.nodes(); ++n) orphans += csr.offset[n + 1] == csr.offset[n];
    std::printf("{\"cmd\":\"check\",\"nodes\":%zu,\"elems\":%zu,\"quads\":%zu,\"trias\":%zu,\"skipped_cards\":%zu,"
                "\"components\":%d,\"edges\":%lld,\"free_edges\":%lld,\"nonmanifold_edges\":%lld,\"orphan_nodes\":%lld,"
                "\"coincident_pairs\":%zu,\"merge_tol\":%g,\"components_if_merged\":%d,"
                "\"secs\":{\"read\":%.3f,\"csr\":%.3f,\"components\":%.3f,\"edges\":%.3f,\"coincident\":%.3f}}\n",
                m.nodes(), m.elems(), st.quads, st.trias, [&] { size_t s = 0; for (auto& [k, v] : st.skipped) s += v; return s; }(),
                c.count, static_cast<long long>(es.edges), static_cast<long long>(es.free),
                static_cast<long long>(es.nonmanifold), static_cast<long long>(orphans), pairs.size(), tol, merged.count,
                t_read, t_csr, t_comp, t_edge, t_coin);
    return 0;
}

int cmd_connect(const Args& a) {
    if (a.pos.size() < 2) throw std::runtime_error("connect <mesh.bdf> <welds.csv> [--tol --max-sheets --threads --out --report]");
    set_threads(static_cast<int>(a.num("threads", 0)));
    auto t0 = Clock::now();
    Mesh m = read_nastran_file(a.pos[0]);
    auto pts = read_points(a.pos[1]);
    double t_read = secs(t0);
    t0 = Clock::now();
    Bvh bvh(m);
    double t_bvh = secs(t0);
    ConnectOptions o;
    o.tol = a.num("tol", 1.2);
    o.max_sheets = static_cast<int>(a.num("max-sheets", 3));
    t0 = Clock::now();
    auto res = connect_all(m, &bvh, pts, o);
    double t_conn = secs(t0);

    size_t ok = 0, one = 0, none = 0, three = 0, off_edge = 0;
    for (const auto& r : res) {
        ok += r.status == WeldStatus::Ok;
        one += r.status == WeldStatus::OneSheet;
        none += r.status == WeldStatus::NoSheet;
        three += r.nsheets == 3;
        for (int s = 0; s < r.nsheets; ++s) off_edge += !r.sheet[s].normal;
    }
    std::vector<WeldResult> good;
    for (const auto& r : res)
        if (r.status == WeldStatus::Ok) good.push_back(r);
    Components before = components(m), after = components(m, good);
    // Everything not in the largest component after welding is "floating".
    int32_t main_c = static_cast<int32_t>(std::max_element(after.elems.begin(), after.elems.end()) - after.elems.begin());
    std::vector<int32_t> floating;
    for (int32_t k = 0; k < after.count; ++k)
        if (k != main_c) floating.insert(floating.end(), after.pids[k].begin(), after.pids[k].end());
    std::sort(floating.begin(), floating.end());

    if (a.opt.count("out")) {
        int64_t first = 1;
        for (auto id : m.node_id) first = std::max(first, id + 1);
        for (auto id : m.elem_id) first = std::max(first, id + 1);
        std::ofstream f(a.opt.at("out"));
        write_cwelds(f, m, pts, res, first, static_cast<int32_t>(a.num("pwid", 9000)));
    }
    if (a.opt.count("report")) {
        std::ofstream f(a.opt.at("report"));
        f << "weld,status,nsheets,pid_a,elem_a,dist_a,pid_b,elem_b,dist_b,pid_c,elem_c,dist_c\n";
        for (size_t i = 0; i < res.size(); ++i) {
            f << i << ',' << static_cast<int>(res[i].status) << ',' << res[i].nsheets;
            for (int s = 0; s < 3; ++s) {
                if (s < res[i].nsheets)
                    f << ',' << res[i].sheet[s].pid << ',' << m.elem_id[res[i].sheet[s].elem] << ','
                      << fmt_real(res[i].sheet[s].dist);
                else
                    f << ",,,";
            }
            f << '\n';
        }
    }
    std::printf("{\"cmd\":\"connect\",\"welds\":%zu,\"ok\":%zu,\"three_sheet\":%zu,\"one_sheet\":%zu,\"no_sheet\":%zu,"
                "\"off_edge_hits\":%zu,\"components_before\":%d,\"components_after\":%d,\"floating_pids\":%s,"
                "\"threads\":%d,\"secs\":{\"read\":%.3f,\"bvh\":%.3f,\"connect\":%.3f}}\n",
                res.size(), ok, three, one, none, off_edge, before.count, after.count, json_pids(floating).c_str(),
                max_threads(), t_read, t_bvh, t_conn);
    return 0;
}

int cmd_eval(const Args& a) {
    if (a.pos.size() < 2) throw std::runtime_error("eval <mesh.bdf> <labels.csv>");
    Mesh m = read_nastran_file(a.pos[0]);
    auto labels = read_labels(a.pos[1]);
    Bvh bvh(m);
    ConnectOptions o;
    o.tol = a.num("tol", 1.2);
    std::vector<Vec3> pts;
    for (auto& l : labels) pts.push_back(l.p);
    auto res = connect_all(m, &bvh, pts, o);

    std::map<std::string, std::pair<size_t, size_t>> by_kind;  // kind -> (correct, total)
    double max_err = 0;
    size_t centroid_checked = 0, centroid_wrong = 0;
    std::vector<size_t> wrong;
    for (size_t i = 0; i < labels.size(); ++i) {
        const Label& l = labels[i];
        const WeldResult& r = res[i];
        bool want_ok = l.kind != "1sheet";
        bool good = (r.status == WeldStatus::Ok) == want_ok && r.nsheets == static_cast<int>(l.pids.size());
        std::vector<int32_t> got;
        for (int s = 0; s < r.nsheets; ++s) got.push_back(r.sheet[s].pid);
        std::sort(got.begin(), got.end());
        good = good && got == l.pids;
        for (int s = 0; good && s < r.nsheets; ++s) {
            size_t k = std::find(l.pids.begin(), l.pids.end(), r.sheet[s].pid) - l.pids.begin();
            Vec3 truth{l.p.x, l.p.y, l.zs[k]};
            double err = norm(r.sheet[s].point - truth);
            max_err = std::max(max_err, err);
            good = good && err < 1e-9 && r.sheet[s].normal;
        }
        if (good && r.nsheets == 3) {  // middle of the stack must be the middle sheet
            std::vector<double> z = l.zs;
            std::sort(z.begin(), z.end());
            good = std::fabs(r.sheet[1].point.z - z[1]) < 1e-9;
        }
        auto& bk = by_kind[l.kind];
        bk.first += good;
        bk.second++;
        if (!good && wrong.size() < 10) wrong.push_back(i);

        // Baseline: "nearest element centroid" per part. Does that element
        // actually contain the true projection?
        if (want_ok)
            for (size_t k = 0; k < l.pids.size(); ++k) {
                Vec3 truth{l.p.x, l.p.y, l.zs[k]};
                int32_t best = -1;
                double bd = 1e300;
                bvh.query(truth, 10.0, [&](int32_t e, const Aabb&) {
                    if (m.pid[e] != l.pids[k]) return;
                    Vec3 c{};
                    for (int j = 0; j < m.nn(e); ++j) c = c + m.xyz[m.conn[e][j]];
                    c = (1.0 / m.nn(e)) * c;
                    double d = norm2(c - truth);
                    if (d < bd || (d == bd && e < best)) bd = d, best = e;
                });
                Projection pr = project(m, best, truth);
                centroid_checked++;
                centroid_wrong += !(pr.normal && pr.dist2 < 1e-18);
            }
    }
    size_t correct = 0, total = 0;
    std::string kinds;
    for (auto& [k, v] : by_kind) {
        correct += v.first, total += v.second;
        kinds += (kinds.empty() ? "" : ",") + ("\"" + k + "\":[" + std::to_string(v.first) + "," + std::to_string(v.second) + "]");
    }
    std::string wl;
    for (size_t i : wrong) wl += (wl.empty() ? "" : ",") + std::to_string(i);
    std::printf("{\"cmd\":\"eval\",\"correct\":%zu,\"total\":%zu,\"by_kind\":{%s},\"max_projection_error\":%.3g,"
                "\"centroid_baseline\":{\"sheet_hits\":%zu,\"wrong_element\":%zu},\"first_wrong\":[%s]}\n",
                correct, total, kinds.c_str(), max_err, centroid_checked, centroid_wrong, wl.c_str());
    return correct == total ? 0 : 1;
}

int cmd_bench(const Args& a) {
    if (a.pos.size() < 2) throw std::runtime_error("bench <mesh.bdf> <welds.csv>");
    auto t0 = Clock::now();
    Mesh m = read_nastran_file(a.pos[0]);
    double t_read = secs(t0);
    double mb = static_cast<double>(fs::file_size(a.pos[0])) / 1e6;
    auto pts = read_points(a.pos[1]);
    t0 = Clock::now();
    Bvh bvh(m);
    double t_bvh = secs(t0);
    ConnectOptions o;
    o.tol = a.num("tol", 1.2);
    int reps = static_cast<int>(a.num("reps", 5));

    // Thread scaling: median of reps, plus a hash to prove identical output.
    std::string scaling;
    uint64_t h0 = 0;
    bool all_same = true;
    double t1 = 0;
    std::vector<WeldResult> ref;
    for (int th : split_ints(a.str("threads", "1,2,3,4"), ',')) {
        set_threads(th);
        std::vector<double> ts;
        uint64_t h = 0;
        for (int r = 0; r < reps; ++r) {
            t0 = Clock::now();
            auto res = connect_all(m, &bvh, pts, o);
            ts.push_back(secs(t0));
            h = hash_results(res);
            if (ref.empty()) ref = std::move(res);
        }
        if (!h0) h0 = h;
        all_same = all_same && h == h0;
        double t = median(ts);
        if (th == 1) t1 = t;
        char buf[200];
        std::snprintf(buf, sizeof buf, "%s{\"threads\":%d,\"secs\":%.4f,\"welds_per_s\":%.0f,\"min\":%.4f,\"max\":%.4f,\"hash\":\"%016llx\"}",
                      scaling.empty() ? "" : ",", th, t, pts.size() / t, *std::min_element(ts.begin(), ts.end()),
                      *std::max_element(ts.begin(), ts.end()), static_cast<unsigned long long>(h));
        scaling += buf;
    }

    // Single-thread ablation, interleaved: each round runs all four variants
    // (rotating the order), so slow drift on a shared machine hits them all
    // alike. Ratios are medians of the per-round paired ratios.
    set_threads(1);
    const char* names[4] = {"full", "no_zorder", "no_prune", "neither"};
    ConnectOptions cfg[4] = {o, o, o, o};
    cfg[1].spatial_order = false;
    cfg[2].prune = false;
    cfg[3].spatial_order = false, cfg[3].prune = false;
    int ab = static_cast<int>(a.num("ab-reps", 11));
    std::vector<double> tv[4];
    for (int r = 0; r < ab; ++r)
        for (int j = 0; j < 4; ++j) {
            int c = (j + r) % 4;
            t0 = Clock::now();
            auto res = connect_all(m, &bvh, pts, cfg[c]);
            tv[c].push_back(secs(t0));
            all_same = all_same && hash_results(res) == h0;
        }
    std::string ablation;
    for (int c = 0; c < 4; ++c) {
        std::vector<double> ratio;
        for (int r = 0; r < ab; ++r) ratio.push_back(tv[c][r] / tv[0][r]);
        char buf[200];
        std::snprintf(buf, sizeof buf, "%s{\"variant\":\"%s\",\"median_secs\":%.4f,\"min\":%.4f,\"max\":%.4f,\"slowdown_vs_full\":%.2f}",
                      c ? "," : "", names[c], median(tv[c]), *std::min_element(tv[c].begin(), tv[c].end()),
                      *std::max_element(tv[c].begin(), tv[c].end()), median(ratio));
        ablation += buf;
    }
    t1 = median(tv[0]);

    // Brute force on a sample: same answers, and how much slower.
    size_t nb = std::min<size_t>(static_cast<size_t>(a.num("brute", 200)), pts.size());
    t0 = Clock::now();
    size_t mismatch = 0;
    for (size_t i = 0; i < nb; ++i) mismatch += !same(connect_one(m, nullptr, pts[i], o), ref[i]);
    double t_brute = secs(t0) / static_cast<double>(nb);
    double t_bvh_one = t1 / static_cast<double>(pts.size());

    // How many elements does a query reach, and what does traversal alone cost?
    size_t tested = 0;
    t0 = Clock::now();
    for (size_t i = 0; i < pts.size(); ++i) bvh.query(pts[i], o.tol, [&](int32_t, const Aabb&) { ++tested; });
    double t_trav = secs(t0);

    std::printf("{\"cmd\":\"bench\",\"nodes\":%zu,\"elems\":%zu,\"welds\":%zu,\"file_mb\":%.1f,\"read_secs\":%.3f,"
                "\"read_mb_per_s\":%.1f,\"bvh_secs\":%.3f,\"bvh_nodes\":%zu,\"bvh_depth\":%d,"
                "\"elems_tested_per_query\":%.1f,\"traversal_only_1t_secs\":%.4f,\"scaling\":[%s],\"identical_across_threads\":%s,"
                "\"brute\":{\"sample\":%zu,\"mismatches\":%zu,\"us_per_weld\":%.1f},\"bvh_us_per_weld_1t\":%.2f,"
                "\"speedup_vs_brute_1t\":%.0f,\"ablation_1t\":[%s],\"hw_threads\":%d}\n",
                m.nodes(), m.elems(), pts.size(), mb, t_read, mb / t_read, t_bvh, bvh.node_count(), bvh.depth(),
                static_cast<double>(tested) / static_cast<double>(pts.size()), t_trav, scaling.c_str(),
                all_same ? "true" : "false", nb, mismatch, t_brute * 1e6, t_bvh_one * 1e6, t_brute / t_bvh_one,
                ablation.c_str(), hw_threads());
    return mismatch || !all_same ? 1 : 0;
}

// Projection accuracy on random warped quads against an analytic answer:
// p = x(xi0,eta0) + d * n(xi0,eta0) projects back to x(xi0,eta0).
int cmd_proj(const Args& a) {
    int n = static_cast<int>(a.num("n", 200000));
    Rng rng(static_cast<uint64_t>(a.num("seed", 7)));
    std::string out;
    for (double warp : {0.0, 0.02, 0.05, 0.1, 0.2}) {
        std::vector<double> en, es, wn, ws;
        int newton_fail = 0;
        for (int i = 0; i < n; ++i) {
            const double h = 2.0;
            std::array<Vec3, 4> x = {Vec3{0, 0, 0}, Vec3{h, 0, 0}, Vec3{h, h, 0}, Vec3{0, h, 0}};
            for (auto& c : x) c = c + Vec3{rng.uni(-0.3, 0.3), rng.uni(-0.3, 0.3), 0};
            double wz = warp * h;  // twist: opposite corners up/down
            x[0].z += wz, x[1].z -= wz, x[2].z += wz, x[3].z -= wz;
            double xi = rng.uni(-0.95, 0.95), eta = rng.uni(-0.95, 0.95);
            Vec3 x0 = quad_point(x, xi, eta);
            Vec3 tx = 0.25 * ((1 - eta) * (x[1] - x[0]) + (1 + eta) * (x[2] - x[3]));
            Vec3 te = 0.25 * ((1 - xi) * (x[3] - x[0]) + (1 + xi) * (x[2] - x[1]));
            Vec3 nrm = cross(tx, te);
            nrm = (1.0 / norm(nrm)) * nrm;
            Vec3 p = x0 + rng.uni(-0.5, 0.5) * nrm;
            auto w0 = quad_shape(xi, eta);
            Projection pn = project_quad(p, x), ps = project_quad_split(p, x);
            newton_fail += !pn.normal;
            en.push_back(norm(pn.point - x0));
            es.push_back(norm(ps.point - x0));
            double dn = 0, ds = 0;
            for (int k = 0; k < 4; ++k) dn = std::max(dn, std::fabs(pn.w[k] - w0[k])), ds = std::max(ds, std::fabs(ps.w[k] - w0[k]));
            wn.push_back(dn);
            ws.push_back(ds);
        }
        auto q = [](std::vector<double> v, double f) {
            std::sort(v.begin(), v.end());
            return v[static_cast<size_t>(f * (v.size() - 1))];
        };
        char buf[400];
        std::snprintf(buf, sizeof buf,
                      "%s{\"warp\":%g,\"newton\":{\"pos_max\":%.2e,\"w_max\":%.2e,\"not_normal\":%d},"
                      "\"split\":{\"pos_p50\":%.2e,\"pos_p99\":%.2e,\"pos_max\":%.2e,\"w_p50\":%.3f,\"w_max\":%.3f}}",
                      out.empty() ? "" : ",", warp, q(en, 1), q(wn, 1), newton_fail, q(es, .5), q(es, .99), q(es, 1),
                      q(ws, .5), q(ws, 1));
        out += buf;
    }
    std::printf("{\"cmd\":\"proj\",\"n\":%d,\"element_mm\":2,\"levels\":[%s]}\n", n, out.c_str());
    return 0;
}

int cmd_lap(const Args& a) {
    if (a.pos.empty()) throw std::runtime_error("lap <out.inp> [--hA --hB --mode shape|nearest|none]");
    LapSpec s;
    s.hA = a.num("hA", 2.5);
    s.hB = a.num("hB", 2.0);
    s.mode = a.str("mode", "shape");
    std::ofstream f(a.pos[0]);
    LapInfo i = write_lap_inp(f, s);
    std::printf("{\"cmd\":\"lap\",\"mode\":\"%s\",\"hA\":%g,\"hB\":%g,\"nodes\":%d,\"elems\":%d,\"welds\":%d,"
                "\"equations\":%d,\"max_snap_mm\":%.3f,\"patch_nodes\":%d}\n",
                s.mode.c_str(), s.hA, s.hB, i.nodes, i.elems, i.welds, i.equations, i.max_snap, i.patch_nodes);
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr,
                     "usage: mesh-connect <gen|check|connect|eval|bench|proj|lap> ...\n"
                     "  gen <dir> [--parts 8 --h 2 --welds 20000 --seed 1]\n"
                     "  check <mesh.bdf> [--merge-tol 1e-3]\n"
                     "  connect <mesh.bdf> <welds.csv> [--tol 1.2 --max-sheets 3 --threads N --out cweld.bdf --report r.csv]\n"
                     "  eval <mesh.bdf> <labels.csv> [--tol 1.2]\n"
                     "  bench <mesh.bdf> <welds.csv> [--threads 1,2,3,4 --reps 5 --brute 200]\n"
                     "  proj [--n 200000]\n"
                     "  lap <out.inp> [--hA 2.5 --hB 2.0 --mode shape|nearest|patch|none]\n");
        return 2;
    }
    std::string c = argv[1];
    Args a = parse(argc, argv);
    try {
        if (c == "gen") return cmd_gen(a);
        if (c == "check") return cmd_check(a);
        if (c == "connect") return cmd_connect(a);
        if (c == "eval") return cmd_eval(a);
        if (c == "bench") return cmd_bench(a);
        if (c == "proj") return cmd_proj(a);
        if (c == "lap") return cmd_lap(a);
        std::fprintf(stderr, "unknown command %s\n", c.c_str());
        return 2;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
}
