// Unit + integration tests. Tiny harness: no external test framework needed.
#include <algorithm>
#include <cmath>
#include <fstream>
#include <cstdio>
#include <cstring>
#include <functional>
#include <sstream>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "mc/bvh.hpp"
#include "mc/c_api.h"
#include "mc/ccx.hpp"
#include "mc/connect.hpp"
#include "mc/gen.hpp"
#include "mc/nastran.hpp"
#include "mc/topology.hpp"

using namespace mc;

static int g_fail = 0, g_checks = 0;
#define CHECK(c)                                                              \
    do {                                                                      \
        ++g_checks;                                                           \
        if (!(c)) {                                                           \
            ++g_fail;                                                         \
            std::printf("  FAIL %s:%d  %s\n", __FILE__, __LINE__, #c);        \
        }                                                                     \
    } while (0)
#define NEAR(a, b, tol) CHECK(std::fabs((a) - (b)) <= (tol))

static std::vector<std::pair<const char*, std::function<void()>>>& tests() {
    static std::vector<std::pair<const char*, std::function<void()>>> t;
    return t;
}
struct Reg {
    Reg(const char* n, std::function<void()> f) { tests().push_back({n, std::move(f)}); }
};
#define TEST(name)                       \
    static void name();                  \
    static Reg reg_##name(#name, name);  \
    static void name()

TEST(triangle_regions) {
    Vec3 a{0, 0, 0}, b{2, 0, 0}, c{0, 2, 0};
    auto in = project_tri({0.5, 0.5, 3}, a, b, c);
    CHECK(in.normal);
    NEAR(in.dist2, 9.0, 1e-12);
    NEAR(in.w[0], 0.5, 1e-12); NEAR(in.w[1], 0.25, 1e-12); NEAR(in.w[2], 0.25, 1e-12);
    auto edge = project_tri({1, -1, 0}, a, b, c);  // below edge ab
    CHECK(!edge.normal);
    NEAR(edge.point.x, 1, 1e-12); NEAR(edge.point.y, 0, 1e-12); NEAR(edge.dist2, 1, 1e-12);
    auto vert = project_tri({-1, -1, 0}, a, b, c);  // vertex a region
    NEAR(vert.w[0], 1, 1e-12); NEAR(vert.dist2, 2, 1e-12);
}

TEST(quad_newton_recovers_parameters_on_warped_quads) {
    Rng rng(3);
    double worst = 0;
    for (int i = 0; i < 20000; ++i) {
        std::array<Vec3, 4> x = {Vec3{0, 0, 0}, Vec3{2, 0, 0}, Vec3{2, 2, 0}, Vec3{0, 2, 0}};
        for (auto& c : x) c = c + Vec3{rng.uni(-.3, .3), rng.uni(-.3, .3), rng.uni(-.3, .3)};
        double xi = rng.uni(-.9, .9), eta = rng.uni(-.9, .9);
        Vec3 x0 = quad_point(x, xi, eta);
        Vec3 tx = 0.25 * ((1 - eta) * (x[1] - x[0]) + (1 + eta) * (x[2] - x[3]));
        Vec3 te = 0.25 * ((1 - xi) * (x[3] - x[0]) + (1 + xi) * (x[2] - x[1]));
        Vec3 n = cross(tx, te);
        Vec3 p = x0 + (0.4 / norm(n)) * n;
        auto pr = project_quad(p, x);
        CHECK(pr.normal);
        worst = std::max(worst, norm(pr.point - x0));
    }
    CHECK(worst < 1e-10);
}

// Regression: an absolute Newton tolerance worked near the origin but not at
// x ~ 1000 mm, where doubles cannot resolve a 1e-14 step.
TEST(quad_projection_far_from_origin) {
    Rng rng(9);
    for (double off : {1e3, 1e5}) {
        int fails = 0;
        for (int i = 0; i < 2000; ++i) {
            std::array<Vec3, 4> x = {Vec3{0, 0, 0}, Vec3{2, 0, 0}, Vec3{2, 2, 0}, Vec3{0, 2, 0}};
            for (auto& c : x) c = c + Vec3{off + rng.uni(-.5, .5), off / 2 + rng.uni(-.5, .5), 0};
            Vec3 p = quad_point(x, rng.uni(-.9, .9), rng.uni(-.9, .9)) + Vec3{0, 0, 0.5};
            auto pr = project_quad(p, x);
            fails += !(pr.normal && std::fabs(std::sqrt(pr.dist2) - 0.5) < 1e-9);
        }
        CHECK(fails == 0);
    }
}

// Regression: full Newton from the centre cycled on this distorted quad
// (weld 54859 of the 1.2M-element benchmark) and the weld fell to an edge.
TEST(quad_projection_distorted_element) {
    std::array<Vec3, 4> x = {Vec3{2912.656513, 255.690950, 0}, Vec3{2914.642957, 254.709952, 0},
                             Vec3{2915.055042, 258.440679, 0}, Vec3{2911.583563, 257.682221, 0}};
    auto pr = project_quad({2914.966091124986, 258.04225481750507, 0.5}, x);
    CHECK(pr.normal);
    NEAR(std::sqrt(pr.dist2), 0.5, 1e-9);
}

TEST(quad_outside_clamps_to_edge) {
    std::array<Vec3, 4> x = {Vec3{0, 0, 0}, Vec3{1, 0, 0}, Vec3{1, 1, 0}, Vec3{0, 1, 0}};
    auto pr = project_quad({2, 0.5, 0}, x);
    CHECK(!pr.normal);
    NEAR(pr.point.x, 1, 1e-12); NEAR(pr.point.y, 0.5, 1e-12);
    NEAR(pr.w[1], 0.5, 1e-12); NEAR(pr.w[2], 0.5, 1e-12);
}

TEST(nastran_formats) {
    // small field, large field with continuation, free field, Nastran exponents,
    // comments, executive/case control noise before BEGIN BULK.
    const char* deck =
        "SOL 101\nCEND\nSUBCASE 1\nBEGIN BULK\n"
        "$ comment\n"
        "PSHELL         1       1    1.25\n"
        "GRID           1              0.      0.      0.\n"
        "GRID*                  2                           1.5-1              0.\n"
        "*                     0.\n"
        "GRID,3,,1.,1.,2.+0 $ inline comment\n"
        "GRID           4             .25     1.0   -1.-3\n"
        "CQUAD4        10       1       1       2       3       4\n"
        "CTRIA3,11,1,1,2,3\n"
        "MAT1           1  2.1+5            .3\n"
        "ENDDATA\n";
    std::istringstream in(deck);
    NastranStats st;
    Mesh m = read_nastran(in, &st);
    CHECK(m.nodes() == 4 && m.elems() == 2);
    NEAR(m.xyz[1].x, 0.15, 1e-15);
    NEAR(m.xyz[2].z, 2.0, 1e-15);
    NEAR(m.xyz[3].x, 0.25, 1e-15); NEAR(m.xyz[3].z, -1e-3, 1e-18);
    CHECK(m.nn(0) == 4 && m.nn(1) == 3);
    NEAR(m.thickness[1], 1.25, 1e-15);
    CHECK(st.skipped["MAT1"] == 1 && st.skipped.count("ENDDATA") == 0);  // ENDDATA ends the deck, is not a card
    NEAR(parse_real("-.5+2"), -50, 1e-12);
    NEAR(parse_real("1.D-2"), 0.01, 1e-15);
    NEAR(parse_real("+3."), 3, 0);

    // round trip is exact
    std::stringstream rt;
    write_nastran(rt, m);
    Mesh m2 = read_nastran(rt);
    CHECK(m2.nodes() == m.nodes() && m2.elems() == m.elems());
    for (size_t i = 0; i < m.nodes(); ++i) CHECK(std::memcmp(&m.xyz[i], &m2.xyz[i], sizeof(Vec3)) == 0);
    CHECK(m2.conn == m.conn && m2.pid == m.pid);
}

TEST(topology_components_edges_coincident) {
    Mesh m;
    // two quads sharing an edge (one part), and a third quad whose nodes
    // coincide with the second's right edge but are separate GRIDs.
    m.add_node(1, {0, 0, 0}), m.add_node(2, {1, 0, 0}), m.add_node(3, {2, 0, 0});
    m.add_node(4, {0, 1, 0}), m.add_node(5, {1, 1, 0}), m.add_node(6, {2, 1, 0});
    m.add_node(7, {2, 0, 0}), m.add_node(8, {3, 0, 0}), m.add_node(9, {3, 1, 0}), m.add_node(10, {2, 1, 1e-5});
    m.add_elem(1, 1, 0, 1, 4, 3);
    m.add_elem(2, 1, 1, 2, 5, 4);
    m.add_elem(3, 2, 6, 7, 8, 9);
    auto c = components(m);
    CHECK(c.count == 2);
    auto es = edge_stats(m);
    CHECK(es.edges == 11 && es.free == 10 && es.nonmanifold == 0);
    auto pairs = coincident_nodes(m, 1e-3);
    CHECK(pairs.size() == 2);
    CHECK(components(m, {}, pairs).count == 1);
    auto csr = node_to_elems(m);
    CHECK(csr.offset[2] - csr.offset[1] == 2);  // node 2 is in elements 1 and 2
}

static Assembly small_assembly() {
    AssemblySpec s;
    s.parts = 4;
    s.h = 5;
    s.welds = 3000;
    s.seed = 11;
    return gen_assembly(s);
}

TEST(bvh_matches_brute_force) {
    Assembly a = small_assembly();
    Bvh bvh(a.mesh);
    ConnectOptions o;
    o.tol = 1.2;
    Rng rng(5);
    int mism = 0;
    for (int i = 0; i < 400; ++i) {
        // half on labelled welds, half random points near the sheets
        Vec3 p = i % 2 ? a.labels[i].p : Vec3{rng.uni(0, 3000), rng.uni(0, 500), rng.uni(-3, 3)};
        auto x = connect_one(a.mesh, &bvh, p, o), y = connect_one(a.mesh, nullptr, p, o);
        bool same = x.status == y.status && x.nsheets == y.nsheets;
        for (int s = 0; same && s < x.nsheets; ++s) same = x.sheet[s].elem == y.sheet[s].elem && x.sheet[s].w == y.sheet[s].w;
        mism += !same;
    }
    CHECK(mism == 0);
}

TEST(labels_all_correct_and_assembly_connects) {
    Assembly a = small_assembly();
    Bvh bvh(a.mesh);
    std::vector<Vec3> pts;
    for (auto& l : a.labels) pts.push_back(l.p);
    ConnectOptions o;
    o.tol = 1.2;
    auto res = connect_all(a.mesh, &bvh, pts, o);
    int bad = 0;
    for (size_t i = 0; i < res.size(); ++i) {
        const Label& l = a.labels[i];
        std::vector<int32_t> got;
        for (int s = 0; s < res[i].nsheets; ++s) got.push_back(res[i].sheet[s].pid);
        std::sort(got.begin(), got.end());
        bool ok = got == l.pids && ((res[i].status == WeldStatus::Ok) == (l.kind != "1sheet"));
        for (int s = 0; ok && s < res[i].nsheets; ++s) {
            size_t k = std::find(l.pids.begin(), l.pids.end(), res[i].sheet[s].pid) - l.pids.begin();
            ok = norm(res[i].sheet[s].point - Vec3{l.p.x, l.p.y, l.zs[k]}) < 1e-9;
        }
        bad += !ok;
    }
    CHECK(bad == 0);
    // 4 parts + 2 patches + 1 stiffener, all separate; welds join everything
    // except the stiffener, which sits outside the tolerance.
    std::vector<WeldResult> good;
    for (auto& r : res)
        if (r.status == WeldStatus::Ok) good.push_back(r);
    CHECK(components(a.mesh).count == 7);
    auto after = components(a.mesh, good);
    CHECK(after.count == 2);
}

TEST(results_identical_across_thread_counts) {
#ifdef _OPENMP
    Assembly a = small_assembly();
    Bvh bvh(a.mesh);
    std::vector<Vec3> pts;
    for (auto& l : a.labels) pts.push_back(l.p);
    ConnectOptions o;
    omp_set_num_threads(1);
    auto r1 = connect_all(a.mesh, &bvh, pts, o);
    omp_set_num_threads(4);
    auto r4 = connect_all(a.mesh, &bvh, pts, o);
    bool same = true;
    for (size_t i = 0; i < r1.size(); ++i)
        for (int s = 0; s < r1[i].nsheets; ++s)
            same = same && r1[i].sheet[s].elem == r4[i].sheet[s].elem &&
                   std::memcmp(&r1[i].sheet[s].point, &r4[i].sheet[s].point, sizeof(Vec3)) == 0;
    CHECK(same);
#endif
}

TEST(three_sheet_stack_is_ordered) {
    Mesh m;
    int64_t nn = 1, ne = 1;
    add_plate(m, {1, 0, 10, 0, 10, 2.0, 1.0}, 1, nn, ne);  // listed top first on purpose
    add_plate(m, {2, 0, 10, 0, 10, 0.0, 1.0}, 2, nn, ne);
    add_plate(m, {3, 0, 10, 0, 10, 1.0, 1.0}, 3, nn, ne);
    ConnectOptions o;
    o.tol = 1.2;
    auto r = connect_one(m, nullptr, {5.3, 5.3, 1.0}, o);
    CHECK(r.nsheets == 3);
    CHECK(r.sheet[1].pid == 3);  // middle sheet in the middle
}

TEST(cweld_output) {
    Assembly a = small_assembly();
    std::vector<Vec3> pts = {a.labels[0].p, a.labels[1].p};
    ConnectOptions o;
    o.tol = 1.2;
    std::vector<WeldResult> r = {connect_one(a.mesh, nullptr, pts[0], o), connect_one(a.mesh, nullptr, pts[1], o)};
    std::ostringstream out;
    size_t n = write_cwelds(out, a.mesh, pts, r, 1000000, 9000);
    CHECK(n >= 1);
    CHECK(out.str().find("ELEMID") != std::string::npos);
}

TEST(lap_inp_writes_equations) {
    std::ostringstream out;
    LapSpec s;
    auto info = write_lap_inp(out, s);
    CHECK(info.welds == 6 && info.equations == 18);
    s.mode = "none";
    std::ostringstream o2;
    CHECK(write_lap_inp(o2, s).equations == 0);
}

TEST(c_api_roundtrip) {
    Assembly a = small_assembly();
    std::string path = "mc_test_mesh.bdf";
    {
        std::ofstream f(path);
        write_nastran(f, a.mesh);
    }
    mc_mesh* m = mc_load_nastran(path.c_str());
    CHECK(m != nullptr);
    if (!m) return;
    CHECK(mc_elem_count(m) == static_cast<int64_t>(a.mesh.elems()));
    std::vector<double> xyz;
    for (int i = 0; i < 100; ++i) xyz.insert(xyz.end(), {a.labels[i].p.x, a.labels[i].p.y, a.labels[i].p.z});
    std::vector<mc_weld> w(100);
    int64_t ok = mc_connect(m, xyz.data(), 100, 1.2, 3, 2, w.data());
    CHECK(ok > 50);
    CHECK(mc_components(m) == 7);
    CHECK(mc_components_with(m, w.data(), 100) < 7);
    CHECK(mc_load_nastran("does/not/exist.bdf") == nullptr && std::strlen(mc_last_error()) > 0);
    mc_free(m);
    std::remove(path.c_str());
}

int main() {
    for (auto& [name, f] : tests()) {
        int before = g_fail;
        f();
        std::printf("%s %s\n", g_fail == before ? "ok  " : "FAIL", name);
    }
    std::printf("%d checks, %d failed\n", g_checks, g_fail);
    return g_fail ? 1 : 0;
}
