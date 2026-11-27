#include "mc/c_api.h"

#include <exception>
#include <fstream>
#include <memory>
#include <string>
#include <unordered_map>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "mc/bvh.hpp"
#include "mc/connect.hpp"
#include "mc/nastran.hpp"
#include "mc/topology.hpp"

struct mc_mesh {
    mc::Mesh m;
    std::unique_ptr<mc::Bvh> bvh;
    std::unordered_map<int64_t, int32_t> elem_index;
};

namespace {

thread_local std::string g_err;

template <class R, class F>
R guard(R fail, F&& f) {
    try {
        g_err.clear();
        return f();
    } catch (const std::exception& e) {
        g_err = e.what();
    } catch (...) {
        g_err = "unknown error";
    }
    return fail;
}

// mc_weld (C, solver ids) -> WeldResult (internal indices)
std::vector<mc::WeldResult> to_internal(const mc_mesh* h, const mc_weld* w, int64_t n) {
    std::vector<mc::WeldResult> r(static_cast<size_t>(n));
    for (int64_t i = 0; i < n; ++i) {
        r[i].status = static_cast<mc::WeldStatus>(w[i].status);
        r[i].nsheets = w[i].nsheets;
        for (int s = 0; s < w[i].nsheets && s < 3; ++s) {
            auto it = h->elem_index.find(w[i].elem_id[s]);
            if (it == h->elem_index.end()) throw std::runtime_error("unknown element id " + std::to_string(w[i].elem_id[s]));
            r[i].sheet[s].elem = it->second;
            r[i].sheet[s].pid = w[i].pid[s];
        }
    }
    return r;
}

}  // namespace

extern "C" {

const char* mc_version(void) { return "mesh-connect 0.1.0"; }
const char* mc_last_error(void) { return g_err.c_str(); }

mc_mesh* mc_load_nastran(const char* path) {
    return guard<mc_mesh*>(nullptr, [&]() -> mc_mesh* {
        auto h = std::make_unique<mc_mesh>();
        h->m = mc::read_nastran_file(path);
        for (size_t e = 0; e < h->m.elems(); ++e) h->elem_index[h->m.elem_id[e]] = static_cast<int32_t>(e);
        return h.release();
    });
}

void mc_free(mc_mesh* m) { delete m; }
int64_t mc_node_count(const mc_mesh* m) { return m ? static_cast<int64_t>(m->m.nodes()) : -1; }
int64_t mc_elem_count(const mc_mesh* m) { return m ? static_cast<int64_t>(m->m.elems()) : -1; }

int32_t mc_components(const mc_mesh* m) {
    return guard<int32_t>(-1, [&] { return mc::components(m->m).count; });
}

int64_t mc_connect(mc_mesh* m, const double* xyz, int64_t n, double tol, int32_t max_sheets, int32_t threads,
                   mc_weld* out) {
    return guard<int64_t>(-1, [&]() -> int64_t {
        if (!m || (!xyz && n) || (!out && n) || n < 0) throw std::runtime_error("bad arguments");
        if (!m->bvh) m->bvh = std::make_unique<mc::Bvh>(m->m);
#ifdef _OPENMP
        if (threads > 0) omp_set_num_threads(threads);
#else
        (void)threads;
#endif
        std::vector<mc::Vec3> pts(static_cast<size_t>(n));
        for (int64_t i = 0; i < n; ++i) pts[i] = {xyz[3 * i], xyz[3 * i + 1], xyz[3 * i + 2]};
        mc::ConnectOptions o;
        o.tol = tol;
        o.max_sheets = max_sheets;
        auto res = mc::connect_all(m->m, m->bvh.get(), pts, o);
        int64_t ok = 0;
        for (int64_t i = 0; i < n; ++i) {
            mc_weld& w = out[i];
            w = mc_weld{};
            w.status = static_cast<int32_t>(res[i].status);
            w.nsheets = res[i].nsheets;
            ok += res[i].status == mc::WeldStatus::Ok;
            for (int s = 0; s < res[i].nsheets; ++s) {
                const auto& h = res[i].sheet[s];
                w.elem_id[s] = m->m.elem_id[h.elem];
                w.pid[s] = h.pid;
                w.dist[s] = h.dist;
                w.point[s][0] = h.point.x, w.point[s][1] = h.point.y, w.point[s][2] = h.point.z;
                for (int k = 0; k < 4; ++k) w.weight[s][k] = h.w[k];
            }
        }
        return ok;
    });
}

int32_t mc_components_with(const mc_mesh* m, const mc_weld* welds, int64_t n) {
    return guard<int32_t>(-1, [&] {
        auto r = to_internal(m, welds, n);
        return mc::components(m->m, r).count;
    });
}

int64_t mc_write_cwelds(const mc_mesh* m, const char* path, const double* xyz, const mc_weld* welds, int64_t n,
                        int64_t first_id, int32_t pwid) {
    return guard<int64_t>(-1, [&]() -> int64_t {
        auto r = to_internal(m, welds, n);
        std::vector<mc::Vec3> pts(static_cast<size_t>(n));
        for (int64_t i = 0; i < n; ++i) pts[i] = {xyz[3 * i], xyz[3 * i + 1], xyz[3 * i + 2]};
        std::ofstream f(path);
        if (!f) throw std::runtime_error(std::string("cannot write ") + path);
        return static_cast<int64_t>(mc::write_cwelds(f, m->m, pts, r, first_id, pwid));
    });
}

}  // extern "C"
