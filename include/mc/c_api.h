/* mesh-connect C ABI: a stable, C-only surface for the shared library
 * (libmeshconnect.so / meshconnect.dll). No C++ types cross it, so it can be
 * loaded from Python (ctypes), C#, or an application built with another
 * compiler. Errors: functions return < 0 / NULL; mc_last_error() says why. */
#ifndef MC_C_API_H
#define MC_C_API_H
#include <stdint.h>

#if defined(_WIN32)
#  if defined(MC_BUILDING_DLL)
#    define MC_API __declspec(dllexport)
#  else
#    define MC_API __declspec(dllimport)
#  endif
#else
#  define MC_API __attribute__((visibility("default")))
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef struct mc_mesh mc_mesh;

typedef struct mc_weld {
    int32_t status;       /* 0 ok, 1 only one sheet, 2 no sheet */
    int32_t nsheets;
    int64_t elem_id[3];   /* solver element ids of the connected elements */
    int32_t pid[3];
    double dist[3];       /* weld point to projected point */
    double point[3][3];   /* projected point on each sheet */
    double weight[3][4];  /* node weights on each element */
} mc_weld;

MC_API const char* mc_version(void);
MC_API const char* mc_last_error(void); /* thread-local */

MC_API mc_mesh* mc_load_nastran(const char* path);
MC_API void mc_free(mc_mesh* m);
MC_API int64_t mc_node_count(const mc_mesh* m);
MC_API int64_t mc_elem_count(const mc_mesh* m);

/* Connected components of the mesh alone (elements sharing nodes). */
MC_API int32_t mc_components(const mc_mesh* m);

/* Realise n welds (xyz packed as 3*n doubles) into out[n]. threads <= 0: default.
 * Returns the number of welds with status 0, or -1 on error. */
MC_API int64_t mc_connect(mc_mesh* m, const double* xyz, int64_t n, double tol, int32_t max_sheets,
                          int32_t threads, mc_weld* out);

/* Components after adding the n connectors from mc_connect. */
MC_API int32_t mc_components_with(const mc_mesh* m, const mc_weld* welds, int64_t n);

/* Write GRID + CWELD cards. Returns CWELD count, or -1 on error. */
MC_API int64_t mc_write_cwelds(const mc_mesh* m, const char* path, const double* xyz, const mc_weld* welds,
                               int64_t n, int64_t first_id, int32_t pwid);

#ifdef __cplusplus
}
#endif
#endif
