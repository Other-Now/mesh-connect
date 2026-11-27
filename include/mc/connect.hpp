#pragma once
// Spot-weld connector realisation: for each weld location, find the closest
// element on each nearby part (within a search tolerance), project the weld
// onto it, and keep the 2..max_sheets closest parts as the welded stack.
#include <array>
#include <cstdint>
#include <span>
#include <vector>

#include "mc/bvh.hpp"
#include "mc/mesh.hpp"

namespace mc {

struct SheetHit {
    int32_t elem = -1;  // element index
    int32_t pid = 0;
    double dist = 0;    // weld point to projected point
    Vec3 point;         // projected point on the element
    std::array<double, 4> w{0, 0, 0, 0};  // node weights at that point
    bool normal = true;  // false: weld lies past the element/part edge
};

enum class WeldStatus : int32_t {
    Ok = 0,        // 2+ sheets found
    OneSheet = 1,  // only one part within tolerance: nothing to connect
    NoSheet = 2,   // no part within tolerance
};

struct WeldResult {
    WeldStatus status = WeldStatus::NoSheet;
    int32_t nsheets = 0;
    // 2 sheets: closest first. 3 sheets: ordered through the stack, so
    // sheet[0]-sheet[1] and sheet[1]-sheet[2] are the welded interfaces.
    std::array<SheetHit, 3> sheet;
};

struct ConnectOptions {
    double tol = 1.0;     // max distance from weld point to a sheet
    int max_sheets = 3;   // 2- or 3-sheet stacks
    bool spatial_order = true;  // process welds in Morton (Z-curve) order for cache locality
    bool prune = true;          // skip elements whose box is farther than the part's best hit
};

// bvh == nullptr: brute force over every element (used as the reference).
WeldResult connect_one(const Mesh& m, const Bvh* bvh, Vec3 p, const ConnectOptions& o);

// All welds; parallel over welds with OpenMP. Each weld is independent and
// written to its own slot, so results do not depend on the thread count.
std::vector<WeldResult> connect_all(const Mesh& m, const Bvh* bvh, std::span<const Vec3> pts,
                                    const ConnectOptions& o);

}  // namespace mc
