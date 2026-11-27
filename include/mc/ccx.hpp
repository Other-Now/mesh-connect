#pragma once
// CalculiX (.inp) export of a two-plate lap joint, to check that connectors
// actually carry load in a finite-element solve, not just that they look
// right geometrically.
#include <iosfwd>
#include <string>

namespace mc {

struct LapSpec {
    double hA = 2.5, hB = 2.0;  // element sizes: different, so the meshes don't line up
    // "shape":   weld tied to each sheet with the element shape functions
    //            (u_A(x_w) = u_B(x_w), interpolated: the CWELD/ELEMID idea)
    // "nearest": tie the nearest node of A to the nearest node of B
    // "patch":   tie the mean displacement of all nodes within the nugget
    //            radius on each sheet (spreads the weld over its diameter)
    // "none":    no connectors (plate A is not attached to anything)
    std::string mode = "shape";
    double load = 1000.0;   // N, total, pulling plate A in -x
    double nugget = 5.0;    // weld diameter for "patch" (mm)
};

struct LapInfo {
    int nodes = 0, elems = 0, welds = 0, equations = 0;
    double max_snap = 0;   // "nearest" mode: largest node-to-weld distance (mm)
    int patch_nodes = 0;   // "patch" mode: most nodes in one weld's patch (both sheets)
};

LapInfo write_lap_inp(std::ostream& out, const LapSpec& s);

}  // namespace mc
