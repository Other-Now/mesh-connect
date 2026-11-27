#pragma once
// Synthetic test assemblies with ground truth. Planar sheets make the true
// projection of a weld analytic, so the labels do not depend on the code
// under test: the right answer is (x, y, z_sheet) on each expected part.
#include <cstdint>
#include <string>
#include <vector>

#include "mc/mesh.hpp"

namespace mc {

struct Plate {
    int32_t pid;
    double x0, x1, y0, y1, z;
    double h;             // target element size
    bool tri = false;     // CTRIA3 (split quads) instead of CQUAD4
    double jitter = 0.0;  // in-plane random offset of interior nodes, fraction of h
    double t = 1.0;       // thickness (PSHELL)
};
void add_plate(Mesh& m, const Plate& p, uint64_t seed, int64_t& next_node, int64_t& next_elem);

struct Label {
    Vec3 p;
    std::string kind;            // "2sheet", "3sheet", "1sheet"
    std::vector<int32_t> pids;   // expected parts, sorted
    std::vector<double> zs;      // sheet z for each expected pid (same order)
};

struct AssemblySpec {
    int parts = 8;       // panels in a row; neighbours overlap in lap-joint flanges
    double h = 2.0;      // base element size (mm); varies per part so meshes don't match
    int welds = 20000;   // labelled weld locations
    uint64_t seed = 1;
    double gap = 1.0;    // mid-surface distance between welded sheets
};
struct Assembly {
    Mesh mesh;
    std::vector<Label> labels;
    std::vector<Plate> plates;
};
Assembly gen_assembly(const AssemblySpec& s);

// Deterministic RNG (splitmix64), so generated meshes are identical everywhere.
struct Rng {
    uint64_t s;
    explicit Rng(uint64_t seed) : s(seed) {}
    uint64_t next() {
        uint64_t z = (s += 0x9E3779B97F4A7C15ull);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
        return z ^ (z >> 31);
    }
    double uni(double a = 0, double b = 1) { return a + (b - a) * (next() >> 11) * 0x1.0p-53; }
};

}  // namespace mc
