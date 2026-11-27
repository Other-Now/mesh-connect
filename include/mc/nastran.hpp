#pragma once
// Nastran bulk-data I/O for shell meshes: GRID, CQUAD4, CTRIA3, PSHELL in
// small-field (8-char), large-field (GRID*, 16-char) and free-field (comma)
// formats, with continuation lines. Other cards are counted and skipped.
#include <iosfwd>
#include <map>
#include <span>
#include <string>
#include <string_view>

#include "mc/connect.hpp"
#include "mc/mesh.hpp"

namespace mc {

struct NastranStats {
    size_t lines = 0, cards = 0, grids = 0, quads = 0, trias = 0, pshells = 0;
    std::map<std::string, size_t> skipped;  // unsupported card name -> count
};

Mesh read_nastran(std::istream& in, NastranStats* st = nullptr);
Mesh read_nastran_file(const std::string& path, NastranStats* st = nullptr);

// Nastran real field: "1.5", "-.5", "1.5-3" (= 1.5e-3), "2.+4", "1.D2", "7".
double parse_real(std::string_view s);

// Free-field, full precision (shortest round-trip), always with a decimal point.
void write_nastran(std::ostream& out, const Mesh& m);

// One GRID at the weld location plus CWELD (TYPE=ELEMID) cards: one per
// adjacent sheet pair (A-B, and B-C for a 3-sheet stack). Failed welds are
// written as comments. Returns the number of CWELD cards written.
size_t write_cwelds(std::ostream& out, const Mesh& m, std::span<const Vec3> pts,
                    std::span<const WeldResult> res, int64_t first_id, int32_t pwid);

std::string fmt_real(double v);

}  // namespace mc
