#include "mc/nastran.hpp"

#include <algorithm>
#include <cctype>
#include <charconv>
#include <cstring>
#include <fstream>
#include <istream>
#include <ostream>
#include <stdexcept>
#include <vector>

namespace mc {
namespace {

std::string_view trim(std::string_view s) {
    while (!s.empty() && std::isspace(static_cast<unsigned char>(s.front()))) s.remove_prefix(1);
    while (!s.empty() && std::isspace(static_cast<unsigned char>(s.back()))) s.remove_suffix(1);
    return s;
}

std::string upper(std::string_view s) {
    std::string r(s);
    for (auto& c : r) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return r;
}

int64_t parse_int(std::string_view s, int64_t dflt = 0) {
    s = trim(s);
    if (s.empty()) return dflt;
    int64_t v = 0;
    auto [p, ec] = std::from_chars(s.data(), s.data() + s.size(), v);
    if (ec != std::errc() || p != s.data() + s.size())
        throw std::runtime_error("bad integer field '" + std::string(s) + "'");
    return v;
}

enum class Fmt { Small, Large, Free };

// Fields of one physical line (excluding the name / continuation field).
void line_fields(std::string_view line, Fmt f, std::vector<std::string>& out) {
    if (f == Fmt::Free) {
        size_t start = line.find(',');  // skip name / continuation marker
        while (start != std::string_view::npos) {
            size_t next = line.find(',', start + 1);
            out.emplace_back(trim(line.substr(start + 1, next == std::string_view::npos ? next : next - start - 1)));
            start = next;
        }
        return;
    }
    const size_t w = f == Fmt::Large ? 16 : 8, n = f == Fmt::Large ? 4 : 8;
    for (size_t i = 0; i < n; ++i) {
        size_t b = 8 + i * w;
        out.emplace_back(b < line.size() ? trim(line.substr(b, w)) : std::string_view{});
    }
}

bool is_continuation(std::string_view line) {
    return !line.empty() && (line[0] == '+' || line[0] == '*' || line[0] == ',' ||
                             (line.size() > 8 && trim(line.substr(0, 8)).empty() && !trim(line).empty()));
}

struct RawElem {
    int64_t id, pid, n[4];
    int nn;
};

}  // namespace

double parse_real(std::string_view s) {
    s = trim(s);
    if (s.empty()) return 0;
    char buf[64];
    size_t n = 0;
    for (size_t i = 0; i < s.size() && n < sizeof(buf) - 2; ++i) {
        char c = s[i];
        if (c == 'D' || c == 'd') c = 'E';
        // implicit exponent: a sign after a digit or '.' that is not after 'E'
        if ((c == '+' || c == '-') && i > 0 && n > 0 && buf[n - 1] != 'E' && buf[n - 1] != 'e')
            buf[n++] = 'E';
        buf[n++] = c;
    }
    const char* b = buf;
    if (*b == '+') ++b;  // from_chars rejects a leading '+'
    double v = 0;
    auto [p, ec] = std::from_chars(b, buf + n, v);
    if (ec != std::errc() || p != buf + n)
        throw std::runtime_error("bad real field '" + std::string(s) + "'");
    return v;
}

Mesh read_nastran(std::istream& in, NastranStats* st) {
    NastranStats local;
    NastranStats& S = st ? *st : local;
    Mesh m;
    std::vector<RawElem> raw;
    std::vector<std::string> f;
    std::string name, line;

    auto flush = [&]() {
        if (name.empty()) return;
        S.cards++;
        auto F = [&](size_t i) -> std::string_view { return i < f.size() ? std::string_view(f[i]) : std::string_view{}; };
        if (name == "GRID") {
            int64_t id = parse_int(F(0));
            if (parse_int(F(1)) != 0) throw std::runtime_error("GRID " + std::to_string(id) + ": CP != 0 not supported");
            m.add_node(id, {parse_real(F(2)), parse_real(F(3)), parse_real(F(4))});
            S.grids++;
        } else if (name == "CQUAD4" || name == "CTRIA3") {
            RawElem r{};
            r.id = parse_int(F(0));
            r.pid = parse_int(F(1), r.id);  // blank PID defaults to EID
            r.nn = name == "CQUAD4" ? 4 : 3;
            for (int i = 0; i < r.nn; ++i) r.n[i] = parse_int(F(2 + i));
            raw.push_back(r);
            (r.nn == 4 ? S.quads : S.trias)++;
        } else if (name == "PSHELL") {
            m.thickness[static_cast<int32_t>(parse_int(F(0)))] = parse_real(F(2));
            S.pshells++;
        } else {
            S.skipped[name]++;
        }
        name.clear();
        f.clear();
    };

    while (std::getline(in, line)) {
        S.lines++;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (auto c = line.find('$'); c != std::string::npos) line.resize(c);
        std::string_view lv = line;
        if (trim(lv).empty()) continue;
        Fmt fmt = lv.find(',') != std::string_view::npos ? Fmt::Free : Fmt::Small;
        if (is_continuation(lv)) {
            if (name.empty()) continue;  // continuation of a skipped header line
            if (fmt == Fmt::Small && lv[0] == '*') fmt = Fmt::Large;
            line_fields(lv, fmt, f);
            continue;
        }
        flush();
        std::string_view head = fmt == Fmt::Free ? lv.substr(0, lv.find(',')) : lv.substr(0, std::min<size_t>(8, lv.size()));
        head = trim(head);
        if (!head.empty() && head.back() == '*') {
            head.remove_suffix(1);
            if (fmt == Fmt::Small) fmt = Fmt::Large;
        }
        name = upper(head);
        if (name.find(' ') != std::string::npos || name == "BEGIN" || name == "ENDDATA") {
            if (name == "ENDDATA") {
                name.clear();
                break;
            }
            name.clear();  // executive/case control or "BEGIN BULK"
            continue;
        }
        line_fields(lv, fmt, f);
    }
    flush();

    m.elem_id.reserve(raw.size());
    m.pid.reserve(raw.size());
    m.conn.reserve(raw.size());
    for (const auto& r : raw) {
        int32_t idx[4] = {-1, -1, -1, -1};
        for (int i = 0; i < r.nn; ++i) {
            auto it = m.node_index.find(r.n[i]);
            if (it == m.node_index.end())
                throw std::runtime_error("element " + std::to_string(r.id) + " references missing GRID " +
                                         std::to_string(r.n[i]));
            idx[i] = it->second;
        }
        m.add_elem(r.id, static_cast<int32_t>(r.pid), idx[0], idx[1], idx[2], idx[3]);
    }
    return m;
}

Mesh read_nastran_file(const std::string& path, NastranStats* st) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("cannot open " + path);
    std::vector<char> buf(1 << 20);
    in.rdbuf()->pubsetbuf(buf.data(), static_cast<std::streamsize>(buf.size()));
    return read_nastran(in, st);
}

std::string fmt_real(double v) {
    char b[40];
    auto [p, ec] = std::to_chars(b, b + sizeof(b), v);  // shortest round-trip
    std::string s(b, p);
    if (s.find_first_of(".eEni") == std::string::npos) s += '.';
    else if (s.find('.') == std::string::npos && s.find('e') != std::string::npos)
        s.insert(s.find('e'), ".");
    return s;
}

void write_nastran(std::ostream& out, const Mesh& m) {
    out << "BEGIN BULK\n";
    std::vector<int32_t> pids(m.pid);
    std::sort(pids.begin(), pids.end());
    pids.erase(std::unique(pids.begin(), pids.end()), pids.end());
    for (int32_t p : pids) {
        auto it = m.thickness.find(p);
        if (it != m.thickness.end()) out << "PSHELL," << p << ",1," << fmt_real(it->second) << "\n";
    }
    std::string s;
    for (size_t i = 0; i < m.nodes(); ++i) {
        s.clear();
        s += "GRID,";
        s += std::to_string(m.node_id[i]);
        s += ",,";
        s += fmt_real(m.xyz[i].x) + "," + fmt_real(m.xyz[i].y) + "," + fmt_real(m.xyz[i].z) + "\n";
        out << s;
    }
    for (size_t e = 0; e < m.elems(); ++e) {
        out << (m.nn(e) == 4 ? "CQUAD4," : "CTRIA3,") << m.elem_id[e] << ',' << m.pid[e];
        for (int i = 0; i < m.nn(e); ++i) out << ',' << m.node_id[m.conn[e][i]];
        out << '\n';
    }
    out << "ENDDATA\n";
}

size_t write_cwelds(std::ostream& out, const Mesh& m, std::span<const Vec3> pts,
                    std::span<const WeldResult> res, int64_t first_id, int32_t pwid) {
    size_t n = 0;
    out << "$ mesh-connect spot welds: CWELD TYPE=ELEMID, one per welded sheet pair\n";
    out << "PWELD," << pwid << ",1,5.\n";
    for (size_t i = 0; i < pts.size(); ++i) {
        int64_t gid = first_id + static_cast<int64_t>(i);
        const WeldResult& r = res[i];
        if (r.status != WeldStatus::Ok) {
            out << "$ weld " << i << " not realised: "
                << (r.status == WeldStatus::OneSheet ? "only one sheet within tolerance" : "no sheet within tolerance")
                << '\n';
            continue;
        }
        out << "GRID," << gid << ",," << fmt_real(pts[i].x) << ',' << fmt_real(pts[i].y) << ','
            << fmt_real(pts[i].z) << '\n';
        for (int s = 0; s + 1 < r.nsheets; ++s) {
            int64_t eid = first_id + 2 * static_cast<int64_t>(i) + s;
            out << "CWELD," << eid << ',' << pwid << ',' << gid << ",ELEMID,,,,\n,"
                << m.elem_id[r.sheet[s].elem] << ',' << m.elem_id[r.sheet[s + 1].elem] << '\n';
            ++n;
        }
    }
    return n;
}

}  // namespace mc
