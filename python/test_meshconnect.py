"""Smoke test of the shared library through Python ctypes (runs on Linux and Windows CI)."""
import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import meshconnect  # noqa: E402


def plate(lines, pid, x0, z, nid, eid, n=10, h=1.0):
    """n x n quads, h mm, at height z, starting at x0. Small-field GRID cards."""
    first = nid
    for j in range(n + 1):
        for i in range(n + 1):
            lines.append("GRID    %8d        %8.3f%8.3f%8.3f" % (nid, x0 + i * h, j * h, z))
            nid += 1
    for j in range(n):
        for i in range(n):
            a = first + j * (n + 1) + i
            lines.append("CQUAD4  %8d%8d%8d%8d%8d%8d" % (eid, pid, a, a + 1, a + n + 2, a + n + 1))
            eid += 1
    return nid, eid


def main():
    lines = ["BEGIN BULK"]
    nid, eid = plate(lines, 1, 0.0, 0.0, 1, 1)
    nid, eid = plate(lines, 2, 7.0, 1.0, nid, eid)  # overlaps x in [7, 10]
    lines.append("ENDDATA")
    d = tempfile.mkdtemp()
    path = os.path.join(d, "lap.bdf")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")

    print(meshconnect.version())
    with meshconnect.Mesh(path) as m:
        assert m.nodes == 242 and m.elems == 200, (m.nodes, m.elems)
        assert m.components() == 2
        welds = m.connect([(8.3, 2.7, 0.5), (8.6, 7.1, 0.5), (3.0, 3.0, 0.5)], tol=0.6)
        assert welds[0]["status"] == "ok" and welds[1]["status"] == "ok"
        assert welds[2]["status"] == "one_sheet"
        s = welds[0]["sheets"]
        assert sorted(x["pid"] for x in s) == [1, 2]
        assert all(abs(x["dist"] - 0.5) < 1e-12 for x in s)
        assert all(abs(sum(x["weights"]) - 1) < 1e-12 for x in s)
        assert m.components(welds) == 1
        n = m.write_cwelds(os.path.join(d, "cweld.bdf"), welds, first_id=100000)
        assert n == 2
        try:
            meshconnect.Mesh(os.path.join(d, "missing.bdf"))
            raise AssertionError("expected an error")
        except RuntimeError as e:
            assert "cannot open" in str(e)
    print("python binding ok")


if __name__ == "__main__":
    main()
