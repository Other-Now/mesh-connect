# mesh-connect

Spot-weld connectors for finite-element shell meshes, in C++20. Given a Nastran
mesh of an assembly (separate sheet-metal parts) and a list of weld locations,
it finds the sheets each weld joins, projects the weld onto the right element
of each sheet, and writes the connectors (`CWELD`, `TYPE=ELEMID`). It also
reports how the model is connected before and after: components, floating
parts, free edges, and coincident nodes that were never merged.

The search uses a BVH and is parallel with OpenMP, with output identical at
any thread count. A C ABI is built as a shared library (`.so` / `.dll`) with a
ctypes Python binding. The connectors were checked in a real FE solve
(CalculiX), not just geometrically.

## What is in it

| Piece | What it does |
|---|---|
| [src/nastran.cpp](src/nastran.cpp) | Reads `GRID`, `CQUAD4`, `CTRIA3` and `PSHELL` in small-field (8-char), large-field (`GRID*`, 16-char) and free-field (comma) formats, with continuation lines and Nastran reals (`1.5-3` = 1.5e-3). Writes free-field decks with exact round-trip, plus `CWELD` cards. |
| [src/geom.cpp](src/geom.cpp) | Point-to-element projection. TRIA3: barycentric. QUAD4: inverse bilinear mapping, using damped Gauss-Newton that switches to Newton near the answer, with an exact edge fallback. Returns the node weights the connector ties to. |
| [src/bvh.cpp](src/bvh.cpp) | Bounding-volume hierarchy over element boxes (median split). Element boxes are stored in leaf order. |
| [src/connect.cpp](src/connect.cpp) | For each weld: closest element per part within a tolerance, then the 2–3 closest parts. A 3-sheet stack is ordered through the thickness, so it welds bottom–middle and middle–top, never bottom–top. Welds run in Z-curve order under OpenMP. |
| [src/topology.cpp](src/topology.cpp) | Node→element adjacency (CSR), connected components (union-find) before and after connectors, free and non-manifold edges, and coincident-node detection (spatial hash). |
| [src/c_api.cpp](src/c_api.cpp), [include/mc/c_api.h](include/mc/c_api.h) | C ABI of 10 functions. It is the only surface exported from `libmeshconnect.so` / `meshconnect.dll`, and no C++ types cross it. |
| [python/meshconnect.py](python/meshconnect.py) | Python binding over the C ABI with `ctypes`; no compiled extension. |
| [src/ccx.cpp](src/ccx.cpp), [validate/run_lap.py](validate/run_lap.py) | Lap joint exported to CalculiX `.inp`, with connectors as `*EQUATION` constraints, then solved and read back. |
| [src/gen.cpp](src/gen.cpp) | Test assemblies with ground truth (see below). |

## Test data and ground truth

There is no public OEM body-in-white model, so `mesh-connect gen` builds one
with known answers:

- **Panels:** 8 panels of 1000 × 500 mm in a row, alternating at z = 0 and z = 1 mm, joined by 40 mm lap flanges.
  - Element size changes per part (2.0 / 2.3 / 2.6 mm), so neighbouring meshes never line up.
  - Interior nodes are jittered by ±25% of the element size, so elements are distorted.
  - Every other part is meshed in triangles.
- **Patches:** reinforcement patches on 4 flanges make 3-sheet stacks.
- **Stiffeners:** stiffeners under the other 3 flanges sit 2.5 mm from the weld plane, outside the 1.2 mm tolerance. They must **not** be welded.
- **Weld labels:** 100,000 labelled welds:
  - 84,178 two-sheet,
  - 5,737 three-sheet,
  - 10,085 over a single sheet, which must be refused.

The sheets are flat, so the correct projection is known analytically:
(x, y, z_sheet) on each expected part. The labels therefore don't depend on the
code being tested.

Mesh: **1,203,735 elements** (422,315 CQUAD4 + 781,420 CTRIA3), 818,890 nodes,
15 parts, 90 MB of Nastran.

## Results

Machine: i5-9300H laptop, WSL2 with 4 vCPUs and 2.9 GB RAM, shared with an
IDE and a browser. Every number comes from a file in [results/](results/) and is
reproduced by [scripts/run_all.sh](scripts/run_all.sh). Where two runs
disagreed, the range is shown.

### Correctness on 100,000 labelled welds ([results/eval_1m.json](results/eval_1m.json))

| | Result |
|---|---|
| Welds with the right parts, right elements and right stack order | **100,000 / 100,000** |
| Largest error of a projected point vs the analytic answer | **1.8e-12 mm** |
| Single-sheet welds correctly refused | 10,085 / 10,085 |
| Stiffeners 2.5 mm away (outside tolerance) welded by mistake | 0; they are reported as the only floating parts |
| Baseline: "element with the nearest centroid" (a common shortcut) | **wrong element for 16.3%** of sheet hits (30,229 / 185,567) |

The centroid shortcut fails because, on a distorted or non-uniform mesh, the
nearest centroid often belongs to a neighbour of the element that actually
contains the point.

### Projection accuracy on warped quads ([results/projection_accuracy.json](results/projection_accuracy.json))

Each test uses a random 2 mm quad with opposite corners lifted ±warp × 2 mm.
The point is placed off the surface along the true normal at known (ξ, η), so
the answer is exact. There are 200,000 points per row.

| Warp | Newton: max position error | Split into 2 triangles: p99 position error | Split: max weight error |
|---|---|---|---|
| 0 (flat) | 1.4e-15 mm | 1.3e-15 mm | 0.36 |
| 0.05 | 1.6e-15 mm | 0.105 mm | 0.37 |
| 0.2 | 1.4e-15 mm | 0.37 mm | 0.36 |

Splitting a quad into two triangles is a common shortcut. It puts the point
in the right place only on flat quads, and even there the weights it ties to
are wrong by up to 0.36. The bilinear weights are what the element actually
interpolates.

### Speed: 1.2M elements, 100,000 welds ([results/bench_1m_run1.json](results/bench_1m_run1.json), [run2](results/bench_1m_run2.json))

| Measurement | Result |
|---|---|
| BVH vs brute force (1 thread) | 4.4–5.1 µs vs 102–170 ms per weld: **>20,000×**, with 0 different answers on a 200-weld sample |
| Thread scaling, 100k welds (median of 7) | 1 thread 0.46–0.47 s · 2 threads 0.23–0.27 s · 3 threads 0.19–0.20 s · 4 threads **0.16 s** (≈ 620k welds/s, 2.9×) |
| Same output at 1, 2, 3 and 4 threads | yes: identical hash of every element id, projected point and weight |
| Read the 90 MB deck | 1.1–1.7 s (60–85 MB/s) |
| BVH build (single-threaded) | 1.8–2.6 s; slower than all the queries together |
| Components / CSR / free edges / coincident nodes | 0.04 / 0.03 / 0.16–0.21 / 0.7–0.9 s ([check_1m.json](results/check_1m.json)) |

**What made it faster.** Single thread, interleaved A/B over 11 rounds, two
runs, median slowdown when each change is removed:

| Removed | Slowdown |
|---|---|
| Z-curve (Morton) ordering of the welds | 1.27–1.43× |
| Box lower-bound pruning (skip the projection if the box is farther than the part's best hit) | 1.16–1.18× |
| Both | 1.50–1.59× |

Both changes reduce cache misses. Per weld, BVH traversal alone takes about
2 µs; most of the rest went to reading the corner nodes of candidate elements.
Storing the element boxes in BVH leaf order means a candidate can be rejected
without touching node memory.

### Connectivity, before and after ([results/connect_1m.json](results/connect_1m.json))

```
components before: 15 (every part separate; no coincident nodes, 11,700 free edges, 0 non-manifold)
welds:             89,915 realised (5,737 three-sheet), 10,085 refused (one sheet), 95,652 CWELD cards
components after:  4  -> floating parts [201, 203, 205]  (the three stiffeners, by design)
```

### Do the connectors carry load? CalculiX lap joint ([results/lap_calculix.json](results/lap_calculix.json))

The model is two 100 × 40 × 1 mm steel plates overlapping by 20 mm, joined
by six welds placed off the nodes, with element sizes that don't match
(2.5 / 2.0 mm). Plate B is clamped and plate A is pulled with 1000 N. Each
connector style is solved at 1× to 4× mesh refinement.

| Connector | Reaction at clamp | Pulled-edge displacement vs reference, 1× → 2× → 3× → 4× mesh |
|---|---|---|
| None | **0 N; ccx still exits 0, displacement 5×10¹⁰ mm** | n/a |
| Nearest node of A tied to nearest node of B | 1000.0000 N | +4.5% → +5.6% → +6.6% → +7.1% |
| Point, tied with element shape functions (what `CWELD ELEMID` does) | 1000.0000 N | +1.5% → +2.9% → +3.4% → +3.9% |
| Patch: mean of all nodes within the 5 mm nugget, on each sheet | 1000.0000 N | +0.66% → +0.07% → +0.02% → reference |

- **Every connected model carries the full load.** The reaction equals the 1000 N applied to 4 decimal places.
- **Without connectors, the solver doesn't complain.** It returns exit code 0, "Job finished", and nonsense displacements. A component count catches a floating part before the solve.
- **A point connector is a singularity.** Tying a weld to one node, or to the 4 nodes of one element, gets softer as the mesh is refined and does not settle. Spreading the weld over its nugget diameter does settle (0.02% between the last two meshes). This is why area-type spot welds (ACM-style) exist.

## Bugs the checks found

1. **Newton tolerance below double precision** ([results/bug1_newton_tolerance_before_fix.json](results/bug1_newton_tolerance_before_fix.json)).
   - *Symptom:* the quad projection passed every unit test (coordinates near the origin) but got only **1,841 of 3,000** labelled welds right on the first assembly.
   - *Cause:* the stopping test was an absolute |Δξ| < 1e-14. At x ≈ 1000 mm, doubles cannot resolve a step that small, so Newton never "converged" and the code fell back to the element edges.
   - *Fix:* iterate in coordinates relative to the element and use a reachable tolerance.
   - *Test:* `quad_projection_far_from_origin`, which fails on the old code.
2. **Full Newton cycling on a distorted quad** ([results/bug2_newton_cycle_before_fix.json](results/bug2_newton_cycle_before_fix.json)).
   - *Symptom:* **2 of 100,000** welds landed on an edge 0.502 mm away instead of inside their element 0.500 mm away.
   - *Cause:* far from the answer the Newton Hessian is nearly indefinite. The step was clamped and cycled (2,2) → (2,−2) → (0.09,−0.11) → (2,2) for all 30 iterations.
   - *Fix:* Gauss-Newton with a backtracking line search, switching to full Newton only near the solution.
   - *Test:* `quad_projection_distorted_element` uses that exact element.
3. **`ENDDATA` counted as an unsupported card.** The checker reported "1 skipped card" on its own output. Fixed, with a test.

## Running it

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release && cmake --build build
ctest --test-dir build                          # 14 tests
python3 python/test_meshconnect.py              # ctypes -> libmeshconnect.so

build/mesh-connect gen work --parts 8 --h 2 --welds 100000
build/mesh-connect check   work/mesh.bdf
build/mesh-connect connect work/mesh.bdf work/welds.csv --out cwelds.bdf --report welds_report.csv
build/mesh-connect eval    work/mesh.bdf work/labels.csv
build/mesh-connect bench   work/mesh.bdf work/welds.csv
python3 validate/run_lap.py build/mesh-connect ccx work/lap   # needs CalculiX
scripts/run_all.sh build/mesh-connect /tmp/mc-work           # every number above
```

```python
from meshconnect import Mesh
m = Mesh("model.bdf")
welds = m.connect([(8.3, 2.7, 0.5)], tol=1.2)
print(welds[0], m.components(), m.components(welds))
```

CI ([.github/workflows/ci.yml](.github/workflows/ci.yml)) does the following:

- **Linux (GCC):** builds and runs the tests, the Python binding, the labelled-assembly check (must be 100%), BVH = brute force with identical output across thread counts, and the CalculiX lap joint.
- **Windows (MSVC, OpenMP 2.0):** builds the DLL and runs the tests and the Python binding against `meshconnect.dll`.

## Not done / limits

- **The test data is synthetic.** The assemblies are generated, with flat sheets and known answers; no real vehicle model was used. Curved surfaces are covered only by the warped-quad projection test.
- **Element and card support is narrow.** Shells only: CQUAD4 and CTRIA3, with no CQUAD8/CTRIA6 or solids. `GRID` with a local coordinate system (CP ≠ 0) is rejected rather than transformed. Parts are identified by PID, so a part folded onto itself (a hem) would not be welded to itself.
- **Only spot welds are handled.** There are no seam or adhesive lines and no bolts.
- **The CWELD output was never run through Nastran** (no licence). The FE check uses this tool's own CalculiX `*EQUATION` export of the same projection and weights.
- **The connectors tie translations only**, between mid-surfaces 1 mm apart, so they carry no moment from that offset. The lap joint shows no secondary bending (out-of-plane displacement ≈ 1e-10 mm), which a real lap joint has. Production CWELD/ACM elements span the gap with a beam or hexa for this reason.
- **The CalculiX check stops at 4× refinement.** At 8× (~106k shell nodes) the solver ran out of memory in the 2.9 GB VM.
- **Windows was only checked in CI.** It is built and tested only in CI; Windows App Control on the dev laptop blocks freshly built executables and DLLs. All timings are Linux (WSL2).
- **Thread counts stop at 4.** The machine has 4 vCPUs, so 8 threads were not measured.
- **The BVH build is single-threaded** and takes longer than querying all 100k welds.
