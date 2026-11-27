#!/usr/bin/env python3
"""Lap-joint check in CalculiX: do the connectors actually carry load?

Two 100 x 40 x 1 mm steel plates overlap by 20 mm and are joined by six spot
welds placed off the mesh nodes. Plate B is clamped at its far edge; plate A
is pulled with 1000 N. For each connector mode and mesh size we write the
.inp with `mesh-connect lap`, solve it with ccx, and read back:
  - the total reaction force at the clamp (must equal the 1000 N load),
  - the mean x-displacement of A's loaded edge (joint + plate stiffness).

usage: run_lap.py <mesh-connect binary> <ccx> <work dir> [out.json]
"""
import json
import os
import re
import subprocess
import sys


def solve(binary, ccx, work, name, mode, ha, hb):
    inp = os.path.join(work, name + ".inp")
    info = json.loads(subprocess.run([binary, "lap", inp, "--mode", mode, "--hA", str(ha), "--hB", str(hb)],
                                     check=True, capture_output=True, text=True).stdout)
    log = subprocess.run([ccx, "-i", name], cwd=work, capture_output=True, text=True)
    dat = open(os.path.join(work, name + ".dat")).read()
    rf = [float(v) for v in re.search(r"total force \(fx,fy,fz\) for set FIX.*?\n\s*\n\s*(\S+)\s+(\S+)\s+(\S+)", dat, re.S).groups()]
    disp = re.search(r"displacements \(vx,vy,vz\) for set LOAD.*?\n\s*\n(.*?)(?:\n\s*\n|\Z)", dat, re.S).group(1)
    rows = [[float(x) for x in line.split()[1:4]] for line in disp.strip().splitlines()]
    ux = sum(r[0] for r in rows) / len(rows)
    uz = max(abs(r[2]) for r in rows)
    return {
        "name": name, "mode": mode, "hA": ha, "hB": hb, "nodes": info["nodes"], "equations": info["equations"],
        "max_snap_mm": info["max_snap_mm"], "patch_nodes": info["patch_nodes"], "ccx_exit": log.returncode,
        "ccx_warned": bool(re.search(r"warning|error|singular", log.stdout + log.stderr, re.I)),
        "reaction_fx": rf[0], "edge_ux_mm": ux, "edge_uz_max_mm": uz,
    }


def main():
    binary, ccx, work = sys.argv[1:4]
    out = sys.argv[4] if len(sys.argv) > 4 else None
    os.makedirs(work, exist_ok=True)
    runs = [solve(binary, ccx, work, "none", "none", 2.5, 2.0)]
    scales = [int(x) for x in os.environ.get("LAP_SCALES", "1,2,4,8").split(",")]
    for scale in scales:
        ha, hb = 2.5 / scale, 2.0 / scale
        for mode in ("shape", "nearest", "patch"):
            runs.append(solve(binary, ccx, work, f"{mode}_{scale}", mode, ha, hb))
    ref = next(r for r in runs if r["name"] == f"patch_{scales[-1]}")["edge_ux_mm"]
    for r in runs:
        r["ux_vs_ref_pct"] = 100 * (r["edge_ux_mm"] - ref) / ref
        print(f'{r["name"]:10s} hA={r["hA"]:<6} hB={r["hB"]:<5} RFx={r["reaction_fx"]:12.4f}  '
              f'ux={r["edge_ux_mm"]: .6e}  ({r["ux_vs_ref_pct"]:+.2f}%)  snap={r["max_snap_mm"]:.3f}  '
              f'ccx_exit={r["ccx_exit"]} warned={r["ccx_warned"]}')
    ok = all(abs(r["reaction_fx"] - 1000) < 1e-6 for r in runs if r["mode"] != "none")
    floating = next(r for r in runs if r["mode"] == "none")
    result = {"load_N": 1000, "reference": f"patch_{scales[-1]}", "runs": runs, "all_connected_runs_carry_full_load": ok,
              "unconnected_run": {"ccx_exit": floating["ccx_exit"], "reaction_fx": floating["reaction_fx"],
                                  "edge_ux_mm": floating["edge_ux_mm"]}}
    if out:
        with open(out, "w") as f:
            json.dump(result, f, indent=1)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
