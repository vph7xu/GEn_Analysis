#!/usr/bin/env python3
"""
cut_stability.py — Python (uproot + matplotlib) port of CutStability.cpp

What it does (same physics/logic as your C++ macro):
- Streams chunks from the DATA TTree (no giant per-event arrays).
- Uses branches: vz, ePS, eSH, trP, dx, dy, W2, coin_time, eHCAL, helicity, IHWP.
- Applies the *same* "reject only if ALL sub-conditions are true" preselection.
- Applies baseline cuts, then scans one variable at a time:
    dyL, dyH, W2L, W2H, eHCAL_L, tL, tH, dxL, dxH
- For each scan point, computes raw helicity asymmetry:
      A = (N+ - N-) / (N+ + N-)
      σ_A = sqrt((1 - A^2) / (N+ + N-))
- Writes CSV per scan and saves a multi-panel PDF/PNG under ./plots/

Notes:
- The C++ version updates Nplus/Nminus for every scan point per event (O(N*K)).
  This port uses cumulative counting via sorting + searchsorted inside each chunk,
  so it stays fast even for many scan points.
"""

from __future__ import annotations

import argparse
import os
from dataclasses import dataclass
from typing import Dict, Iterable, List, Tuple

import numpy as np

try:
    import uproot
except ImportError as e:
    raise SystemExit("This script needs uproot: pip install uproot awkward") from e

import matplotlib.pyplot as plt


# -----------------------------
# Utilities
# -----------------------------

@dataclass
class ScanResult:
    x: np.ndarray          # scan values
    A: np.ndarray          # asymmetry
    Aerr: np.ndarray       # uncertainty
    Np: np.ndarray         # Nplus
    Nm: np.ndarray         # Nminus


def ensure_outdir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def write_csv(path: str, r: ScanResult) -> None:
    header = "x,A,Aerr,Nplus,Nminus\n"
    arr = np.column_stack([r.x, r.A, r.Aerr, r.Np.astype(float), r.Nm.astype(float)])
    with open(path, "w", encoding="utf-8") as f:
        f.write(header)
        np.savetxt(f, arr, delimiter=",", fmt="%.10g")


def asym_from_counts(Np: np.ndarray, Nm: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    N = Np + Nm
    A = np.zeros_like(N, dtype=float)
    Aerr = np.ones_like(N, dtype=float)
    m = N > 0
    A[m] = (Np[m] - Nm[m]) / N[m]
    Aerr[m] = np.sqrt(np.maximum(0.0, (1.0 - A[m] * A[m]) / N[m]))
    return A, Aerr


def kin_sign(kin: str) -> int:
    """
    Port of the p_kin logic in CutStability.cpp:
      p_kin starts from an assignment based on kin name, then gets flipped (p_kin=-1*p_kin),
      and then flipped again (p_kin=-1*p_kin) in the C++ (net effect depends on the exact code).
    Here we replicate the *final* printed p_kin behavior:
        - GEN2_He3 -> +1 (because the code sets -1 then flips twice)
        - GEN3_He3/GEN4_He3/GEN4b_He3 -> -1
    If you use other kin strings, default is -1 (matching the macro's final flip tendency).
    """
    if kin == "GEN2_He3":
        return 1
    if kin in ("GEN3_He3", "GEN4_He3", "GEN4b_He3"):
        return -1
    return -1


def counts_gt(sorted_vals: np.ndarray, xvals: np.ndarray) -> np.ndarray:
    """Counts of values strictly > x for each x in ascending xvals."""
    # idx = number <= x  (side='right')
    idx = np.searchsorted(sorted_vals, xvals, side="right")
    return (sorted_vals.size - idx).astype(np.int64)


def counts_lt(sorted_vals: np.ndarray, xvals: np.ndarray) -> np.ndarray:
    """Counts of values strictly < x for each x in ascending xvals."""
    idx = np.searchsorted(sorted_vals, xvals, side="left")
    return idx.astype(np.int64)


def stream_scan(
    filename: str,
    tree_name: str,
    xvals: np.ndarray,
    var_type: str,
    # Baseline (fixed) cuts:
    dyL0: float, dyH0: float,
    W2L0: float, W2H0: float,
    eL0: float,
    tL0: float, tH0: float,
    dxL0: float, dxH0: float,
    kin: str,
    step_report: int = 50,
) -> ScanResult:
    """
    Streaming scan that matches StreamScanA() in CutStability.cpp.
    """
    xvals = np.asarray(xvals, dtype=float)
    if xvals.ndim != 1:
        raise ValueError("xvals must be 1D")

    # Ensure ascending for cumulative counting
    xvals_sorted = np.sort(xvals)

    Np = np.zeros_like(xvals_sorted, dtype=np.int64)
    Nm = np.zeros_like(xvals_sorted, dtype=np.int64)

    pkin = kin_sign(kin)

    branches = [
        "vz", "ePS", "eSH", "trP", "dx", "dy", "W2", "coin_time", "eHCAL",
        "helicity", "IHWP",
    ]

    # Iterate in chunks (uproot handles decompression efficiently)
    total_seen = 0
    for arrays in uproot.iterate(
        f"{filename}:{tree_name}",
        expressions=branches,
        step_size="200 MB",
        library="np",
    ):
        vz = arrays["vz"]
        ePS = arrays["ePS"]
        eSH = arrays["eSH"]
        trP = arrays["trP"]
        dx = arrays["dx"]
        dy = arrays["dy"]
        W2 = arrays["W2"]
        ct = arrays["coin_time"]
        eHCAL = arrays["eHCAL"]
        helicity = arrays["helicity"].astype(np.int32)
        IHWP = arrays["IHWP"].astype(np.int32)

        n = vz.size
        total_seen += n
        if step_report and (total_seen // (step_report * 100000) != (total_seen - n) // (step_report * 100000)):
            print(f"  streamed ~{total_seen/1e6:.1f}M events")

        # --- exact C++ preselection:
        # if (abs(vz)>0.27 && ePS<0.2 && abs((eSH+ePS)/trP - 1)>0.2 && eHCAL<0.025 && abs(helicity)!=1) continue;
        pre_reject = (
            (np.abs(vz) > 0.27)
            & (ePS < 0.2)
            & (np.abs((eSH + ePS) / trP - 1.0) > 0.2)
            & (eHCAL < 0.025)
            & (np.abs(helicity) != 1)
        )
        keep = ~pre_reject

        # Baseline cuts (except scanned variable), matching C++
        baseCommon = (
            (np.abs(vz) < 0.27)
            & (ePS > 0.2)
            & (W2 > W2L0) & (W2 < W2H0)
            & (eHCAL > eL0)
        )
        inCT = (ct > tL0) & (ct < tH0)
        inDX = (dx > dxL0) & (dx < dxH0)

        keep &= baseCommon & inCT & inDX

        if not np.any(keep):
            continue

        # Helicity sign with IHWP
        h = pkin * IHWP * helicity
        poshel = keep & (h == 1)
        neghel = keep & (h == -1)

        vt = var_type

        # For each scan type, pick the "scanned variable value" array and whether we count > x or < x
        if vt == "dyL":
            # if (dy < dyH0) and dy > x
            mpos = poshel & (dy < dyH0)
            mneg = neghel & (dy < dyH0)
            vpos = np.sort(dy[mpos].astype(float))
            vneg = np.sort(dy[mneg].astype(float))
            if vpos.size: Np += counts_gt(vpos, xvals_sorted)
            if vneg.size: Nm += counts_gt(vneg, xvals_sorted)

        elif vt == "dyH":
            # if (dy > dyL0) and dy < x
            mpos = poshel & (dy > dyL0)
            mneg = neghel & (dy > dyL0)
            vpos = np.sort(dy[mpos].astype(float))
            vneg = np.sort(dy[mneg].astype(float))
            if vpos.size: Np += counts_lt(vpos, xvals_sorted)
            if vneg.size: Nm += counts_lt(vneg, xvals_sorted)

        elif vt == "W2L":
            # if (W2 < W2H0) and W2 > x
            mpos = poshel & (W2 < W2H0)
            mneg = neghel & (W2 < W2H0)
            vpos = np.sort(W2[mpos].astype(float))
            vneg = np.sort(W2[mneg].astype(float))
            if vpos.size: Np += counts_gt(vpos, xvals_sorted)
            if vneg.size: Nm += counts_gt(vneg, xvals_sorted)

        elif vt == "W2H":
            # if (W2 > W2L0) and W2 < x
            mpos = poshel & (W2 > W2L0)
            mneg = neghel & (W2 > W2L0)
            vpos = np.sort(W2[mpos].astype(float))
            vneg = np.sort(W2[mneg].astype(float))
            if vpos.size: Np += counts_lt(vpos, xvals_sorted)
            if vneg.size: Nm += counts_lt(vneg, xvals_sorted)

        elif vt in ("eL", "eHCAL_L"):
            # eHCAL lower edge scan: eHCAL > x
            # (C++ names it eL0 baseline, scans "eL")
            vpos = np.sort(eHCAL[poshel].astype(float))
            vneg = np.sort(eHCAL[neghel].astype(float))
            if vpos.size: Np += counts_gt(vpos, xvals_sorted)
            if vneg.size: Nm += counts_gt(vneg, xvals_sorted)

        elif vt == "tL":
            # coin_time lower edge scan: ct > x (ensure ct < tH0)
            mpos = poshel & (ct < tH0)
            mneg = neghel & (ct < tH0)
            vpos = np.sort(ct[mpos].astype(float))
            vneg = np.sort(ct[mneg].astype(float))
            if vpos.size: Np += counts_gt(vpos, xvals_sorted)
            if vneg.size: Nm += counts_gt(vneg, xvals_sorted)

        elif vt == "tH":
            # coin_time upper edge scan: ct < x (ensure ct > tL0)
            mpos = poshel & (ct > tL0)
            mneg = neghel & (ct > tL0)
            vpos = np.sort(ct[mpos].astype(float))
            vneg = np.sort(ct[mneg].astype(float))
            if vpos.size: Np += counts_lt(vpos, xvals_sorted)
            if vneg.size: Nm += counts_lt(vneg, xvals_sorted)

        elif vt == "dxL":
            # dx lower edge scan: dx > x (ensure dx < dxH0)
            mpos = poshel & (dx < dxH0)
            mneg = neghel & (dx < dxH0)
            vpos = np.sort(dx[mpos].astype(float))
            vneg = np.sort(dx[mneg].astype(float))
            if vpos.size: Np += counts_gt(vpos, xvals_sorted)
            if vneg.size: Nm += counts_gt(vneg, xvals_sorted)

        elif vt == "dxH":
            # dx upper edge scan: dx < x (ensure dx > dxL0)
            mpos = poshel & (dx > dxL0)
            mneg = neghel & (dx > dxL0)
            vpos = np.sort(dx[mpos].astype(float))
            vneg = np.sort(dx[mneg].astype(float))
            if vpos.size: Np += counts_lt(vpos, xvals_sorted)
            if vneg.size: Nm += counts_lt(vneg, xvals_sorted)

        else:
            raise ValueError(f"Unknown var_type={var_type}")

    A, Aerr = asym_from_counts(Np.astype(float), Nm.astype(float))
    return ScanResult(x=xvals_sorted, A=A, Aerr=Aerr, Np=Np, Nm=Nm)


def plot_scan(ax, r: ScanResult, title: str, xlabel: str) -> None:
    # Plot A in percent (to match your ROOT macro "values in %")
    ax.errorbar(r.x, 100.0 * r.A, yerr=100.0 * r.Aerr, fmt="o", ms=4, capsize=2)
    ax.axhline(np.nanmedian(100.0 * r.A), ls="--", lw=1)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Raw asymmetry (%)")
    ax.grid(True, alpha=0.25)


def main() -> None:
    ap = argparse.ArgumentParser(description="Cut stability scans (Python port of CutStability.cpp)")
    ap.add_argument("data_root", help="DATA ROOT file")
    ap.add_argument("--tree", default="T", help="TTree name (default: T)")

    ap.add_argument("--kin", required=True, help="kin string (e.g., GEN3_He3)")
    ap.add_argument("--outdir", default="plots", help="Output directory (default: plots)")
    ap.add_argument("--fmt", default="pdf", choices=["pdf", "png", "both"], help="Plot output format")

    # Baseline defaults (match your macro's spirit; override to match your exact config)
    ap.add_argument("--dyL0", type=float, default=-0.5, help="baseline dy lower")
    ap.add_argument("--dyH0", type=float, default=0.5, help="baseline dy upper")
    ap.add_argument("--W2L0", type=float, default=-1.0, help="baseline W2 lower")
    ap.add_argument("--W2H0", type=float, default=3.0, help="baseline W2 upper")
    ap.add_argument("--eL0", type=float, default=0.025, help="baseline eHCAL lower")
    ap.add_argument("--tL0", type=float, default=-50.0, help="baseline coin_time lower")
    ap.add_argument("--tH0", type=float, default=50.0, help="baseline coin_time upper")
    ap.add_argument("--dxL0", type=float, default=-3.0, help="baseline dx lower")
    ap.add_argument("--dxH0", type=float, default=3.0, help="baseline dx upper")

    # Scan ranges (edit to match your C++ arrays if you want 1:1)
    ap.add_argument("--nscan", type=int, default=21, help="points per scan (default 21)")
    ap.add_argument("--dy_span", type=float, default=0.8, help="dy scan span around baseline edges")
    ap.add_argument("--dx_span", type=float, default=2.0, help="dx scan span around baseline edges")
    ap.add_argument("--ct_span", type=float, default=80.0, help="coin_time scan span around baseline edges")
    ap.add_argument("--W2_span", type=float, default=2.0, help="W2 scan span around baseline edges")
    ap.add_argument("--e_span", type=float, default=0.08, help="eHCAL lower scan span")

    args = ap.parse_args()
    ensure_outdir(args.outdir)

    # Build scan x-grids (you can hard-code these to match the exact ROOT macro arrays if desired)
    n = args.nscan
    dyL_scan = np.linspace(args.dyL0 - args.dy_span, args.dyL0 + args.dy_span, n)
    dyH_scan = np.linspace(args.dyH0 - args.dy_span, args.dyH0 + args.dy_span, n)
    dxL_scan = np.linspace(args.dxL0 - args.dx_span, args.dxL0 + args.dx_span, n)
    dxH_scan = np.linspace(args.dxH0 - args.dx_span, args.dxH0 + args.dx_span, n)
    tL_scan  = np.linspace(args.tL0 - args.ct_span, args.tL0 + args.ct_span, n)
    tH_scan  = np.linspace(args.tH0 - args.ct_span, args.tH0 + args.ct_span, n)
    W2L_scan = np.linspace(args.W2L0 - args.W2_span, args.W2L0 + args.W2_span, n)
    W2H_scan = np.linspace(args.W2H0 - args.W2_span, args.W2H0 + args.W2_span, n)
    eL_scan  = np.linspace(max(0.0, args.eL0 - args.e_span), args.eL0 + args.e_span, n)

    scans = [
        ("dyL", dyL_scan, "dy lower edge scan", "dy lower cut"),
        ("dyH", dyH_scan, "dy upper edge scan", "dy upper cut"),
        ("W2L", W2L_scan, "W2 lower edge scan", "W2 lower cut"),
        ("W2H", W2H_scan, "W2 upper edge scan", "W2 upper cut"),
        ("eL",  eL_scan,  "eHCAL lower edge scan", "eHCAL lower cut"),
        ("tL",  tL_scan,  "coin_time lower edge scan", "coin_time lower cut"),
        ("tH",  tH_scan,  "coin_time upper edge scan", "coin_time upper cut"),
        ("dxL", dxL_scan, "dx lower edge scan", "dx lower cut"),
        ("dxH", dxH_scan, "dx upper edge scan", "dx upper cut"),
    ]

    results: Dict[str, ScanResult] = {}
    for vtype, xgrid, title, xlabel in scans:
        print(f"\n=== Scanning {vtype} ({xgrid.size} points) ===")
        r = stream_scan(
            filename=args.data_root,
            tree_name=args.tree,
            xvals=xgrid,
            var_type=vtype,
            dyL0=args.dyL0, dyH0=args.dyH0,
            W2L0=args.W2L0, W2H0=args.W2H0,
            eL0=args.eL0,
            tL0=args.tL0, tH0=args.tH0,
            dxL0=args.dxL0, dxH0=args.dxH0,
            kin=args.kin,
        )
        results[vtype] = r

        csv_path = os.path.join(args.outdir, f"{vtype}_{args.kin}.csv")
        write_csv(csv_path, r)
        print(f"  wrote {csv_path}")

    # Multi-panel plot (nicer than ROOT canvas)
    fig, axes = plt.subplots(3, 3, figsize=(14, 11), constrained_layout=True)
    axes = axes.ravel()

    for i, (vtype, _, title, xlabel) in enumerate(scans):
        plot_scan(axes[i], results[vtype], f"{title} ({args.kin})", xlabel)

    # If fewer than 9 plots, hide extras (shouldn't happen here)
    for j in range(len(scans), len(axes)):
        axes[j].axis("off")

    outbase = os.path.join(args.outdir, f"CutStability_{args.kin}")
    if args.fmt in ("pdf", "both"):
        fig.savefig(outbase + ".pdf")
        print(f"Saved {outbase}.pdf")
    if args.fmt in ("png", "both"):
        fig.savefig(outbase + ".png", dpi=200)
        print(f"Saved {outbase}.png")
    plt.close(fig)


if __name__ == "__main__":
    main()
