#!/usr/bin/env python3
"""
cut_stability_width.py — Python (uproot + matplotlib) port of CutStabilityWidth.cpp

What it does:
- Streams DATA TTree in chunks.
- For each variable (dx, dy, coin_time), fixes the midpoint and varies the *half-width*.
- Applies baseline cuts (vz, ePS, W2 window, eHCAL threshold, etc.) and then applies the
  symmetric window for the scanned variable:
      mid - width < var < mid + width
- Computes raw helicity asymmetry and uncertainty at each width.
- Writes CSVs and saves plots under ./plots/

This is designed to be a drop-in nicer-plot replacement for the ROOT macro.
"""

from __future__ import annotations

import argparse
import os
from dataclasses import dataclass
from typing import Tuple

import numpy as np

try:
    import uproot
except ImportError as e:
    raise SystemExit("This script needs uproot: pip install uproot awkward") from e

import matplotlib.pyplot as plt


@dataclass
class ScanResult:
    w: np.ndarray
    A: np.ndarray
    Aerr: np.ndarray
    Np: np.ndarray
    Nm: np.ndarray


def ensure_outdir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def write_csv(path: str, r: ScanResult) -> None:
    header = "width,A,Aerr,Nplus,Nminus\n"
    arr = np.column_stack([r.w, r.A, r.Aerr, r.Np.astype(float), r.Nm.astype(float)])
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
    if kin == "GEN2_He3":
        return 1
    if kin in ("GEN3_He3", "GEN4_He3", "GEN4b_He3"):
        return -1
    return -1


def stream_width_scan(
    filename: str,
    tree_name: str,
    var_name: str,
    widths: np.ndarray,
    mid: float,
    # Baseline cuts:
    W2L0: float, W2H0: float,
    eL0: float,
    tL0: float, tH0: float,
    dxL0: float, dxH0: float,
    dyL0: float, dyH0: float,
    kin: str,
) -> ScanResult:
    widths = np.asarray(widths, dtype=float)
    widths = np.sort(widths)

    Np = np.zeros_like(widths, dtype=np.int64)
    Nm = np.zeros_like(widths, dtype=np.int64)

    pkin = kin_sign(kin)

    branches = [
        "vz", "ePS", "eSH", "trP", "dx", "dy", "W2", "coin_time", "eHCAL",
        "helicity", "IHWP", var_name,
    ]
    # iterate
    for arrays in uproot.iterate(
        f"{filename}:{tree_name}",
        expressions=list(dict.fromkeys(branches)),  # unique
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
        v = arrays[var_name].astype(float)

        # preselection (same pattern as other macro)
        pre_reject = (
            (np.abs(vz) > 0.27)
            & (ePS < 0.2)
            & (np.abs((eSH + ePS) / trP - 1.0) > 0.2)
            & (eHCAL < 0.025)
            & (np.abs(helicity) != 1)
        )
        keep = ~pre_reject

        baseCommon = (
            (np.abs(vz) < 0.27)
            & (ePS > 0.2)
            & (W2 > W2L0) & (W2 < W2H0)
            & (eHCAL > eL0)
        )
        inCT = (ct > tL0) & (ct < tH0)
        inDX = (dx > dxL0) & (dx < dxH0)
        inDY = (dy > dyL0) & (dy < dyH0)

        keep &= baseCommon & inCT & inDX & inDY
        if not np.any(keep):
            continue

        h = pkin * IHWP * helicity
        pos = keep & (h == 1)
        neg = keep & (h == -1)

        # For symmetric window around mid: |v-mid| < width
        # We want counts for each width; do cumulative counting on |v-mid|.
        dv_pos = np.abs(v[pos] - mid)
        dv_neg = np.abs(v[neg] - mid)

        dv_pos.sort()
        dv_neg.sort()

        # Counts with dv < width:
        # side='left' gives number < width
        if dv_pos.size:
            Np += np.searchsorted(dv_pos, widths, side="left").astype(np.int64)
        if dv_neg.size:
            Nm += np.searchsorted(dv_neg, widths, side="left").astype(np.int64)

    A, Aerr = asym_from_counts(Np.astype(float), Nm.astype(float))
    return ScanResult(w=widths, A=A, Aerr=Aerr, Np=Np, Nm=Nm)


def plot_width(ax, r: ScanResult, title: str, xlabel: str) -> None:
    ax.errorbar(r.w, 100.0 * r.A, yerr=100.0 * r.Aerr, fmt="o", ms=4, capsize=2)
    ax.axhline(np.nanmedian(100.0 * r.A), ls="--", lw=1)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Raw asymmetry (%)")
    ax.grid(True, alpha=0.25)


def main() -> None:
    ap = argparse.ArgumentParser(description="Cut width stability (Python port of CutStabilityWidth.cpp)")
    ap.add_argument("data_root", help="DATA ROOT file")
    ap.add_argument("--tree", default="T", help="TTree name (default: T)")
    ap.add_argument("--kin", required=True, help="kin string (e.g., GEN3_He3)")
    ap.add_argument("--outdir", default="plots", help="Output directory")
    ap.add_argument("--fmt", default="pdf", choices=["pdf", "png", "both"])

    # Baseline cuts (override to match your exact analysis)
    ap.add_argument("--W2L0", type=float, default=-1.0)
    ap.add_argument("--W2H0", type=float, default=3.0)
    ap.add_argument("--eL0", type=float, default=0.025)
    ap.add_argument("--tL0", type=float, default=-50.0)
    ap.add_argument("--tH0", type=float, default=50.0)
    ap.add_argument("--dxL0", type=float, default=-3.0)
    ap.add_argument("--dxH0", type=float, default=3.0)
    ap.add_argument("--dyL0", type=float, default=-0.5)
    ap.add_argument("--dyH0", type=float, default=0.5)

    # Scan controls
    ap.add_argument("--nscan", type=int, default=25)
    ap.add_argument("--dx_mid", type=float, default=0.0)
    ap.add_argument("--dy_mid", type=float, default=0.0)
    ap.add_argument("--ct_mid", type=float, default=0.0)
    ap.add_argument("--dx_wmax", type=float, default=4.0)
    ap.add_argument("--dy_wmax", type=float, default=1.5)
    ap.add_argument("--ct_wmax", type=float, default=120.0)

    args = ap.parse_args()
    ensure_outdir(args.outdir)

    widths_dx = np.linspace(0.1, args.dx_wmax, args.nscan)
    widths_dy = np.linspace(0.05, args.dy_wmax, args.nscan)
    widths_ct = np.linspace(1.0, args.ct_wmax, args.nscan)

    # dx width scan
    print("\n=== dx width scan ===")
    r_dx = stream_width_scan(
        filename=args.data_root, tree_name=args.tree, var_name="dx",
        widths=widths_dx, mid=args.dx_mid,
        W2L0=args.W2L0, W2H0=args.W2H0,
        eL0=args.eL0,
        tL0=args.tL0, tH0=args.tH0,
        dxL0=args.dxL0, dxH0=args.dxH0,
        dyL0=args.dyL0, dyH0=args.dyH0,
        kin=args.kin,
    )
    write_csv(os.path.join(args.outdir, f"dx_width_{args.kin}.csv"), r_dx)

    # dy width scan
    print("\n=== dy width scan ===")
    r_dy = stream_width_scan(
        filename=args.data_root, tree_name=args.tree, var_name="dy",
        widths=widths_dy, mid=args.dy_mid,
        W2L0=args.W2L0, W2H0=args.W2H0,
        eL0=args.eL0,
        tL0=args.tL0, tH0=args.tH0,
        dxL0=args.dxL0, dxH0=args.dxH0,
        dyL0=args.dyL0, dyH0=args.dyH0,
        kin=args.kin,
    )
    write_csv(os.path.join(args.outdir, f"dy_width_{args.kin}.csv"), r_dy)

    # coin_time width scan
    print("\n=== coin_time width scan ===")
    r_ct = stream_width_scan(
        filename=args.data_root, tree_name=args.tree, var_name="coin_time",
        widths=widths_ct, mid=args.ct_mid,
        W2L0=args.W2L0, W2H0=args.W2H0,
        eL0=args.eL0,
        tL0=args.tL0, tH0=args.tH0,
        dxL0=args.dxL0, dxH0=args.dxH0,
        dyL0=args.dyL0, dyH0=args.dyH0,
        kin=args.kin,
    )
    write_csv(os.path.join(args.outdir, f"coin_time_width_{args.kin}.csv"), r_ct)

    # Plot
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8), constrained_layout=True)
    plot_width(axes[0], r_dx, f"dx width scan ({args.kin})", "half-width (dx)")
    plot_width(axes[1], r_dy, f"dy width scan ({args.kin})", "half-width (dy)")
    plot_width(axes[2], r_ct, f"coin_time width scan ({args.kin})", "half-width (coin_time)")

    outbase = os.path.join(args.outdir, f"CutStabilityWidth_{args.kin}")
    if args.fmt in ("pdf", "both"):
        fig.savefig(outbase + ".pdf")
        print(f"Saved {outbase}.pdf")
    if args.fmt in ("png", "both"):
        fig.savefig(outbase + ".png", dpi=200)
        print(f"Saved {outbase}.png")
    plt.close(fig)


if __name__ == "__main__":
    main()
