#!/usr/bin/env python3
"""
GEN_plotter.py

Python version of your ROOT macro behavior:
- World data from DB/World_Data.dat (grouped by experiment name)
- Ye global-fit central + 1σ band for BOTH ratio and GE from neutron_lookup.dat
- RCQM (Miller) curve from DB/F1u.csv, F1d.csv, F2u.csv, F2d.csv (same math as C++)
- DSE (Roberts) curve from PSM_theory() (same function)
- "This Work" points with stat error bars + systematic boxes
- Watermark "Exploratory"

Outputs:
  <outdir>/kin23Ratio.png, .pdf
  <outdir>/kin23GEn.png,   .pdf
"""

from __future__ import annotations

import argparse
import os
import math
from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch


# ----------------------------- constants -------------------------------------
MU_N = -1.9130427  # neutron magnetic moment (same sign convention as your code uses in places)
KAPPA_D = -2.03
KAPPA_U =  1.67


# ----------------------- THIS WORK (edit if needed) --------------------------
# Q2 points:
Q2_THIS = np.array([2.99, 6.82, 9.78])

# Panel 1 (ratio): mu_n * GE/GM  (your spreadsheet row "mu_n*GEn/G_Mn")
R_THIS      = np.array([0.4677, 0.7886, 0.6936])
R_THIS_STAT = np.array([0.0387, 0.1311, 0.2508])
R_THIS_SYS  = np.array([0.0345, 0.0342, 0.0722])

# Panel 2 (GE): GE^n
GE_THIS      = np.array([0.0172, 0.0059, 0.0023])
GE_THIS_STAT = np.array([0.0014, 0.0010, 0.0008])
GE_THIS_SYS  = np.array([0.0013, 0.0003, 0.0003])
# -----------------------------------------------------------------------------


# ----------------------------- helpers ---------------------------------------
def GD(Q2: np.ndarray) -> np.ndarray:
    return (1.0 + Q2 / 0.71) ** (-2)


def set_style():
    """Apply a consistent style to every figure.

    The legend font size can be overridden by passing a value
    through the command line (``--legend-fontsize``).  We keep
    the default at 26 for backwards compatibility but users who
    want smaller legends can now easily dial it down.
    """

    plt.rcParams.update({
        "figure.dpi": 140,
        "savefig.dpi": 300,
        "font.size": 22,
        "axes.labelsize": 34,
        "axes.titlesize": 40,
        # legend font size is controlled by callers via an argument
        "legend.fontsize": set_style.legend_fontsize,
        "xtick.labelsize": 28,
        "ytick.labelsize": 28,
        "axes.spines.top": True,
        "axes.spines.right": True,
    })

# default value for the legend font size; callers may override this
# before invoking ``set_style`` (see ``main``).
set_style.legend_fontsize = 20


def add_watermark(ax, text="Exploratory", fontsize: float = 110):
    """Draw a semi‑transparent watermark on an axes.

    ``fontsize`` defaults to the old hard‑coded value (110) but can
    be overridden from the command line via the ``--watermark-font``
    option.
    """

    ax.text(
        0.55, 0.45, text,
        transform=ax.transAxes,
        ha="center", va="center",
        fontsize=fontsize,
        color="red",
        alpha=0.25,
        fontweight="bold",
        zorder=0
    )


def draw_sys_boxes(ax, x, y, ysys, dx=0.32, alpha=0.25, zorder=1):
    for xi, yi, si in zip(np.asarray(x), np.asarray(y), np.asarray(ysys)):
        if not np.isfinite(si) or si <= 0:
            continue
        rect = Rectangle(
            (xi - dx, yi - si), 2*dx, 2*si,
            facecolor="red", edgecolor="red",
            alpha=alpha, linewidth=1.0, zorder=zorder
        )
        ax.add_patch(rect)


# ----------------------------- World data ------------------------------------
@dataclass
class WorldGroup:
    name: str
    q2: np.ndarray
    ratio: np.ndarray
    ratio_err: np.ndarray
    ge: np.ndarray
    ge_err: np.ndarray


def load_world_data_dat(path: str) -> List[WorldGroup]:
    """
    Matches your C++ World_Data.dat reading:
    Each non-comment line should look like:
      Q2  NAME  (mu_n GE/GM)  err_ratio  GE  err_GE
    but we robustly skip any line whose first token isn't numeric
    (this fixes your 'end' parsing crash).
    """
    groups: Dict[str, Dict[str, List[float]]] = {}

    with open(path, "r") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            if len(parts) < 6:
                continue

            # robust numeric check for Q2
            try:
                q2 = float(parts[0])
            except ValueError:
                continue  # e.g. "end"

            name = parts[1]
            try:
                ratio = float(parts[2])
                ratio_err = float(parts[3])
                ge = float(parts[4])
                ge_err = float(parts[5])
            except ValueError:
                continue

            if name not in groups:
                groups[name] = {"q2": [], "ratio": [], "ratio_err": [], "ge": [], "ge_err": []}

            groups[name]["q2"].append(q2)
            groups[name]["ratio"].append(ratio)
            groups[name]["ratio_err"].append(ratio_err)
            groups[name]["ge"].append(ge)
            groups[name]["ge_err"].append(ge_err)

    out: List[WorldGroup] = []
    for name, d in groups.items():
        out.append(WorldGroup(
            name=name,
            q2=np.array(d["q2"], dtype=float),
            ratio=np.array(d["ratio"], dtype=float),
            ratio_err=np.array(d["ratio_err"], dtype=float),
            ge=np.array(d["ge"], dtype=float),
            ge_err=np.array(d["ge_err"], dtype=float),
        ))
    return out


def world_marker_style(name: str):
    """
    Make it look like your ROOT plot (common ones).
    Fallback cycles if unknown.
    """
    # ROOT-like mapping from your legend image
    if "Madey" in name:
        return dict(marker="s", mfc="black", mec="black", color="black")
    if "Geis" in name:
        return dict(marker="^", mfc="limegreen", mec="limegreen", color="limegreen")
    if "Riordan" in name:
        return dict(marker="v", mfc="blue", mec="blue", color="blue")
    if "Schlimme" in name:
        return dict(marker="o", mfc="none", mec="magenta", color="magenta")
    return dict(marker="o", mfc="none", mec="black", color="black")


# ----------------------------- Ye band (FROM neutron_lookup.dat) -------------
def load_neutron_lookup(path: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    neutron_lookup.dat columns (from your file header):
      Q2, GEn/GD, dGEn/GD, dGEn_Par/GD, GMn/muGD, dGMn/muGD, dGMn_Par/muGD

    We build:
      GE(Q2) and sigma_GE
      ratio(Q2) = (mu_n*GE/GMn) and sigma_ratio
    using the same implied Ye-style fit content.

    Key identities:
      GE = (GEn/GD) * GD
      GMn = (GMn/muGD) * (mu_n) * GD
      ratio = mu_n*GE/GMn = (GEn/GD) / (GMn/muGD)

    For uncertainties we use the provided TOTAL columns:
      d(GEn/GD) = dGEn/GD
      d(GMn/muGD) = dGMn/muGD
    and propagate assuming uncorrelated.
    """
    rows = []
    with open(path, "r") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            if len(parts) < 7:
                continue
            try:
                Q2 = float(parts[0])
                A = float(parts[1])   # GEn/GD
                sA = float(parts[2])  # dGEn/GD
                B = float(parts[4])   # GMn/muGD
                sB = float(parts[5])  # dGMn/muGD
            except ValueError:
                continue
            rows.append((Q2, A, sA, B, sB))

    arr = np.array(rows, dtype=float)
    Q2 = arr[:, 0]
    A  = arr[:, 1]
    sA = arr[:, 2]
    B  = arr[:, 3]
    sB = arr[:, 4]

    gd = GD(Q2)

    GE = A * gd
    sGE = sA * gd

    ratio = np.zeros_like(Q2)
    sratio = np.zeros_like(Q2)

    good = (np.abs(B) > 0) & (np.abs(A) > 0)
    ratio[good] = A[good] / B[good]
    # sigma(r) for r=A/B
    sratio[good] = np.abs(ratio[good]) * np.sqrt((sA[good] / A[good])**2 + (sB[good] / B[good])**2)

    return Q2, GE, sGE, ratio, sratio, gd


# ----------------------------- RCQM curve (FROM F1/F2 CSVs) -------------------
def load_two_col_csv(path: str) -> Tuple[np.ndarray, np.ndarray]:
    xs, ys = [], []
    with open(path, "r") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = [p.strip() for p in s.split(",")]
            if len(parts) < 2:
                continue
            try:
                xs.append(float(parts[0]))
                ys.append(float(parts[1]))
            except ValueError:
                continue
    return np.array(xs, dtype=float), np.array(ys, dtype=float)


def interp1(x: np.ndarray, y: np.ndarray, xq: np.ndarray) -> np.ndarray:
    # safe linear interpolation
    return np.interp(xq, x, y)


def f1f2_to_ratio_or_ge(F1: np.ndarray, F2: np.ndarray, Q2: np.ndarray, want_ratio: bool) -> np.ndarray:
    # This matches your C++ helper F1F2_to_GEGM() logic
    Mn = 0.939565  # GeV
    tau = Q2 / (4.0 * Mn * Mn)
    GE = F1 - tau * F2
    GM = F1 + F2
    if want_ratio:
        # mu_n * GE / GM
        out = np.zeros_like(Q2)
        good = np.abs(GM) > 0
        out[good] = MU_N * GE[good] / GM[good]
        return out
    else:
        return GE  # (your C++ returns GE for the GE panel)


def rcqm_curve(dbdir: str, want_ratio: bool) -> Tuple[np.ndarray, np.ndarray]:
    xF1u, yF1u = load_two_col_csv(os.path.join(dbdir, "F1u.csv"))
    xF1d, yF1d = load_two_col_csv(os.path.join(dbdir, "F1d.csv"))
    xF2u, yF2u = load_two_col_csv(os.path.join(dbdir, "F2u.csv"))
    xF2d, yF2d = load_two_col_csv(os.path.join(dbdir, "F2d.csv"))

    # Q2 grid similar to your macro
    Q2max = 12.0
    nQ2 = 30
    Q2_plot = np.array([0.0] + [i*(Q2max/nQ2) + 0.15 for i in range(nQ2+1)], dtype=float)

    # Eval like ROOT TGraph::Eval
    # and apply the same /Q^4 and kappas
    eps = 1e-12
    Q2_safe = np.maximum(Q2_plot, eps)

    F1u = interp1(xF1u, yF1u, Q2_plot) / (Q2_safe**4)
    F1d = interp1(xF1d, yF1d, Q2_plot) / (Q2_safe**4)
    F2u = interp1(xF2u, yF2u, Q2_plot) * KAPPA_U / (Q2_safe**4)
    F2d = interp1(xF2d, yF2d, Q2_plot) * KAPPA_D / (Q2_safe**4)

    # your C++ combination:
    F1 = 2.0/3.0 * F1d - 1.0/3.0 * F1u
    F2 = 2.0/3.0 * F2d - 1.0/3.0 * F2u

    y = f1f2_to_ratio_or_ge(F1, F2, Q2_plot, want_ratio=want_ratio)
    return Q2_plot, y


# ----------------------------- DSE curve -------------------------------------
def PSM_theory(Q2: np.ndarray) -> np.ndarray:
    # exact same as your TF1 in C++
    x = Q2
    num = (0.4133*x - 0.1038*x**2 + 0.01252*x**3 + 0.007748*x**4)
    den = (1 + 0.9771*x - 0.3506*x**2 + 0.09863*x**3 + 0.009927*x**4 + 0.000002173*x**5)
    return num / den


# ----------------------------- plotting --------------------------------------
def rooty_ticks(ax):
    ax.minorticks_on()
    ax.tick_params(axis="both", which="major", direction="in", length=10)
    ax.tick_params(axis="both", which="minor", direction="in", length=6)


def plot_ratio(dbdir: str, world_dat: str, neutron_lookup: str, outdir: str,
               watermark: str, watermark_font: float,
               world_marker_size: float, this_marker_size: float):
    groups = load_world_data_dat(world_dat)
    Q2ye, _, _, Rye, sRye, _ = load_neutron_lookup(neutron_lookup)

    # RCQM + DSE
    q_rcqm, y_rcqm = rcqm_curve(dbdir, want_ratio=True)
    q_grid = np.linspace(1e-3, 12.0, 400)
    y_dse = PSM_theory(q_grid)

    fig, ax = plt.subplots(figsize=(12.5, 8.8))
    ax.set_title(r"$G_E/G_M$ Neutron Results", pad=18)
    ax.set_xlabel(r"$Q^2\;(\mathrm{GeV}^2)$")
    ax.set_ylabel(r"$\mu_n\,G_E^n/G_M^n$")
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 1.8)

    # Ye band + central (from neutron_lookup.dat)
    ax.fill_between(Q2ye, Rye - sRye, Rye + sRye, alpha=0.25, linewidth=0, label=r"Global Fit - Ye (1$\sigma$)")
    ax.plot(Q2ye, Rye, linewidth=2.5, label="Global Fit - Ye (central)")

    # RCQM and DSE
    ax.plot(q_rcqm, y_rcqm, linestyle=":", linewidth=2.8, label="RCQM - Miller")
    ax.plot(q_grid, y_dse, linestyle="--", linewidth=2.6, label="DSE - Roberts")

    # World points
    for g in groups:
        st = world_marker_style(g.name)
        ax.errorbar(
            g.q2, g.ratio, yerr=g.ratio_err,
            fmt=st["marker"],
            markersize=world_marker_size,
            markerfacecolor=st["mfc"],
            markeredgecolor=st["mec"],
            ecolor=st["color"],
            elinewidth=1.8,
            capsize=0,
            linestyle="none",
            label=g.name
        )

    # This work syst boxes + stat bars
    draw_sys_boxes(ax, Q2_THIS, R_THIS, R_THIS_SYS, dx=0.32, alpha=0.25, zorder=1)
    ax.errorbar(
        Q2_THIS, R_THIS, yerr=R_THIS_STAT,
        fmt="o", markersize=this_marker_size,
        markerfacecolor="red", markeredgecolor="red",
        ecolor="red", elinewidth=2.0, capsize=0,
        linestyle="none",
        label="This Work",
        zorder=3
    )

    if watermark:
        add_watermark(ax, watermark, fontsize=watermark_font)

    rooty_ticks(ax)

    # legend: add syst patch explicitly
    sys_patch = Patch(facecolor="red", edgecolor="red", alpha=0.25, label="This Work (syst.)")
    handles, labels = ax.get_legend_handles_labels()
    if "This Work (syst.)" not in labels:
        handles.append(sys_patch)
        labels.append("This Work (syst.)")

    ax.legend(handles, labels, frameon=False, loc="upper left", ncol=2)

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "kin234Ratio.png"))
    fig.savefig(os.path.join(outdir, "kin234Ratio.pdf"))
    plt.close(fig)


def plot_ge(world_dat: str, neutron_lookup: str, outdir: str,
            watermark: str, watermark_font: float,
            world_marker_size: float, this_marker_size: float):
    groups = load_world_data_dat(world_dat)
    Q2ye, GEye, sGEye, _, _, _ = load_neutron_lookup(neutron_lookup)

    fig, ax = plt.subplots(figsize=(12.5, 8.8))
    ax.set_title(r"$G_E^n$ Neutron Results", pad=18)
    ax.set_xlabel(r"$Q^2\;(\mathrm{GeV}^2)$")
    ax.set_ylabel(r"$G_E^n$")
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 0.08)

    # Ye band + central for GE (THIS is what you were missing)
    ax.fill_between(Q2ye, GEye - sGEye, GEye + sGEye, alpha=0.25, linewidth=0, label=r"Global Fit - Ye (1$\sigma$)")
    ax.plot(Q2ye, GEye, linewidth=2.5, label="Global Fit - Ye (central)")

    # World points (GE columns)
    for g in groups:
        st = world_marker_style(g.name)
        ax.errorbar(
            g.q2, g.ge, yerr=g.ge_err,
            fmt=st["marker"],
            markersize=world_marker_size,
            markerfacecolor=st["mfc"],
            markeredgecolor=st["mec"],
            ecolor=st["color"],
            elinewidth=1.8,
            capsize=0,
            linestyle="none",
            label=g.name
        )

    # This work syst boxes + stat bars
    draw_sys_boxes(ax, Q2_THIS, GE_THIS, GE_THIS_SYS, dx=0.32, alpha=0.25, zorder=1)
    ax.errorbar(
        Q2_THIS, GE_THIS, yerr=GE_THIS_STAT,
        fmt="o", markersize=this_marker_size,
        markerfacecolor="red", markeredgecolor="red",
        ecolor="red", elinewidth=2.0, capsize=0,
        linestyle="none",
        label="This Work",
        zorder=3
    )

    if watermark:
        add_watermark(ax, watermark, fontsize=watermark_font)

    rooty_ticks(ax)

    sys_patch = Patch(facecolor="red", edgecolor="red", alpha=0.25, label="This Work (syst.)")
    handles, labels = ax.get_legend_handles_labels()
    if "This Work (syst.)" not in labels:
        handles.append(sys_patch)
        labels.append("This Work (syst.)")

    ax.legend(handles, labels, frameon=False, loc="upper left", ncol=2)

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "kin234GEn.png"))
    fig.savefig(os.path.join(outdir, "kin234GEn.pdf"))
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dbdir", default="DB", help="DB directory (default: DB)")
    ap.add_argument("--world-dat", default=None, help="World_Data.dat (default: DB/World_Data.dat)")
    ap.add_argument("--neutron-lookup", default=None, help="neutron_lookup.dat (default: DB/neutron_lookup.dat)")
    ap.add_argument("--outdir", default="plots", help="Output directory (default: plots)")
    ap.add_argument("--watermark", default="Exploratory", help='Watermark text (default: "Exploratory")')
    # choose somewhat smaller default sizes than before so that the
    # legend and watermark aren't overwhelmingly large on a standard
    # figure; the user can still override these values explicitly.
    ap.add_argument("--legend-fontsize", type=float, default=20,
                    help="Font size for the legend text (default: 20)")
    ap.add_argument("--watermark-font", type=float, default=80,
                    help="Font size for the watermark text (default: 80)")
    ap.add_argument("--world-marker-size", type=float, default=6,
                    help="Marker size for world data points (default: 10)")
    ap.add_argument("--this-marker-size", type=float, default=6,
                    help="Marker size for the ``This Work`` points (default: 11)")
    args = ap.parse_args()

    dbdir = args.dbdir
    world_dat = args.world_dat or os.path.join(dbdir, "World_Data.dat")
    neutron_lookup = args.neutron_lookup or os.path.join(dbdir, "neutron_lookup.dat")

    os.makedirs(args.outdir, exist_ok=True)
    # stash requested legend size on the function for use when
    # ``set_style`` is invoked; this avoids changing the signature
    # of ``set_style`` throughout the module.
    set_style.legend_fontsize = args.legend_fontsize
    set_style()

    # generate the two panels
    plot_ratio(dbdir, world_dat, neutron_lookup, args.outdir,
               args.watermark, watermark_font=args.watermark_font,
               world_marker_size=args.world_marker_size,
               this_marker_size=args.this_marker_size)
    plot_ge(world_dat, neutron_lookup, args.outdir,
            args.watermark, watermark_font=args.watermark_font,
            world_marker_size=args.world_marker_size,
            this_marker_size=args.this_marker_size)

    print("[OK] Wrote:")
    print(" ", os.path.join(args.outdir, "kin234Ratio.png"))
    print(" ", os.path.join(args.outdir, "kin234Ratio.pdf"))
    print(" ", os.path.join(args.outdir, "kin234GEn.png"))
    print(" ", os.path.join(args.outdir, "kin234GEn.pdf"))


if __name__ == "__main__":
    main()