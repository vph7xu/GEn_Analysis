#!/usr/bin/env python3
"""
neutron_overlay_ratio_and_GEn_v7.py

Two overlays (ratio + GEn):
  - previous measurements
  - this work
  - Ye global-fit central curve + 1 sigma band
  - model curves:
      RCQM
      DSE
      GPD-like smooth comparison curve
      VMD = Lomon GKex(02S) exact analytic implementation
      pQCD Lambda=150 MeV
      pQCD Lambda=300 MeV

Run:
  python3 neutron_overlay_ratio_and_GEn_v7.py \
    --dbdir DB \
    --measurements DB/GEn_measurements_grouped_by_method.csv \
    --lookup DB/neutron_lookup.dat \
    --outdir out_plots \
    --watermark Exploratory
"""

from __future__ import annotations

import argparse
import os
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
from matplotlib.lines import Line2D

# ----------------------- Physics constants -----------------------------------
MU_N = -1.9130427
M_N  = 0.9395654133  # GeV

# ----------------------- "This Work" defaults --------------------------------
Q2_THIS = np.array([2.99, 6.81, 9.81], dtype=float)

R_THIS      = np.array([0.4865, 0.7587, 0.6112], dtype=float)
R_THIS_STAT = np.array([0.0404, 0.1247, 0.2590], dtype=float)
R_THIS_SYS  = np.array([0.0361, 0.0321, 0.0670], dtype=float)

GE_THIS      = np.array([0.0179, 0.0057, 0.0020], dtype=float)
GE_THIS_STAT = np.array([0.0015, 0.0009, 0.0009], dtype=float)
GE_THIS_SYS  = np.array([0.0013, 0.0003, 0.0003], dtype=float)


def set_style():
    plt.rcParams.update({
        "figure.dpi": 140,
        "savefig.dpi": 300,
        "font.size": 18,
        "axes.labelsize": 28,
        "axes.titlesize": 34,
        "legend.fontsize": 14,
        "xtick.labelsize": 22,
        "ytick.labelsize": 22,
        "axes.grid": False,
        "axes.spines.top": True,
        "axes.spines.right": True,
    })


def add_watermark(ax, text: str):
    ax.text(
        0.5, 0.50, text,
        transform=ax.transAxes,
        ha="center", va="center",
        fontsize=90, color="red", alpha=0.25,
        fontweight="bold", zorder=0,
    )


def draw_sys_boxes(ax, x, y, y_sys, dx=0.28, alpha=0.20, color="magenta", zorder=1):
    for xi, yi, si in zip(np.asarray(x), np.asarray(y), np.asarray(y_sys)):
        if not np.isfinite(si) or si <= 0:
            continue
        rect = Rectangle(
            (xi - dx, yi - si), 2 * dx, 2 * si,
            facecolor=color, edgecolor=color,
            alpha=alpha, linewidth=1.0, zorder=zorder,
        )
        ax.add_patch(rect)


def load_measurements_csv(path: str) -> pd.DataFrame:
    """
    Expected columns:
      Paper, Q2, ratio, stat, sys, GEn, stat.1, sys.1

    Returns normalized columns:
      paper, Q2, ratio, ratio_stat, ratio_sys, GEn, GEn_stat, GEn_sys
    """
    df = pd.read_csv(path)

    def must(col: str) -> str:
        if col not in df.columns:
            raise ValueError(f"Measurements CSV missing '{col}'. Got: {list(df.columns)}")
        return col

    out = pd.DataFrame({
        "paper": df[must("Paper")].astype(str),
        "Q2": pd.to_numeric(df[must("Q2")], errors="coerce"),
        "ratio": pd.to_numeric(df[must("ratio")], errors="coerce"),
        "ratio_stat": pd.to_numeric(df[must("stat")], errors="coerce"),
        "ratio_sys": pd.to_numeric(df[must("sys")], errors="coerce"),
        "GEn": pd.to_numeric(df[must("GEn")], errors="coerce"),
        "GEn_stat": pd.to_numeric(df[must("stat.1")], errors="coerce"),
        "GEn_sys": pd.to_numeric(df[must("sys.1")], errors="coerce"),
    })

    return out[np.isfinite(out["Q2"])].copy()


def load_ye_bands_from_neutron_lookup(path: str):
    """
    Read neutron_lookup.dat and build Ye bands for:
      - ratio = mu_n GEn / GMn
      - GEn
      - GMn

    Uses *_Par columns as 1 sigma band widths.
    """
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        header_line = ""
        for ln in f:
            if ln.strip():
                header_line = ln.strip()
                break
    hdr = header_line.lstrip("#").strip().split()

    df = pd.read_csv(path, sep=r"\s+", comment="#", header=None, names=hdr)

    q2 = pd.to_numeric(df["Q2"], errors="coerce").to_numpy(float)
    gen_over_gd = pd.to_numeric(df["GEn/GD"], errors="coerce").to_numpy(float)
    dgen_par_over_gd = pd.to_numeric(df["dGEn_Par/GD"], errors="coerce").to_numpy(float)

    gmn_over_mu_gd = pd.to_numeric(df["GMn/muGD"], errors="coerce").to_numpy(float)
    dgmn_par_over_mu_gd = pd.to_numeric(df["dGMn_Par/muGD"], errors="coerce").to_numpy(float)

    gd = (1.0 + q2 / 0.71) ** (-2)

    GE = gen_over_gd * gd
    sGE = dgen_par_over_gd * gd

    GM = gmn_over_mu_gd * (MU_N * gd)
    sGM = dgmn_par_over_mu_gd * (abs(MU_N) * gd)

    ratio = np.full_like(q2, np.nan, dtype=float)
    sratio = np.full_like(q2, np.nan, dtype=float)

    m = (
        np.isfinite(q2) & np.isfinite(GE) & np.isfinite(sGE) &
        np.isfinite(GM) & np.isfinite(sGM) & (GM != 0)
    )

    ratio[m] = MU_N * GE[m] / GM[m]
    term1 = (MU_N / GM[m]) * sGE[m]
    term2 = (MU_N * GE[m] / (GM[m] * GM[m])) * sGM[m]
    sratio[m] = np.sqrt(term1**2 + term2**2)

    order = np.argsort(q2[m])
    q2s = q2[m][order]

    ratio_band = (q2s, ratio[m][order], sratio[m][order])
    ge_band = (q2s, GE[m][order], sGE[m][order])
    gm_band = (q2s, GM[m][order], sGM[m][order])
    return ratio_band, ge_band, gm_band


def _read_two_col_csv(path: str) -> tuple[np.ndarray, np.ndarray]:
    df = pd.read_csv(path, comment="#", header=None)
    if df.shape[1] < 2:
        raise ValueError(f"{path}: expected at least 2 columns, got {df.shape[1]}")
    x = pd.to_numeric(df.iloc[:, 0], errors="coerce").to_numpy(float)
    y = pd.to_numeric(df.iloc[:, 1], errors="coerce").to_numpy(float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    order = np.argsort(x)
    return x[order], y[order]


# ----------------------- RCQM model ------------------------------------------
KAPPA_D = -2.03
KAPPA_U = 1.67

def rcqm_curve(dbdir: str, is_ratio: bool) -> tuple[np.ndarray, np.ndarray]:
    """
    RCQM-style curve using the same flavor-form-factor CSV inputs
    used in your earlier script.
    """
    f1u_x, f1u_y = _read_two_col_csv(os.path.join(dbdir, "F1u.csv"))
    f1d_x, f1d_y = _read_two_col_csv(os.path.join(dbdir, "F1d.csv"))
    f2u_x, f2u_y = _read_two_col_csv(os.path.join(dbdir, "F2u.csv"))
    f2d_x, f2d_y = _read_two_col_csv(os.path.join(dbdir, "F2d.csv"))

    q2 = np.linspace(0.02, 12.0, 500)

    def interp(xgrid, ygrid, xq):
        return np.interp(xq, xgrid, ygrid, left=ygrid[0], right=ygrid[-1])

    out = np.zeros_like(q2)

    for i, Q2 in enumerate(q2):
        F1u = interp(f1u_x, f1u_y, Q2) / (Q2**4)
        F1d = interp(f1d_x, f1d_y, Q2) / (Q2**4)
        F2u = interp(f2u_x, f2u_y, Q2) * KAPPA_U / (Q2**4)
        F2d = interp(f2d_x, f2d_y, Q2) * KAPPA_D / (Q2**4)

        F1 = (2.0 / 3.0) * F1d - (1.0 / 3.0) * F1u
        F2 = (2.0 / 3.0) * F2d - (1.0 / 3.0) * F2u

        tau = Q2 / (4.0 * M_N * M_N)
        GE = F1 - tau * F2
        GM = F1 + F2

        out[i] = (MU_N * GE / GM) if is_ratio else GE

    return q2, out


# ----------------------- DSE-like curve --------------------------------------
def dse_ratio_curve(Q2: np.ndarray) -> np.ndarray:
    """
    DSE-style closed-form comparison curve carried over from your earlier script.
    """
    x = np.asarray(Q2, dtype=float)
    num = (0.4133 * x - 0.1038 * x**2 + 0.01252 * x**3 + 0.007748 * x**4)
    den = (1 + 0.9771 * x - 0.3506 * x**2 + 0.09863 * x**3 + 0.009927 * x**4 + 0.000002173 * x**5)
    y = np.full_like(x, np.nan, dtype=float)
    m = den != 0
    y[m] = num[m] / den[m]
    return y

from math import gamma
from functools import lru_cache

# ----------------------- GPD model: Diehl et al. -----------------------------
# Hybrid paper-based implementation:
#
# H_v^q(x,t):
#   H_v^q(x,t) = q_v(x) * exp[t f_q(x)]
#
# with default-fit profile:
#   f_q(x) = alpha_p * (1-x)^3 * log(1/x) + B_q (1-x)^3 + A_q x (1-x)^2
#
# Parameters from the H_v default fit:
#   alpha_p = 0.9 GeV^-2
#   A_u = 1.22, A_d = 2.59
#   B_u = B_d = 0.59
#
# E_v^q(x,t):
#   E_v^q(x,t) = e_v^q(x) * exp[t g_q(x)]
#
# with
#   e_v^q(x) = N_q * kappa_q * x^{-alpha} * (1-x)^{beta_q}
#
# and
#   g_q(x) = alpha_p * (1-x)^3 * log(1/x) + D_q (1-x)^3 + C_q x (1-x)^2
#
# Parameters from the alternative fit:
#   alpha = 0.55
#   beta_u = 3.99
#   beta_d = 5.59
#   C_u = 1.22, C_d = 2.59
#   D_u = 0.38, D_d = -0.75
#
# This gives a much more stable and physically sensible GPD curve than the
# current Guidal-R2 implementation in this script.

E_U =  2.0 / 3.0
E_D = -1.0 / 3.0

KAPPA_U_GPD = 1.67*1.48
KAPPA_D_GPD = -2.03*1.48

DIEHL_H_PARAMS = {
    "alpha_p": 0.9,   # GeV^-2
    "A_u": 1.22,
    "A_d": 2.59,
    "B_u": 0.59,
    "B_d": 0.59,
}

DIEHL_E_PARAMS = {
    "alpha": 0.55,
    "alpha_p": 0.9,   # GeV^-2
    "beta_u": 3.99,
    "beta_d": 5.59,
    "C_u": 1.22,
    "C_d": 2.59,
    "D_u": 0.38,
    "D_d": -0.75,
}

def uv_mrst2002(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    return 0.262 * x**(-0.69) * (1.0 - x)**3.50 * (1.0 + 3.83 * x**0.5 + 37.65 * x)

def dv_mrst2002(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    return 0.061 * x**(-0.65) * (1.0 - x)**4.03 * (1.0 + 49.05 * x**0.5 + 8.65 * x)

def _xgrid_for_gpd(n: int = 4000) -> np.ndarray:
    return np.linspace(1e-5, 1.0 - 1e-5, n)

# ---- H_v part ----
def diehl_f_u(x: np.ndarray, p: dict = DIEHL_H_PARAMS) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    return (
        p["alpha_p"] * (1.0 - x)**3 * np.log(1.0 / x)
        + p["B_u"] * (1.0 - x)**3
        + p["A_u"] * x * (1.0 - x)**2
    )

def diehl_f_d(x: np.ndarray, p: dict = DIEHL_H_PARAMS) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    return (
        p["alpha_p"] * (1.0 - x)**3 * np.log(1.0 / x)
        + p["B_d"] * (1.0 - x)**3
        + p["A_d"] * x * (1.0 - x)**2
    )

def H_u_diehl(x: np.ndarray, Q2: float, p: dict = DIEHL_H_PARAMS) -> np.ndarray:
    return uv_mrst2002(x) * np.exp(-Q2 * diehl_f_u(x, p))

def H_d_diehl(x: np.ndarray, Q2: float, p: dict = DIEHL_H_PARAMS) -> np.ndarray:
    return dv_mrst2002(x) * np.exp(-Q2 * diehl_f_d(x, p))

# ---- E_v part ----
def diehl_Nq(alpha: float, beta_q: float) -> float:
    return gamma(2.0 - alpha + beta_q) / (gamma(1.0 - alpha) * gamma(1.0 + beta_q))

def diehl_ev_forward_u(x: np.ndarray, p: dict = DIEHL_E_PARAMS) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    Nq = diehl_Nq(p["alpha"], p["beta_u"])
    return Nq * KAPPA_U_GPD * x**(-p["alpha"]) * (1.0 - x)**p["beta_u"]

def diehl_ev_forward_d(x: np.ndarray, p: dict = DIEHL_E_PARAMS) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    Nq = diehl_Nq(p["alpha"], p["beta_d"])
    return Nq * KAPPA_D_GPD * x**(-p["alpha"]) * (1.0 - x)**p["beta_d"]

def diehl_g_u(x: np.ndarray, p: dict = DIEHL_E_PARAMS) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    return (
        p["alpha_p"] * (1.0 - x)**3 * np.log(1.0 / x)
        + p["D_u"] * (1.0 - x)**3
        + p["C_u"] * x * (1.0 - x)**2
    )

def diehl_g_d(x: np.ndarray, p: dict = DIEHL_E_PARAMS) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    return (
        p["alpha_p"] * (1.0 - x)**3 * np.log(1.0 / x)
        + p["D_d"] * (1.0 - x)**3
        + p["C_d"] * x * (1.0 - x)**2
    )

def E_u_diehl(x: np.ndarray, Q2: float, p: dict = DIEHL_E_PARAMS) -> np.ndarray:
    return diehl_ev_forward_u(x, p) * np.exp(-Q2 * diehl_g_u(x, p))

def E_d_diehl(x: np.ndarray, Q2: float, p: dict = DIEHL_E_PARAMS) -> np.ndarray:
    return diehl_ev_forward_d(x, p) * np.exp(-Q2 * diehl_g_d(x, p))

# ---- flavor form factors ----
def F1u_gpd(Q2: float) -> float:
    x = _xgrid_for_gpd()
    return float(np.trapz(H_u_diehl(x, Q2), x))

def F1d_gpd(Q2: float) -> float:
    x = _xgrid_for_gpd()
    return float(np.trapz(H_d_diehl(x, Q2), x))

def F2u_gpd(Q2: float) -> float:
    x = _xgrid_for_gpd()
    return float(np.trapz(E_u_diehl(x, Q2), x))

def F2d_gpd(Q2: float) -> float:
    x = _xgrid_for_gpd()
    return float(np.trapz(E_d_diehl(x, Q2), x))

def F1n_gpd(Q2: float) -> float:
    return E_U * F1d_gpd(Q2) + E_D * F1u_gpd(Q2)

def F2n_gpd(Q2: float) -> float:
    return E_U * F2d_gpd(Q2) + E_D * F2u_gpd(Q2)

def gpd_GEn_curve(Q2: np.ndarray) -> np.ndarray:
    q2 = np.asarray(Q2, dtype=float)
    out = np.full_like(q2, np.nan, dtype=float)
    for i, qq in enumerate(q2):
        tau = qq / (4.0 * M_N * M_N)
        out[i] = F1n_gpd(float(qq)) - tau * F2n_gpd(float(qq))
    return out

def gpd_GMn_curve(Q2: np.ndarray) -> np.ndarray:
    q2 = np.asarray(Q2, dtype=float)
    out = np.full_like(q2, np.nan, dtype=float)
    for i, qq in enumerate(q2):
        out[i] = F1n_gpd(float(qq)) + F2n_gpd(float(qq))
    return out

def gpd_ratio_curve_exact(Q2: np.ndarray) -> np.ndarray:
    ge = gpd_GEn_curve(Q2)
    gm = gpd_GMn_curve(Q2)
    out = np.full_like(np.asarray(Q2, dtype=float), np.nan, dtype=float)
    m = np.isfinite(gm) & (np.abs(gm) > 1e-10)
    out[m] = MU_N * ge[m] / gm[m]
    return out


# ----------------------- Exact VMD: Lomon GKex(02S) --------------------------
# Functional form from:
#   E. L. Lomon, Phys. Rev. C 66, 045501 (2002)
#
# How it is built:
#   1) The paper gives F1^iv, F2^iv, F1^is, F2^is in Eq. (6),
#      with meson and direct terms.
#   2) The form factors F1^a, F2^a, F1^D, F2^D and the phi-specific pieces
#      are given in Eq. (7), with Q^2-tilde enforcing the pQCD limit.
#   3) Neutron Dirac/Pauli FFs follow from isospin decomposition:
#        F1n = (F1^is - F1^iv)/2
#        F2n = (F2^is - F2^iv)/2
#   4) Sachs FFs:
#        GEn = F1n - tau F2n
#        GMn = F1n + F2n
#      tau = Q2/(4 M_N^2)
#
# Parameter set below is GKex(02S), Table I of the same paper.

KAPPA_V = 3.706
KAPPA_S = -0.120

M_RHO    = 0.776
M_OMEGA  = 0.784
M_PHI    = 1.019
M_RHOP   = 1.45
M_OMEGAP = 1.419

GKEX02S = {
    "g_rhop_over_f_rhop": 0.0401,
    "kappa_rhop": 6.8190,

    "g_omega_over_f_omega": 0.6739,
    "kappa_omega": 0.8762,

    "g_phi_over_f_phi": -0.1676,
    "kappa_phi": 7.0172,
    "mu_phi": 0.8544,

    "g_omegap_over_f_omegap": 0.2552,
    "kappa_omegap": 1.4916,

    "Lambda1": 0.9407,
    "LambdaD": 1.2111,
    "Lambda2": 2.7891,
    "LambdaQCD": 0.150,

    "N": 1.0,
}

GKEX05 = {
    # 2006 updated fit parameters (same GKex functional form, new fit)
    "g_rhop_over_f_rhop": 0.007208,
    "kappa_rhop": 12.0,

    "g_omega_over_f_omega": 0.7021,
    "kappa_omega": 0.4027,

    "g_phi_over_f_phi": -0.1711,
    "kappa_phi": 0.01,
    "mu_phi": 0.2,

    "g_omegap_over_f_omegap": 0.164,
    "kappa_omegap": -2.973,

    "Lambda1": 0.93088,
    "LambdaD": 1.181,
    "Lambda2": 2.6115,
    "LambdaQCD": 0.150,

    "N": 1.0,
}


def gkex_q2_tilde(Q2: np.ndarray, LambdaD: float, LambdaQCD: float) -> np.ndarray:
    Q2 = np.asarray(Q2, dtype=float)
    return Q2 * np.log((LambdaD**2 + Q2) / (LambdaQCD**2)) / np.log(LambdaD**2 / (LambdaQCD**2))


def gkex_F1_a(Q2: np.ndarray, L1: float, L2: float, LD: float, LQCD: float) -> np.ndarray:
    qt2 = gkex_q2_tilde(Q2, LD, LQCD)
    return (L1**2 / (L1**2 + qt2)) * (L2**2 / (L2**2 + qt2))


def gkex_F2_a(Q2: np.ndarray, L1: float, L2: float, LD: float, LQCD: float) -> np.ndarray:
    qt2 = gkex_q2_tilde(Q2, LD, LQCD)
    return (L1**2 / (L1**2 + qt2)) * (L2**2 / (L2**2 + qt2))**2


def gkex_F1_D(Q2: np.ndarray, LD: float, L2: float, LQCD: float) -> np.ndarray:
    qt2 = gkex_q2_tilde(Q2, LD, LQCD)
    return (LD**2 / (LD**2 + qt2)) * (L2**2 / (L2**2 + qt2))


def gkex_F2_D(Q2: np.ndarray, LD: float, L2: float, LQCD: float) -> np.ndarray:
    qt2 = gkex_q2_tilde(Q2, LD, LQCD)
    return (LD**2 / (LD**2 + qt2)) * (L2**2 / (L2**2 + qt2))**2


def gkex_F1_phi(Q2: np.ndarray, L1: float, L2: float, LD: float, LQCD: float) -> np.ndarray:
    Q2 = np.asarray(Q2, dtype=float)
    base = gkex_F1_a(Q2, L1, L2, LD, LQCD)
    return base * (Q2 / (L1**2 + Q2))**1.5


def gkex_F2_phi(Q2: np.ndarray, L1: float, L2: float, LD: float, LQCD: float, mu_phi: float) -> np.ndarray:
    Q2 = np.asarray(Q2, dtype=float)
    base = gkex_F2_a(Q2, L1, L2, LD, LQCD)
    factor = ((L1**2 / mu_phi**2) * ((Q2 + mu_phi**2) / (L1**2 + Q2)))**1.5
    return base * factor


def gkex_rho_prefactor_F1(Q2: np.ndarray, N: float) -> np.ndarray:
    Q2 = np.asarray(Q2, dtype=float)
    return (N / 2.0) * ((1.0317 + 0.0875 * (1.0 + Q2 / 0.3176)**(-2.0)) / (1.0 + Q2 / 0.5496))


def gkex_rho_prefactor_F2(Q2: np.ndarray, N: float) -> np.ndarray:
    Q2 = np.asarray(Q2, dtype=float)
    return (N / 2.0) * ((5.7824 + 0.3907 * (1.0 + Q2 / 0.1422)**(-1.0)) / (1.0 + Q2 / 0.5362))


def gkex_F1_iv(Q2: np.ndarray, p: dict = GKEX05) -> np.ndarray:
    grhop = p["g_rhop_over_f_rhop"]
    L1, LD, L2, LQCD = p["Lambda1"], p["LambdaD"], p["Lambda2"], p["LambdaQCD"]
    N = p["N"]

    F1r = gkex_F1_a(Q2, L1, L2, LD, LQCD)

    term_rho  = gkex_rho_prefactor_F1(Q2, N) * F1r
    term_rhop = grhop * (M_RHOP**2 / (M_RHOP**2 + Q2)) * F1r
    term_D    = (1.0 - 1.1192 * N / 2.0 - grhop) * gkex_F1_D(Q2, LD, L2, LQCD)

    return term_rho + term_rhop + term_D


def gkex_F2_iv(Q2: np.ndarray, p: dict = GKEX05) -> np.ndarray:
    grhop = p["g_rhop_over_f_rhop"]
    krhop = p["kappa_rhop"]
    L1, LD, L2, LQCD = p["Lambda1"], p["LambdaD"], p["Lambda2"], p["LambdaQCD"]
    N = p["N"]

    F2r = gkex_F2_a(Q2, L1, L2, LD, LQCD)

    term_rho  = gkex_rho_prefactor_F2(Q2, N) * F2r
    term_rhop = krhop * grhop * (M_RHOP**2 / (M_RHOP**2 + Q2)) * F2r
    term_D    = (KAPPA_V - 6.1731 * N / 2.0 - krhop * grhop) * gkex_F2_D(Q2, LD, L2, LQCD)

    return term_rho + term_rhop + term_D


def gkex_F1_is(Q2: np.ndarray, p: dict = GKEX05) -> np.ndarray:
    gow   = p["g_omega_over_f_omega"]
    gop   = p["g_omegap_over_f_omegap"]
    gphi  = p["g_phi_over_f_phi"]
    L1, LD, L2, LQCD = p["Lambda1"], p["LambdaD"], p["Lambda2"], p["LambdaQCD"]

    F1om = gkex_F1_a(Q2, L1, L2, LD, LQCD)
    F1ph = gkex_F1_phi(Q2, L1, L2, LD, LQCD)

    term_omega  = gow  * (M_OMEGA**2  / (M_OMEGA**2  + Q2)) * F1om
    term_omegap = gop  * (M_OMEGAP**2 / (M_OMEGAP**2 + Q2)) * F1om
    term_phi    = gphi * (M_PHI**2    / (M_PHI**2    + Q2)) * F1ph
    term_D      = (1.0 - gow - gop) * gkex_F1_D(Q2, LD, L2, LQCD)

    return term_omega + term_omegap + term_phi + term_D


def gkex_F2_is(Q2: np.ndarray, p: dict = GKEX05) -> np.ndarray:
    gow   = p["g_omega_over_f_omega"]
    kow   = p["kappa_omega"]
    gop   = p["g_omegap_over_f_omegap"]
    kop   = p["kappa_omegap"]
    gphi  = p["g_phi_over_f_phi"]
    kphi  = p["kappa_phi"]
    mu_phi = p["mu_phi"]
    L1, LD, L2, LQCD = p["Lambda1"], p["LambdaD"], p["Lambda2"], p["LambdaQCD"]

    F2om = gkex_F2_a(Q2, L1, L2, LD, LQCD)
    F2ph = gkex_F2_phi(Q2, L1, L2, LD, LQCD, mu_phi)

    term_omega  = kow  * gow  * (M_OMEGA**2  / (M_OMEGA**2  + Q2)) * F2om
    term_omegap = kop  * gop  * (M_OMEGAP**2 / (M_OMEGAP**2 + Q2)) * F2om
    term_phi    = kphi * gphi * (M_PHI**2    / (M_PHI**2    + Q2)) * F2ph
    term_D      = (KAPPA_S - kow * gow - kop * gop - kphi * gphi) * gkex_F2_D(Q2, LD, L2, LQCD)

    return term_omega + term_omegap + term_phi + term_D


def gkex_F1n(Q2: np.ndarray, p: dict = GKEX05) -> np.ndarray:
    return 0.5 * (gkex_F1_is(Q2, p) - gkex_F1_iv(Q2, p))


def gkex_F2n(Q2: np.ndarray, p: dict = GKEX05) -> np.ndarray:
    return 0.5 * (gkex_F2_is(Q2, p) - gkex_F2_iv(Q2, p))


def lomon_gkex_GEn(Q2: np.ndarray, p: dict = GKEX05) -> np.ndarray:
    Q2 = np.asarray(Q2, dtype=float)
    tau = Q2 / (4.0 * M_N * M_N)
    return gkex_F1n(Q2, p) - tau * gkex_F2n(Q2, p)


def lomon_gkex_GMn(Q2: np.ndarray, p: dict = GKEX05) -> np.ndarray:
    return gkex_F1n(Q2, p) + gkex_F2n(Q2, p)


def lomon_gkex_ratio_curve(Q2: np.ndarray, p: dict = GKEX05) -> np.ndarray:
    ge = lomon_gkex_GEn(Q2, p)
    gm = lomon_gkex_GMn(Q2, p)

    out = np.full_like(np.asarray(Q2, dtype=float), np.nan, dtype=float)
    m = np.isfinite(gm) & (gm != 0.0)
    out[m] = MU_N * ge[m] / gm[m]
    return out


# ----------------------- Exact pQCD plotting prescription --------------------
def pqcd_ratio_curve_exact(Q2: np.ndarray, lam_gev: float, q2_ref: float = 1.5, y_ref: float = 0.30) -> np.ndarray:
    """
    pQCD-inspired curve used in Riordan-style comparison plots.

    Shape:
      log^2(Q2/Lambda^2) / Q2

    Then normalized so the plotted ratio passes through y_ref at q2_ref.
    """
    x = np.asarray(Q2, dtype=float)

    def shape(z):
        z = np.asarray(z, dtype=float)
        out = np.full_like(z, np.nan, dtype=float)
        m = z > 0.0
        out[m] = np.log(z[m] / (lam_gev * lam_gev))**2 / z[m]
        return out

    y = shape(x)
    y0 = float(shape(np.array([q2_ref]))[0])

    if not np.isfinite(y0) or y0 == 0.0:
        return np.full_like(x, np.nan, dtype=float)

    return y_ref * y / y0


def ratio_to_gen(ratio: np.ndarray, q2: np.ndarray, gm_q2: np.ndarray, gm_vals: np.ndarray) -> np.ndarray:
    """
    Convert ratio = mu_n GEn / GMn into GEn using an input GMn(Q2) curve.
    """
    interp_gm = np.interp(q2, gm_q2, gm_vals, left=np.nan, right=np.nan)
    out = np.full_like(q2, np.nan, dtype=float)
    m = np.isfinite(ratio) & np.isfinite(interp_gm)
    out[m] = ratio[m] * interp_gm[m] / MU_N
    return out


def marker_map_from_papers(papers: list[str]) -> dict[str, str]:
    markers = ["o", "s", "^", "v", "D", "P", "X", "<", ">", "h", "H", "p", "*"]
    return {p: markers[i % len(markers)] for i, p in enumerate(sorted(set(papers)))}


def color_map_from_papers(papers: list[str]) -> dict[str, str]:
    colors = ["red", "black", "blue"]
    return {p: colors[i % len(colors)] for i, p in enumerate(sorted(set(papers)))}


def write_plotted_values_csv(
    out_csv: str,
    df: pd.DataFrame,
    ratio_band: tuple[np.ndarray, np.ndarray, np.ndarray],
    ge_band: tuple[np.ndarray, np.ndarray, np.ndarray],
    this_q2: np.ndarray,
    this_ratio: np.ndarray,
    this_ratio_stat: np.ndarray,
    this_ratio_sys: np.ndarray,
    this_ge: np.ndarray,
    this_ge_stat: np.ndarray,
    this_ge_sys: np.ndarray,
):
    rows = []

    def combine_errors(stat, sys):
        stat_finite = np.isfinite(stat)
        sys_finite = np.isfinite(sys)
        if stat_finite and sys_finite:
            return float(np.sqrt(stat**2 + sys**2))
        if stat_finite:
            return float(stat)
        if sys_finite:
            return float(sys)
        return np.nan

    for _, r in df.iterrows():
        q2 = pd.to_numeric(r.get("Q2"), errors="coerce")
        paper = str(r.get("paper", ""))

        ratio = pd.to_numeric(r.get("ratio"), errors="coerce")
        ratio_stat = pd.to_numeric(r.get("ratio_stat"), errors="coerce")
        ratio_sys = pd.to_numeric(r.get("ratio_sys"), errors="coerce")

        if np.isfinite(q2) and np.isfinite(ratio):
            rows.append({
                "dataset": "previous_measurement",
                "panel": "ratio",
                "paper": paper,
                "Q2": float(q2),
                "value": float(ratio),
                "stat_err": float(ratio_stat) if np.isfinite(ratio_stat) else np.nan,
                "sys_err": float(ratio_sys) if np.isfinite(ratio_sys) else np.nan,
                "total_err": combine_errors(ratio_stat, ratio_sys),
            })

        ge = pd.to_numeric(r.get("GEn"), errors="coerce")
        ge_stat = pd.to_numeric(r.get("GEn_stat"), errors="coerce")
        ge_sys = pd.to_numeric(r.get("GEn_sys"), errors="coerce")

        if np.isfinite(q2) and np.isfinite(ge):
            rows.append({
                "dataset": "previous_measurement",
                "panel": "GEn",
                "paper": paper,
                "Q2": float(q2),
                "value": float(ge),
                "stat_err": float(ge_stat) if np.isfinite(ge_stat) else np.nan,
                "sys_err": float(ge_sys) if np.isfinite(ge_sys) else np.nan,
                "total_err": combine_errors(ge_stat, ge_sys),
            })

    for q2, y, est, esys in zip(this_q2, this_ratio, this_ratio_stat, this_ratio_sys):
        rows.append({
            "dataset": "this_work",
            "panel": "ratio",
            "paper": "This Work",
            "Q2": float(q2),
            "value": float(y),
            "stat_err": float(est) if np.isfinite(est) else np.nan,
            "sys_err": float(esys) if np.isfinite(esys) else np.nan,
            "total_err": combine_errors(est, esys),
        })

    for q2, y, est, esys in zip(this_q2, this_ge, this_ge_stat, this_ge_sys):
        rows.append({
            "dataset": "this_work",
            "panel": "GEn",
            "paper": "This Work",
            "Q2": float(q2),
            "value": float(y),
            "stat_err": float(est) if np.isfinite(est) else np.nan,
            "sys_err": float(esys) if np.isfinite(esys) else np.nan,
            "total_err": combine_errors(est, esys),
        })

    rq2, ry, rs = ratio_band
    for q2, y, s in zip(rq2, ry, rs):
        rows.append({
            "dataset": "global_fit_ye",
            "panel": "ratio",
            "paper": "Ye",
            "Q2": float(q2),
            "value": float(y),
            "stat_err": np.nan,
            "sys_err": float(s) if np.isfinite(s) else np.nan,
            "total_err": float(s) if np.isfinite(s) else np.nan,
        })

    gq2, gy, gs = ge_band
    for q2, y, s in zip(gq2, gy, gs):
        rows.append({
            "dataset": "global_fit_ye",
            "panel": "GEn",
            "paper": "Ye",
            "Q2": float(q2),
            "value": float(y),
            "stat_err": np.nan,
            "sys_err": float(s) if np.isfinite(s) else np.nan,
            "total_err": float(s) if np.isfinite(s) else np.nan,
        })

    out_df = pd.DataFrame(rows, columns=[
        "dataset", "panel", "paper", "Q2", "value",
        "stat_err", "sys_err", "total_err"
    ])
    out_df.to_csv(out_csv, index=False)
    print(f"[OK] Wrote plotted values CSV: {out_csv}")


def plot_overlay(
    *,
    df: pd.DataFrame,
    band_q2: np.ndarray,
    band_y: np.ndarray,
    band_s: np.ndarray,
    out_png: str,
    out_pdf: str,
    title: str,
    ylabel: str,
    y_col: str,
    xlim=(0, 12),
    ylim=None,
    watermark: str | None,
    show_models: bool,
    dbdir: str,
    gm_band: tuple[np.ndarray, np.ndarray] | None = None,
    legend_loc: str = "upper left",
    show_extra_models: bool = True,
):
    # Fill missing ratio values from measured GEn and global-fit GMn,
    # propagating both stat and sys separately.
    if y_col == "ratio" and gm_band is not None:
        gm_q2, gm_vals = gm_band
        mask = (~np.isfinite(df["ratio"])) & np.isfinite(df["GEn"]) & np.isfinite(df["Q2"])
        if mask.any():
            q2_fill = df.loc[mask, "Q2"].to_numpy(float)
            interp_gm = np.interp(q2_fill, gm_q2, gm_vals, left=np.nan, right=np.nan)

            valid = np.isfinite(interp_gm) & (interp_gm != 0.0)
            nfilled = np.count_nonzero(valid)

            if nfilled > 0:
                ge_vals = df.loc[mask, "GEn"].to_numpy(float)

                ratio_vals = np.full_like(ge_vals, np.nan, dtype=float)
                ratio_vals[valid] = MU_N * ge_vals[valid] / interp_gm[valid]
                df.loc[mask, "ratio"] = ratio_vals

                ge_stat = df.loc[mask, "GEn_stat"].to_numpy(float)
                ratio_stat = np.full_like(ge_vals, np.nan, dtype=float)
                good_stat = valid & np.isfinite(ge_stat)
                ratio_stat[good_stat] = np.abs(MU_N / interp_gm[good_stat]) * ge_stat[good_stat]
                df.loc[mask, "ratio_stat"] = ratio_stat

                ge_sys = df.loc[mask, "GEn_sys"].to_numpy(float)
                ratio_sys = np.full_like(ge_vals, np.nan, dtype=float)
                good_sys = valid & np.isfinite(ge_sys)
                ratio_sys[good_sys] = np.abs(MU_N / interp_gm[good_sys]) * ge_sys[good_sys]
                df.loc[mask, "ratio_sys"] = ratio_sys

                warnings.warn(f"Filled {nfilled} missing ratio values from GE+GM global fit.")

    fig, ax = plt.subplots(figsize=(11.0, 8.0))

    ax.set_title(title, pad=18)
    ax.set_xlabel(r"$Q^2\;(\mathrm{GeV}^2)$")
    ax.set_ylabel(ylabel)
    ax.set_xlim(*xlim)
    if ylim is not None:
        ax.set_ylim(*ylim)

    ye_band = ax.fill_between(band_q2, band_y - band_s, band_y + band_s, alpha=0.25, linewidth=0)
    ye_line, = ax.plot(band_q2, band_y, linewidth=2.0)

    rcqm_line = None
    dse_line = None
    gpd_line = None
    vmd_line = None
    pqcd150_line = None
    pqcd300_line = None

    if show_models:
        q2_model = np.linspace(0.02, xlim[1], 500)

        # RCQM
        try:
            q2m, ym = rcqm_curve(dbdir, is_ratio=(y_col == "ratio"))
            rcqm_line, = ax.plot(q2m, ym, linestyle=":", linewidth=2.0, color="limegreen")
        except Exception as e:
            warnings.warn(f"RCQM model unavailable ({e}); skipping RCQM curve.")

        # DSE
        if y_col == "ratio":
            y_dse = dse_ratio_curve(q2_model)
        else:
            if gm_band is not None:
                gm_q2, gm_vals = gm_band
                y_dse = ratio_to_gen(dse_ratio_curve(q2_model), q2_model, gm_q2, gm_vals)
            else:
                y_dse = None
        if y_dse is not None:
            dse_line, = ax.plot(q2_model, y_dse, linestyle=(0, (10, 4)), linewidth=1.4, color="red")

        # Exact VMD = Lomon GKex(02S)
        if y_col == "ratio":
            y_vmd = lomon_gkex_ratio_curve(q2_model)
        else:
            y_vmd = lomon_gkex_GEn(q2_model)
        vmd_line, = ax.plot(q2_model, y_vmd, linestyle="--", linewidth=1.8, color="blue")

        # GPD-like comparison
        # Exact GPD model from Guidal et al. R2
        if y_col == "ratio":
            y_gpd = gpd_ratio_curve_exact(q2_model)
        else:
            y_gpd = gpd_GEn_curve(q2_model)

        gpd_line, = ax.plot(q2_model, y_gpd, linestyle="-", linewidth=1.4, color="red")

        if show_extra_models:

            # pQCD Lambda = 150 MeV
            y_pqcd150_ratio = pqcd_ratio_curve_exact(q2_model, lam_gev=0.15, q2_ref=1.5, y_ref=0.30)
            # pQCD Lambda = 300 MeV
            y_pqcd300_ratio = pqcd_ratio_curve_exact(q2_model, lam_gev=0.30, q2_ref=1.5, y_ref=0.30)

            if y_col == "ratio":
                y_pqcd150 = y_pqcd150_ratio
                y_pqcd300 = y_pqcd300_ratio
            else:
                if gm_band is not None:
                    gm_q2, gm_vals = gm_band
                    y_pqcd150 = ratio_to_gen(y_pqcd150_ratio, q2_model, gm_q2, gm_vals)
                    y_pqcd300 = ratio_to_gen(y_pqcd300_ratio, q2_model, gm_q2, gm_vals)
                else:
                    y_pqcd150 = None
                    y_pqcd300 = None

            if y_pqcd150 is not None:
                pqcd150_line, = ax.plot(q2_model, y_pqcd150, linestyle="-.", linewidth=1.4, color="magenta")
            if y_pqcd300 is not None:
                pqcd300_line, = ax.plot(q2_model, y_pqcd300, linestyle=(0, (14, 8)), linewidth=1.4, color="black")

    papers = df["paper"].dropna().astype(str).tolist()
    p2m = marker_map_from_papers(papers)
    p2c = color_map_from_papers(papers)

    for paper, d in df.groupby("paper", sort=True):
        d = d[np.isfinite(d["Q2"]) & np.isfinite(d[y_col])]
        if d.empty:
            continue

        if y_col == "ratio":
            ystat = d["ratio_stat"].to_numpy(float)
            ysys = d["ratio_sys"].to_numpy(float)
        else:
            ystat = d["GEn_stat"].to_numpy(float)
            ysys = d["GEn_sys"].to_numpy(float)

        ystat = np.where(np.isfinite(ystat), ystat, 0.0)
        ysys = np.where(np.isfinite(ysys), ysys, 0.0)
        yerr = np.sqrt(ystat**2 + ysys**2)

        col = p2c.get(str(paper), "black")
        ax.errorbar(
            d["Q2"].to_numpy(float),
            d[y_col].to_numpy(float),
            yerr=yerr,
            fmt=p2m.get(str(paper), "o"),
            markersize=5,
            markerfacecolor=col,
            markeredgecolor=col,
            ecolor=col,
            elinewidth=1.4,
            capsize=0,
            linestyle="none",
            zorder=2,
        )

    if y_col == "ratio":
        draw_sys_boxes(ax, Q2_THIS, R_THIS, R_THIS_SYS, dx=0.28, alpha=0.15, color="magenta", zorder=1)
        total_err = np.sqrt(R_THIS_STAT**2 + R_THIS_SYS**2)
        ax.errorbar(
            Q2_THIS, R_THIS, yerr=total_err,
            fmt="o", markersize=6,
            markerfacecolor="magenta", markeredgecolor="magenta",
            ecolor="magenta", elinewidth=2.0, capsize=0,
            linestyle="none",
            zorder=3,
        )
    else:
        draw_sys_boxes(ax, Q2_THIS, GE_THIS, GE_THIS_SYS, dx=0.28, alpha=0.15, color="magenta", zorder=1)
        total_err = np.sqrt(GE_THIS_STAT**2 + GE_THIS_SYS**2)
        ax.errorbar(
            Q2_THIS, GE_THIS, yerr=total_err,
            fmt="o", markersize=6,
            markerfacecolor="magenta", markeredgecolor="magenta",
            ecolor="magenta", elinewidth=2.0, capsize=0,
            linestyle="none",
            zorder=3,
        )

    if watermark:
        add_watermark(ax, watermark)

    ax.minorticks_on()
    ax.tick_params(axis="both", which="major", direction="in", length=9)
    ax.tick_params(axis="both", which="minor", direction="in", length=5)

    handles = [
        Patch(facecolor=ye_band.get_facecolor()[0], edgecolor="none", alpha=0.25),
        Line2D([0], [0], color=ye_line.get_color(), lw=2.0),
    ]
    labels = [
        "Global Fit - Ye (1$\\sigma$)",
        "Global Fit - Ye (central)",
    ]

    if show_models and rcqm_line is not None:
        handles.append(Line2D([0], [0], color="limegreen", lw=2.0, ls=":"))
        labels.append("RCQM")
    if show_models and gpd_line is not None:
        handles.append(Line2D([0], [0], color="red", lw=1.4, ls="-"))
        labels.append("GPD (Diehl)")
    if show_models and vmd_line is not None:
        handles.append(Line2D([0], [0], color="blue", lw=1.8, ls="--"))
        labels.append("VMD (Lomon GKex05)")
    if show_models and dse_line is not None:
        handles.append(Line2D([0], [0], color="red", lw=1.4, ls=(0, (10, 4))))
        labels.append("DSE")
    if show_models and pqcd150_line is not None:
        handles.append(Line2D([0], [0], color="magenta", lw=1.4, ls="-."))
        labels.append(r"pQCD, $\Lambda=150$ MeV")
    if show_models and pqcd300_line is not None:
        handles.append(Line2D([0], [0], color="black", lw=1.4, ls=(0, (14, 8))))
        labels.append(r"pQCD, $\Lambda=300$ MeV")

    for paper in sorted(set(papers)):
        col = p2c.get(paper, "black")
        handles.append(Line2D(
            [0], [0],
            marker=p2m.get(paper, "o"),
            linestyle="none",
            markersize=6,
            markerfacecolor=col,
            markeredgecolor=col,
            color=col,
        ))
        labels.append(paper)

    handles.append(Line2D([0], [0], marker="o", linestyle="none", markersize=6,
                          markerfacecolor="magenta", markeredgecolor="magenta", color="magenta"))
    labels.append("This Work")

    handles.append(Patch(facecolor="magenta", edgecolor="magenta", alpha=0.15))
    labels.append("This Work (syst.)")

    ax.legend(handles, labels, frameon=False, loc=legend_loc, ncol=2, fontsize=14)

    fig.tight_layout()
    fig.savefig(out_png)
    fig.savefig(out_pdf)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dbdir", default="DB/", help="DB directory")
    ap.add_argument("--measurements", default=None, help="Measurements CSV")
    ap.add_argument("--lookup", default=None, help="neutron_lookup.dat")
    ap.add_argument("--outdir", default="out_plots", help="Output directory")
    ap.add_argument("--watermark", default="Exploratory")
    args = ap.parse_args()

    set_style()

    meas_path = args.measurements or os.path.join(args.dbdir, "GEn_measurements_grouped_by_method.csv")
    lookup_path = args.lookup or os.path.join(args.dbdir, "neutron_lookup.dat")

    os.makedirs(args.outdir, exist_ok=True)

    df = load_measurements_csv(meas_path)
    (rq2, ry, rs), (gq2, gy, gs), (gmq2, gm, sgm) = load_ye_bands_from_neutron_lookup(lookup_path)

    plot_overlay(
        df=df,
        band_q2=rq2,
        band_y=ry,
        band_s=rs,
        out_png=os.path.join(args.outdir, "ratio_overlay_2026FEBstyle.png"),
        out_pdf=os.path.join(args.outdir, "ratio_overlay_2026FEBstyle.pdf"),
        title=r"$G_E/G_M$ Neutron Results",
        ylabel=r"$\mu_n\,G_E^n/G_M^n$",
        y_col="ratio",
        xlim=(0, 12),
        ylim=(0, 1.8),
        watermark=args.watermark,
        show_models=True,
        dbdir=args.dbdir,
        gm_band=(gmq2, gm),
        legend_loc="upper left",
        show_extra_models=False,
    )

    plot_overlay(
        df=df,
        band_q2=gq2,
        band_y=gy,
        band_s=gs,
        out_png=os.path.join(args.outdir, "GEn_overlay_2026FEBstyle.png"),
        out_pdf=os.path.join(args.outdir, "GEn_overlay_2026FEBstyle.pdf"),
        title=r"$G_E^n$ Neutron Results",
        ylabel=r"$G_E^n$",
        y_col="GEn",
        xlim=(0, 12),
        ylim=(0, 0.08),
        watermark=args.watermark,
        show_models=True,
        dbdir=args.dbdir,
        gm_band=(gmq2, gm),
        legend_loc="upper right",
        show_extra_models=False,
    )

    write_plotted_values_csv(
        out_csv=os.path.join(args.outdir, "plotted_values_and_errors.csv"),
        df=df,
        ratio_band=(rq2, ry, rs),
        ge_band=(gq2, gy, gs),
        this_q2=Q2_THIS,
        this_ratio=R_THIS,
        this_ratio_stat=R_THIS_STAT,
        this_ratio_sys=R_THIS_SYS,
        this_ge=GE_THIS,
        this_ge_stat=GE_THIS_STAT,
        this_ge_sys=GE_THIS_SYS,
    )

    print("[OK] Wrote:")
    print("  ", os.path.join(args.outdir, "ratio_overlay_2026FEBstyle.pdf"))
    print("  ", os.path.join(args.outdir, "GEn_overlay_2026FEBstyle.pdf"))
    print("  ", os.path.join(args.outdir, "plotted_values_and_errors.csv"))


if __name__ == "__main__":
    main()