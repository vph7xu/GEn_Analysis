#!/usr/bin/env python3
"""
optimize_cuts_dxdy_ellipse_opt.py

Cut optimization for GEn-style quasi-elastic selection with:
  - Chunked I/O with uproot.iterate (avoids OOM)
  - Same C++-like preselection logic you used before (kept exactly as in v2)
  - coin_time cut applied to DATA only (SIM/BKG assumed in-time)
  - dx template fit with NNLS on unit-area PDFs (proton, neutron, optional bkg)
  - Dilution fn computed from the fit inside a dx window
  - Objective: minimize sigma_phys ~= sigma_Araw / fn
      where sigma_Araw is pure counting-stat error from DATA helicity counts

NEW vs your current script:
  ✅ dx/dy cut is now an ellipse centered at (dx0,dy0) (default 0,0)
       ((dx-dx0)/rdx)^2 + ((dy-dy0)/rdy)^2 <= 1
     and we optimize rdx and rdy (and optionally dx0/dy0 if you enable it).
  ✅ Random and Optuna samplers both operate on a discrete grid (default 0.01).
     - random: sample then snap to grid
     - optuna: uses suggest_float(..., step=grid) (Optuna <=3 style)

Outputs:
  - best_cuts JSON (includes rdx/rdy, W2/eHCAL/coin_time ranges)
  - QA plots for dx fit + residual/pull + variable distributions
  - 1D scans around the final best point for (rdx,rdy,W2L,W2H,eL,tL,tH)

Example:
  python3 -u optimize_cuts_dxdy_ellipse_opt.py \
    --data DATA.root --sim SIM_elastic.root --bkg SIM_inelastic.root \
    --method optuna --n-trials 10000 \
    --grid-step 0.01 --cut-step 0.01 \
    --tL-min 113 --tL-max 118 --tH-min 122 --tH-max 127 \
    --out best_cuts_GEN3.json --qa-dir qa_GEN3_best

Notes:
  - This script assumes DATA has branch "helicity" with values +/-1 (or generally sign indicates helicity).
  - For SIM templates, fnucl is used to split neutron/proton (defaults fnucl_n=0, fnucl_p=1).
"""

import argparse
import json
import math
import os
import sys
import time
from dataclasses import dataclass
from typing import Dict, Optional, Tuple

import numpy as np
import uproot

# Optional deps
try:
    import optuna
except Exception:
    optuna = None

try:
    from scipy.optimize import nnls as scipy_nnls
except Exception:
    scipy_nnls = None

import matplotlib.pyplot as plt


# -----------------------------
# Data structures
# -----------------------------
@dataclass
class CutParams:
    # ellipse parameters (dx/dy)
    rdx: float
    rdy: float
    dx0: float
    dy0: float
    # other cuts
    W2L: float
    W2H: float
    eL: float
    tL: float
    tH: float

    def as_dict(self) -> Dict[str, float]:
        return dict(rdx=self.rdx, rdy=self.rdy, dx0=self.dx0, dy0=self.dy0,
                    W2L=self.W2L, W2H=self.W2H, eL=self.eL, tL=self.tL, tH=self.tH)


@dataclass
class ObjectiveConfig:
    dx_window: Tuple[float, float]
    min_data_events: int
    min_template_events: int
    rdx_min: float
    rdy_min: float
    W2_min_width: float
    t_min_width: float
    penalty_tight: float


# -----------------------------
# Progress bar
# -----------------------------
def print_progress(i: int, n: int, t0: float, best_sigma: float, prefix: str = "", bar_width: int = 30, every: int = 50):
    if n <= 0:
        return
    if (i % every) != 0 and (i + 1) != n:
        return

    frac = (i + 1) / n
    filled = int(bar_width * frac)
    bar = "#" * filled + "-" * (bar_width - filled)

    elapsed = time.time() - t0
    rate = (i + 1) / elapsed if elapsed > 0 else 0.0
    eta = (n - (i + 1)) / rate if rate > 0 else float("inf")

    def fmt_time(x):
        if not np.isfinite(x):
            return "?:??"
        m = int(x // 60)
        s = int(x % 60)
        return f"{m}:{s:02d}"

    msg = (
        f"\r{prefix}[{bar}] {100*frac:6.2f}% "
        f"{i+1}/{n}  "
        f"elapsed {fmt_time(elapsed)}  "
        f"eta {fmt_time(eta)}  "
        f"best sigma_phys {best_sigma:.6g}"
    )
    sys.stdout.write(msg)
    sys.stdout.flush()
    if (i + 1) == n:
        sys.stdout.write("\n")
        sys.stdout.flush()


# -----------------------------
# Numeric helpers
# -----------------------------
def hist1d(x: np.ndarray, w: np.ndarray, bins: np.ndarray) -> np.ndarray:
    h, _ = np.histogram(x, bins=bins, weights=w)
    return h.astype(np.float64)


def nnls_fit(data_h: np.ndarray, templates: np.ndarray) -> np.ndarray:
    if scipy_nnls is not None:
        coeff, _ = scipy_nnls(templates, data_h)
        return coeff
    coeff = np.linalg.lstsq(templates, data_h, rcond=None)[0]
    return np.maximum(coeff, 0.0)


def snap(x: float, step: float) -> float:
    return round(x / step) * step


def snap_params(c: CutParams, step: float) -> CutParams:
    return CutParams(
        rdx=max(step, snap(c.rdx, step)),
        rdy=max(step, snap(c.rdy, step)),
        dx0=snap(c.dx0, step),
        dy0=snap(c.dy0, step),
        W2L=snap(c.W2L, step),
        W2H=snap(c.W2H, step),
        eL=snap(c.eL, step),
        tL=snap(c.tL, step),
        tH=snap(c.tH, step),
    )


# -----------------------------
# Preselection (kept as in v2)
# -----------------------------
def preselect_cpp_like_data(vz, ePS, eSH, trP, eHCAL, helicity) -> np.ndarray:
    reject = (np.abs(vz) > 0.27) & (ePS < 0.2) & (np.abs((eSH + ePS) / trP - 1.0) > 0.2) & (eHCAL < 0.025) & (np.abs(helicity) != 1)
    return ~reject


def preselect_cpp_like_sim(vz, ePS, eSH, trP, eHCAL) -> np.ndarray:
    reject = (np.abs(vz) > 0.27) & (ePS < 0.2) & (np.abs((eSH + ePS) / trP - 1.0) > 0.2) & (eHCAL < 0.025)
    return ~reject


# -----------------------------
# Cut masks
# -----------------------------
def ellipse_mask(dx: np.ndarray, dy: np.ndarray, c: CutParams) -> np.ndarray:
    # avoid divide-by-zero
    rdx = max(c.rdx, 1e-9)
    rdy = max(c.rdy, 1e-9)
    u = (dx - c.dx0) / rdx
    v = (dy - c.dy0) / rdy
    return (u*u + v*v) <= 1.0


def build_cut_mask(dx, dy, W2, eHCAL, ct, cuts: CutParams, apply_ct_cut: bool) -> np.ndarray:
    m = ellipse_mask(dx, dy, cuts) & (W2 > cuts.W2L) & (W2 < cuts.W2H) & (eHCAL > cuts.eL)
    if apply_ct_cut:
        m &= (ct > cuts.tL) & (ct < cuts.tH)
    return m


# -----------------------------
# Asymmetry error helper
# -----------------------------
def asymmetry_and_error_from_counts(Np: float, Nm: float):
    """
    Raw helicity asymmetry and counting-statistical error:
      A = (N+ - N-) / (N+ + N-)
      sigma_A = 2*sqrt(N+*N-) / (N+ + N-)^(3/2)
    """
    Ntot = Np + Nm
    if Ntot <= 0:
        return 0.0, float("inf")
    A = (Np - Nm) / Ntot
    if Np <= 0 or Nm <= 0:
        return float(A), float("inf")
    sigma = 2.0 * math.sqrt(Np * Nm) / (Ntot ** 1.5)
    return float(A), float(sigma)


# -----------------------------
# Objective eval
# -----------------------------
def evaluate(
    bins: np.ndarray,
    # DATA
    dxD, dyD, W2D, ctD, eD, helD, wD,
    # neutron template
    dxN, dyN, W2N, ctN, eN, wN,
    # proton template
    dxP, dyP, W2P, ctP, eP, wP,
    # optional bkg template
    bkg: Optional[Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]],
    cuts: CutParams,
    cfg: ObjectiveConfig,
):
    """
    Returns (score, details) where score = -sigma_phys (maximize score).
    """
    # basic validity
    if cuts.W2L >= cuts.W2H or cuts.tL >= cuts.tH:
        return -np.inf, None
    if cuts.rdx < cfg.rdx_min or cuts.rdy < cfg.rdy_min:
        return -np.inf, None
    if (cuts.W2H - cuts.W2L) < cfg.W2_min_width:
        return -np.inf, None
    if (cuts.tH - cuts.tL) < cfg.t_min_width:
        return -np.inf, None

    # DATA uses coin_time cut; templates DO NOT
    mD = build_cut_mask(dxD, dyD, W2D, eD, ctD, cuts, apply_ct_cut=True)
    mN = build_cut_mask(dxN, dyN, W2N, eN, ctN, cuts, apply_ct_cut=False)
    mP = build_cut_mask(dxP, dyP, W2P, eP, ctP, cuts, apply_ct_cut=False)

    nD = int(mD.sum())
    nN = int(mN.sum())
    nP = int(mP.sum())

    if nD < cfg.min_data_events:
        return -np.inf, None
    if nN < cfg.min_template_events or nP < cfg.min_template_events:
        return -np.inf, None

    # raw asymmetry counting error from DATA after cuts
    hel_sel = helD[mD]
    Nplus = float(np.sum(hel_sel > 0))
    Nminus = float(np.sum(hel_sel < 0))
    Araw, sigma_Araw = asymmetry_and_error_from_counts(Nplus, Nminus)
    if not np.isfinite(sigma_Araw):
        return -np.inf, None

    # dx fit (unit-area PDFs)
    hD = hist1d(dxD[mD], wD[mD], bins)
    hN = hist1d(dxN[mN], wN[mN], bins)
    hP = hist1d(dxP[mP], wP[mP], bins)

    sumN = float(hN.sum())
    sumP = float(hP.sum())
    if sumN <= 0 or sumP <= 0:
        return -np.inf, None

    n_pdf = hN / sumN
    p_pdf = hP / sumP

    has_bkg = bkg is not None
    b_pdf = None
    nB = 0
    if has_bkg:
        dxB, dyB, W2B, ctB, eB, wB = bkg
        mB = build_cut_mask(dxB, dyB, W2B, eB, ctB, cuts, apply_ct_cut=False)
        nB = int(mB.sum())
        if nB < cfg.min_template_events:
            return -np.inf, None
        hB = hist1d(dxB[mB], wB[mB], bins)
        sumB = float(hB.sum())
        if sumB <= 0:
            return -np.inf, None
        b_pdf = hB / sumB

    templates = [p_pdf, n_pdf] + ([b_pdf] if has_bkg else [])
    T = np.vstack(templates).T
    coeff = nnls_fit(hD, T)

    cP = float(coeff[0])
    cN = float(coeff[1])
    cB = float(coeff[2]) if has_bkg else 0.0
    if cP <= 0:
        return -np.inf, None

    N = cP
    R = cN / cP
    Nbg = (cB / cP) if has_bkg else 0.0

    modelP = N * p_pdf
    modelN = (N * R) * n_pdf
    modelB = (N * Nbg) * b_pdf if has_bkg else (0.0 * p_pdf)
    modelTot = modelP + modelN + (modelB if has_bkg else 0.0)

    lo, hi = cfg.dx_window
    centers = 0.5 * (bins[:-1] + bins[1:])
    win = (centers >= lo) & (centers <= hi)

    Nn_win = float(modelN[win].sum())
    Ntot_win = float(modelTot[win].sum())
    if Ntot_win <= 0 or Nn_win <= 0:
        return -np.inf, None

    fn = Nn_win / Ntot_win
    if fn <= 0:
        return -np.inf, None

    sigma_phys = float(sigma_Araw / fn)

    if cfg.penalty_tight > 0:
        # optional penalty to avoid absurdly-tight cuts
        inv_width = (
            1.0 / max(cuts.rdx, 1e-9) +
            1.0 / max(cuts.rdy, 1e-9) +
            1.0 / max(cuts.W2H - cuts.W2L, 1e-9) +
            1.0 / max(cuts.tH - cuts.tL, 1e-9)
        )
        sigma_phys = sigma_phys + cfg.penalty_tight * inv_width

    score = -sigma_phys

    details = dict(
        score=float(score),
        sigma_phys=float(sigma_phys),
        Araw=float(Araw),
        sigma_Araw=float(sigma_Araw),
        Nplus=float(Nplus),
        Nminus=float(Nminus),
        fn=float(fn),
        cuts=cuts.as_dict(),
        params=dict(N=N, R=R, Nbg=Nbg),
        coeff=dict(cP=cP, cN=cN, cB=cB),
        counts=dict(nData=nD, nN=nN, nP=nP, nB=nB),
        yields=dict(Nn_window=Nn_win, Ntot_window=Ntot_win),
        hists=dict(
            hD=hD,
            p_pdf=p_pdf, n_pdf=n_pdf, b_pdf=b_pdf,
            modelP=modelP, modelN=modelN, modelB=modelB, hTot=modelTot
        ),
        bins=bins,
        centers=centers,
    )
    return float(score), details


# -----------------------------
# Chunked loading
# -----------------------------
def iterate_filtered(
    rootfile: str,
    treename: str,
    branches: Tuple[str, ...],
    keep_fn,
    step_size: str,
    rng: np.random.Generator,
    downsample: float,
    max_keep: Optional[int],
    dtype: np.dtype,
) -> Dict[str, np.ndarray]:
    out = {br: [] for br in branches}
    kept = 0
    t0 = time.time()
    chunk_i = 0

    for chunk in uproot.iterate(f"{rootfile}:{treename}", list(branches), library="np", step_size=step_size):
        chunk_i += 1
        keep = keep_fn(chunk)
        if downsample < 1.0:
            keep &= (rng.random(keep.shape[0]) < downsample)

        n_keep = int(keep.sum())
        if n_keep == 0:
            continue

        if max_keep is not None and kept + n_keep > max_keep:
            idx = np.flatnonzero(keep)
            need = max_keep - kept
            idx = idx[:need]
            for br in branches:
                out[br].append(chunk[br][idx].astype(dtype, copy=False))
            kept += len(idx)
            break
        else:
            for br in branches:
                out[br].append(chunk[br][keep].astype(dtype, copy=False))
            kept += n_keep

        if chunk_i % 10 == 0:
            sys.stdout.write(f"\rLoading {os.path.basename(rootfile)}: kept {kept} rows ...")
            sys.stdout.flush()

    sys.stdout.write(f"\rLoading {os.path.basename(rootfile)}: kept {kept} rows. Done in {time.time()-t0:.1f}s\n")
    sys.stdout.flush()

    out2 = {}
    for br, pieces in out.items():
        out2[br] = np.concatenate(pieces) if len(pieces) else np.array([], dtype=dtype)
    return out2


# -----------------------------
# QA plots
# -----------------------------
def ensure_dir(d: str):
    if d and not os.path.isdir(d):
        os.makedirs(d, exist_ok=True)


def save_fig(path_png: str, path_pdf: str):
    plt.savefig(path_png, dpi=160, bbox_inches="tight")
    plt.savefig(path_pdf, bbox_inches="tight")
    plt.close()


def qa_dx_fit(details: dict, qa_dir: str, dx_window: Tuple[float, float]):
    ensure_dir(qa_dir)
    hD = details["hists"]["hD"]
    modelP = details["hists"]["modelP"]
    modelN = details["hists"]["modelN"]
    modelB = details["hists"]["modelB"]
    hTot = details["hists"]["hTot"]
    bins = details["bins"]
    centers = details["centers"]
    lo, hi = dx_window

    N = details["params"]["N"]
    R = details["params"]["R"]
    Nbg = details["params"]["Nbg"]

    plt.figure()
    yerr = np.sqrt(np.maximum(hD, 1.0))
    plt.errorbar(centers, hD, yerr=yerr, fmt=".", label="Data")
    plt.step(bins[:-1], hTot, where="post", label="Fit total")
    plt.step(bins[:-1], modelP, where="post", label=f"N*p (proton), N={N:.3g}")
    plt.step(bins[:-1], modelN, where="post", label=f"N*R*n (neutron), R={R:.3g}")
    if details["hists"].get("b_pdf", None) is not None:
        plt.step(bins[:-1], modelB, where="post", label=f"N*Nbg*b (bkg), Nbg={Nbg:.3g}")
    plt.axvspan(lo, hi, alpha=0.15, label=f"dx window [{lo},{hi}]")
    plt.xlabel("dx")
    plt.ylabel("Counts")
    plt.title(f"dx fit @ best cuts | sigma_phys={details['sigma_phys']:.4g} | fn={details['fn']:.3f}")
    plt.legend()
    save_fig(os.path.join(qa_dir, "qa_dx_fit_best.png"),
             os.path.join(qa_dir, "qa_dx_fit_best.pdf"))



def qa_W2_fit(
    qa_dir: str,
    W2_bins: np.ndarray,
    # data
    W2D: np.ndarray, wD: np.ndarray, mD: np.ndarray,
    # neutron template
    W2N: np.ndarray, wN: np.ndarray, mN: np.ndarray,
    # proton template
    W2P: np.ndarray, wP: np.ndarray, mP: np.ndarray,
    # optional bkg template (full tuple) + its mask
    bkg_tuple,
    mB: Optional[np.ndarray] = None,
    title_extra: str = ""
):
    """Template fit of W2 after final cuts (NNLS on unit-area PDFs)."""
    ensure_dir(qa_dir)

    hD = hist1d(W2D[mD], wD[mD], W2_bins)
    yerr = np.sqrt(np.maximum(hD, 1.0))

    hN = hist1d(W2N[mN], wN[mN], W2_bins)
    hP = hist1d(W2P[mP], wP[mP], W2_bins)

    sumN = float(hN.sum()); sumP = float(hP.sum()); sumD = float(hD.sum())
    if sumD <= 0 or sumN <= 0 or sumP <= 0:
        return None

    n_pdf = hN / sumN
    p_pdf = hP / sumP

    has_bkg = bkg_tuple is not None
    b_pdf = None
    if has_bkg:
        dxB, dyB, W2B, ctB, eB, wB = bkg_tuple
        if mB is None:
            mB = np.ones_like(W2B, dtype=bool)
        hB = hist1d(W2B[mB], wB[mB], W2_bins)
        sumB = float(hB.sum())
        if sumB > 0:
            b_pdf = hB / sumB
        else:
            has_bkg = False

    templates = [p_pdf, n_pdf] + ([b_pdf] if has_bkg else [])
    T = np.vstack(templates).T
    coeff = nnls_fit(hD, T)

    cP = float(coeff[0])
    cN = float(coeff[1])
    cB = float(coeff[2]) if has_bkg else 0.0

    modelP = cP * p_pdf
    modelN = cN * n_pdf
    modelB = cB * b_pdf if has_bkg else (0.0 * p_pdf)
    modelTot = modelP + modelN + (modelB if has_bkg else 0.0)

    centers = 0.5 * (W2_bins[:-1] + W2_bins[1:])

    # --- Plot ---
    plt.figure()
    plt.errorbar(centers, hD, yerr=yerr, fmt=".", label="Data")
    plt.step(W2_bins[:-1], modelTot, where="post", label="Fit total")
    plt.step(W2_bins[:-1], modelP, where="post", label=f"p template, cP={cP:.3g}")
    plt.step(W2_bins[:-1], modelN, where="post", label=f"n template, cN={cN:.3g}")
    if has_bkg:
        plt.step(W2_bins[:-1], modelB, where="post", label=f"bkg template, cB={cB:.3g}")
    plt.xlabel("W2")
    plt.ylabel("Counts")
    plt.title("W2 template fit (after best cuts)" + (f" | {title_extra}" if title_extra else ""))
    plt.legend()
    save_fig(os.path.join(qa_dir, "qa_W2_fit_best.png"),
             os.path.join(qa_dir, "qa_W2_fit_best.pdf"))

    # --- Residual/pull ---
    resid = hD - modelTot
    pull = resid / np.sqrt(np.maximum(hD, 1.0))

    plt.figure()
    plt.plot(centers, resid, ".", label="Data - Fit")
    plt.axhline(0.0)
    plt.xlabel("W2")
    plt.ylabel("Residual")
    plt.title("W2 residual (Data - Fit)")
    plt.legend()
    save_fig(os.path.join(qa_dir, "qa_W2_residual_best.png"),
             os.path.join(qa_dir, "qa_W2_residual_best.pdf"))

    plt.figure()
    plt.plot(centers, pull, ".", label="Pull")
    plt.axhline(0.0)
    plt.axhline(3.0, linestyle="--")
    plt.axhline(-3.0, linestyle="--")
    plt.xlabel("W2")
    plt.ylabel("(Data - Fit)/sqrt(Data)")
    plt.title("W2 pull")
    plt.legend()
    save_fig(os.path.join(qa_dir, "qa_W2_pull_best.png"),
             os.path.join(qa_dir, "qa_W2_pull_best.pdf"))

    return dict(
        coeff=dict(cP=cP, cN=cN, cB=cB),
        counts=dict(
            data=float(hD.sum()),
            model_p=float(modelP.sum()),
            model_n=float(modelN.sum()),
            model_b=float(modelB.sum()) if has_bkg else 0.0,
        ),
    )

def qa_dx_residual(details: dict, qa_dir: str):
    ensure_dir(qa_dir)
    hD = details["hists"]["hD"]
    hTot = details["hists"]["hTot"]
    centers = details["centers"]

    resid = hD - hTot
    pull = resid / np.sqrt(np.maximum(hD, 1.0))

    plt.figure()
    plt.plot(centers, resid, ".", label="Data - Fit")
    plt.axhline(0.0)
    plt.xlabel("dx")
    plt.ylabel("Residual")
    plt.title("dx residual (Data - Fit)")
    plt.legend()
    save_fig(os.path.join(qa_dir, "qa_dx_residual_best.png"),
             os.path.join(qa_dir, "qa_dx_residual_best.pdf"))

    plt.figure()
    plt.plot(centers, pull, ".", label="Pull")
    plt.axhline(0.0)
    plt.axhline(3.0, linestyle="--")
    plt.axhline(-3.0, linestyle="--")
    plt.xlabel("dx")
    plt.ylabel("(Data - Fit)/sqrt(Data)")
    plt.title("dx pull")
    plt.legend()
    save_fig(os.path.join(qa_dir, "qa_dx_pull_best.png"),
             os.path.join(qa_dir, "qa_dx_pull_best.pdf"))


def qa_data_vars_before_after(dxD, dyD, W2D, ctD, eHD, mask_best, qa_dir: str):
    ensure_dir(qa_dir)
    vars_ = [
        ("dx", dxD),
        ("dy", dyD),
        ("W2", W2D),
        ("eHCAL", eHD),
        ("coin_time", ctD),
    ]

    plt.figure(figsize=(12, 8))
    for i, (name, arr) in enumerate(vars_, 1):
        plt.subplot(2, 3, i)
        lo = np.nanpercentile(arr, 1)
        hi = np.nanpercentile(arr, 99)
        if not np.isfinite(lo) or not np.isfinite(hi) or lo == hi:
            lo, hi = float(np.min(arr)), float(np.max(arr))
        bins = 120
        plt.hist(arr, bins=bins, range=(lo, hi), histtype="step", label="before")
        plt.hist(arr[mask_best], bins=bins, range=(lo, hi), histtype="step", label="after")
        plt.xlabel(name)
        plt.ylabel("Counts")
        plt.legend()

    plt.tight_layout()
    save_fig(os.path.join(qa_dir, "qa_data_vars_before_after.png"),
             os.path.join(qa_dir, "qa_data_vars_before_after.pdf"))


def scan_1d_sigma(eval_one, best: CutParams, args, qa_dir: str):
    ensure_dir(qa_dir)

    params = [
        ("rdx", lambda c: c.rdx, lambda c,v: CutParams(v, c.rdy, c.dx0, c.dy0, c.W2L, c.W2H, c.eL, c.tL, c.tH), args.rdx_min, args.rdx_max),
        ("rdy", lambda c: c.rdy, lambda c,v: CutParams(c.rdx, v, c.dx0, c.dy0, c.W2L, c.W2H, c.eL, c.tL, c.tH), args.rdy_min, args.rdy_max),
        ("W2L", lambda c: c.W2L, lambda c,v: CutParams(c.rdx, c.rdy, c.dx0, c.dy0, v, c.W2H, c.eL, c.tL, c.tH), args.W2L_min, args.W2L_max),
        ("W2H", lambda c: c.W2H, lambda c,v: CutParams(c.rdx, c.rdy, c.dx0, c.dy0, c.W2L, v, c.eL, c.tL, c.tH), args.W2H_min, args.W2H_max),
        ("eL",  lambda c: c.eL,  lambda c,v: CutParams(c.rdx, c.rdy, c.dx0, c.dy0, c.W2L, c.W2H, v, c.tL, c.tH),  args.eL_min,  args.eL_max),
        ("tL",  lambda c: c.tL,  lambda c,v: CutParams(c.rdx, c.rdy, c.dx0, c.dy0, c.W2L, c.W2H, c.eL, v, c.tH),  args.tL_min,  args.tL_max),
        ("tH",  lambda c: c.tH,  lambda c,v: CutParams(c.rdx, c.rdy, c.dx0, c.dy0, c.W2L, c.W2H, c.eL, c.tL, v),  args.tH_min,  args.tH_max),
    ]

    step = args.grid_step
    for name, getv, setv, vmin, vmax in params:
        xs = np.arange(snap(vmin, step), snap(vmax, step) + 0.5*step, step)
        ys = []
        for x in xs:
            c = setv(best, float(x))
            c = snap_params(c, step)
            # keep ordering constraints
            if not (c.W2L < c.W2H and c.tL < c.tH and c.rdx > 0 and c.rdy > 0):
                ys.append(np.nan); continue
            score, det = eval_one(c)
            if det is None or not np.isfinite(det["sigma_phys"]):
                ys.append(np.nan)
            else:
                ys.append(det["sigma_phys"])

        plt.figure()
        plt.plot(xs, ys, ".-")
        plt.axvline(getv(best), linestyle="--", label="best")
        plt.xlabel(name)
        plt.ylabel("sigma_phys = sigma_Araw / fn")
        plt.title(f"1D scan: {name} (others fixed)")
        plt.legend()
        save_fig(os.path.join(qa_dir, f"scan1d_{name}.png"),
                 os.path.join(qa_dir, f"scan1d_{name}.pdf"))


# -----------------------------
# CLI
# -----------------------------
def parse_args():
    ap = argparse.ArgumentParser(description="Optimize QE cuts using dx template fit and asymmetry statistical error, with ellipse dx/dy cut.")

    ap.add_argument("--data", required=True)
    ap.add_argument("--sim", required=True)
    ap.add_argument("--bkg", default="")
    ap.add_argument("--tree", default="Tout")

    ap.add_argument("--br-vz", default="vz")
    ap.add_argument("--br-ePS", default="ePS")
    ap.add_argument("--br-eSH", default="eSH")
    ap.add_argument("--br-trP", default="trP")
    ap.add_argument("--br-dx", default="dx")
    ap.add_argument("--br-dy", default="dy")
    ap.add_argument("--br-W2", default="W2")
    ap.add_argument("--br-ct", default="coin_time")
    ap.add_argument("--br-eHCAL", default="eHCAL")
    ap.add_argument("--br-helicity", default="helicity")

    ap.add_argument("--br-weight", default="weight")
    ap.add_argument("--br-fnucl", default="fnucl")
    ap.add_argument("--fnucl-n", type=float, default=0.0)
    ap.add_argument("--fnucl-p", type=float, default=1.0)

    ap.add_argument("--dx-min", type=float, default=-4.0)
    ap.add_argument("--dx-max", type=float, default=4.0)
    ap.add_argument("--dx-bins", type=int, default=100)
    ap.add_argument("--dx-win-low", type=float, default=-0.4)
    ap.add_argument("--dx-win-high", type=float, default=0.4)

    ap.add_argument("--method", choices=["random", "optuna"], default="random")
    ap.add_argument("--n-trials", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=1)

    # Ellipse radii bounds (optimize these)
    ap.add_argument("--rdx-min", type=float, default=0.30)
    ap.add_argument("--rdx-max", type=float, default=0.70)
    ap.add_argument("--rdy-min", type=float, default=0.30)
    ap.add_argument("--rdy-max", type=float, default=0.70)

    # Center (fixed by default at 0,0; can allow scanning with --opt-center)
    ap.add_argument("--dx0", type=float, default=0.0)
    ap.add_argument("--dy0", type=float, default=0.0)
    ap.add_argument("--opt-center", action="store_true", help="Also optimize dx0 and dy0 on the grid within bounds below.")
    ap.add_argument("--dx0-min", type=float, default=-0.20)
    ap.add_argument("--dx0-max", type=float, default=0.20)
    ap.add_argument("--dy0-min", type=float, default=-0.20)
    ap.add_argument("--dy0-max", type=float, default=0.20)

    # W2 window
    ap.add_argument("--W2L-min", type=float, default=-2.0)
    ap.add_argument("--W2L-max", type=float, default=0.0)
    ap.add_argument("--W2H-min", type=float, default=1.2)
    ap.add_argument("--W2H-max", type=float, default=2.0)

    # W2 QA/fit histogram (after best cuts)
    ap.add_argument("--W2-plot-min", type=float, default=-2.0, help="W2 histogram min for QA fit plot")
    ap.add_argument("--W2-plot-max", type=float, default=4.0, help="W2 histogram max for QA fit plot")
    ap.add_argument("--W2-plot-bins", type=int, default=120, help="W2 histogram bins for QA fit plot")

    # eHCAL threshold
    ap.add_argument("--eL-min", type=float, default=0.025)
    ap.add_argument("--eL-max", type=float, default=0.35)

    # DATA only coin-time window
    ap.add_argument("--tL-min", type=float, default=170.0)
    ap.add_argument("--tL-max", type=float, default=190.0)
    ap.add_argument("--tH-min", type=float, default=170.0)
    ap.add_argument("--tH-max", type=float, default=190.0)

    ap.add_argument("--min-data-events", type=int, default=200)
    ap.add_argument("--min-template-events", type=int, default=200)
    ap.add_argument("--W2-min-width", type=float, default=0.30)
    ap.add_argument("--t-min-width", type=float, default=3.0)
    ap.add_argument("--penalty-tight", type=float, default=0.0)

    ap.add_argument("--step-size", default="200 MB")
    ap.add_argument("--dtype", choices=["float32", "float64"], default="float32")
    ap.add_argument("--max-data", type=int, default=0)
    ap.add_argument("--max-sim", type=int, default=0)
    ap.add_argument("--max-bkg", type=int, default=0)
    ap.add_argument("--downsample-data", type=float, default=1.0)
    ap.add_argument("--downsample-sim", type=float, default=1.0)
    ap.add_argument("--downsample-bkg", type=float, default=1.0)

    # grid control
    ap.add_argument("--grid-step", type=float, default=0.01, help="Resolution used for BOTH random and optuna proposals (and for 1D scans).")

    ap.add_argument("--out", default="best_cuts.json")
    ap.add_argument("--qa-dir", default="qa_best")
    ap.add_argument("--no-1d-scan", action="store_true", help="Disable 1D scan plots.")
    return ap.parse_args()


# -----------------------------
# Main
# -----------------------------
def main():
    args = parse_args()
    rng = np.random.default_rng(args.seed)
    dtype = np.float32 if args.dtype == "float32" else np.float64

    bins = np.linspace(args.dx_min, args.dx_max, args.dx_bins + 1)

    # -------- load DATA --------
    data_branches = (args.br_dx, args.br_dy, args.br_W2, args.br_ct, args.br_eHCAL,
                     args.br_vz, args.br_ePS, args.br_eSH, args.br_trP, args.br_helicity)

    def keep_data(chunk):
        return preselect_cpp_like_data(chunk[args.br_vz], chunk[args.br_ePS], chunk[args.br_eSH],
                                       chunk[args.br_trP], chunk[args.br_eHCAL], chunk[args.br_helicity])

    max_data = args.max_data if args.max_data > 0 else None
    D = iterate_filtered(args.data, args.tree, data_branches, keep_data, args.step_size, rng, args.downsample_data, max_data, dtype)

    dxD = D[args.br_dx]; dyD = D[args.br_dy]; W2D = D[args.br_W2]; ctD = D[args.br_ct]
    eHD = D[args.br_eHCAL]; helD = D[args.br_helicity]
    wD = np.ones_like(dxD, dtype=np.float64)

    # -------- load SIM --------
    sim_branches = (args.br_dx, args.br_dy, args.br_W2, args.br_eHCAL,
                    args.br_vz, args.br_ePS, args.br_eSH, args.br_trP,
                    args.br_weight, args.br_fnucl)

    def keep_sim(chunk):
        return preselect_cpp_like_sim(chunk[args.br_vz], chunk[args.br_ePS], chunk[args.br_eSH],
                                      chunk[args.br_trP], chunk[args.br_eHCAL])

    max_sim = args.max_sim if args.max_sim > 0 else None
    S = iterate_filtered(args.sim, args.tree, sim_branches, keep_sim, args.step_size, rng, args.downsample_sim, max_sim, dtype)

    dxS = S[args.br_dx]; dyS = S[args.br_dy]; W2S = S[args.br_W2]
    eHS = S[args.br_eHCAL]; wS = S[args.br_weight]; fnucl = S[args.br_fnucl]
    ctS = np.zeros_like(dxS, dtype=dtype)

    ok_np = (fnucl == args.fnucl_n) | (fnucl == args.fnucl_p)
    dxS, dyS, W2S, ctS, eHS, wS, fnucl = dxS[ok_np], dyS[ok_np], W2S[ok_np], ctS[ok_np], eHS[ok_np], wS[ok_np], fnucl[ok_np]

    isN = (fnucl == args.fnucl_n)
    isP = (fnucl == args.fnucl_p)
    if isN.sum() == 0 or isP.sum() == 0:
        raise RuntimeError(f"Sim split produced empty sample(s): n={isN.sum()} p={isP.sum()} (check fnucl values).")

    dxN, dyN, W2N, ctN, eHN, wN = dxS[isN], dyS[isN], W2S[isN], ctS[isN], eHS[isN], wS[isN]
    dxP, dyP, W2P, ctP, eHP, wP = dxS[isP], dyS[isP], W2S[isP], ctS[isP], eHS[isP], wS[isP]

    # -------- load BKG (optional) --------
    bkg_tuple = None
    if args.bkg:
        bkg_branches = (args.br_dx, args.br_dy, args.br_W2, args.br_eHCAL,
                        args.br_vz, args.br_ePS, args.br_eSH, args.br_trP,
                        args.br_weight)

        def keep_bkg(chunk):
            return preselect_cpp_like_sim(chunk[args.br_vz], chunk[args.br_ePS], chunk[args.br_eSH],
                                          chunk[args.br_trP], chunk[args.br_eHCAL])

        max_bkg = args.max_bkg if args.max_bkg > 0 else None
        B = iterate_filtered(args.bkg, args.tree, bkg_branches, keep_bkg, args.step_size, rng, args.downsample_bkg, max_bkg, dtype)

        dxB = B[args.br_dx]; dyB = B[args.br_dy]; W2B = B[args.br_W2]; eHB = B[args.br_eHCAL]
        wB = B[args.br_weight]
        ctB = np.zeros_like(dxB, dtype=dtype)
        bkg_tuple = (dxB, dyB, W2B, ctB, eHB, wB)

    print("Finished loading + preselection. Starting optimization...")

    cfg = ObjectiveConfig(
        dx_window=(args.dx_win_low, args.dx_win_high),
        min_data_events=args.min_data_events,
        min_template_events=args.min_template_events,
        rdx_min=args.rdx_min,
        rdy_min=args.rdy_min,
        W2_min_width=args.W2_min_width,
        t_min_width=args.t_min_width,
        penalty_tight=args.penalty_tight,
    )

    def eval_one(cuts: CutParams):
        return evaluate(
            bins,
            dxD, dyD, W2D, ctD, eHD, helD, wD,
            dxN, dyN, W2N, ctN, eHN, wN,
            dxP, dyP, W2P, ctP, eHP, wP,
            bkg_tuple,
            cuts,
            cfg,
        )

    best_score = -np.inf
    best_cuts: Optional[CutParams] = None
    best_details = None

    step = args.grid_step

    if args.method == "random":
        t0 = time.time()
        for i in range(args.n_trials):
            dx0 = float(rng.uniform(args.dx0_min, args.dx0_max)) if args.opt_center else float(args.dx0)
            dy0 = float(rng.uniform(args.dy0_min, args.dy0_max)) if args.opt_center else float(args.dy0)

            cuts = CutParams(
                rdx=float(rng.uniform(args.rdx_min, args.rdx_max)),
                rdy=float(rng.uniform(args.rdy_min, args.rdy_max)),
                dx0=dx0,
                dy0=dy0,
                W2L=float(rng.uniform(args.W2L_min, args.W2L_max)),
                W2H=float(rng.uniform(args.W2H_min, args.W2H_max)),
                eL=float(rng.uniform(args.eL_min, args.eL_max)),
                tL=float(rng.uniform(args.tL_min, args.tL_max)),
                tH=float(rng.uniform(args.tH_min, args.tH_max)),
            )
            cuts = snap_params(cuts, step)
            score, det = eval_one(cuts)
            if score > best_score:
                best_score = score
                best_cuts = cuts
                best_details = det
                print(f"\n[trial {i+1}/{args.n_trials}] NEW BEST sigma_phys={det['sigma_phys']:.6g} cuts={best_cuts.as_dict()} fn={det['fn']:.4g} sigma_Araw={det['sigma_Araw']:.4g}")
            best_sigma = (best_details["sigma_phys"] if best_details is not None else float("inf"))
            print_progress(i, args.n_trials, t0, best_sigma, prefix="Optimize ", every=50)

    else:
        if optuna is None:
            raise RuntimeError("optuna not installed. pip install optuna")
        if scipy_nnls is None:
            print("Warning: scipy not found; NNLS uses fallback. Recommended: pip install scipy")

        def objective(trial):
            if args.opt_center:
                dx0 = trial.suggest_float("dx0", args.dx0_min, args.dx0_max, step=step)
                dy0 = trial.suggest_float("dy0", args.dy0_min, args.dy0_max, step=step)
            else:
                dx0 = float(args.dx0)
                dy0 = float(args.dy0)

            cuts = CutParams(
                rdx=trial.suggest_float("rdx", args.rdx_min, args.rdx_max, step=step),
                rdy=trial.suggest_float("rdy", args.rdy_min, args.rdy_max, step=step),
                dx0=dx0,
                dy0=dy0,
                W2L=trial.suggest_float("W2L", args.W2L_min, args.W2L_max, step=step),
                W2H=trial.suggest_float("W2H", args.W2H_min, args.W2H_max, step=step),
                eL=trial.suggest_float("eL", args.eL_min, args.eL_max, step=step),
                tL=trial.suggest_float("tL", args.tL_min, args.tL_max, step=step),
                tH=trial.suggest_float("tH", args.tH_min, args.tH_max, step=step),
            )
            score, _ = eval_one(cuts)
            return score

        study = optuna.create_study(direction="maximize", sampler=optuna.samplers.TPESampler(seed=args.seed))
        study.optimize(objective, n_trials=args.n_trials)

        # reconstruct best cut params
        bp = study.best_params
        best_cuts = CutParams(
            rdx=float(bp["rdx"]),
            rdy=float(bp["rdy"]),
            dx0=float(bp["dx0"]) if args.opt_center else float(args.dx0),
            dy0=float(bp["dy0"]) if args.opt_center else float(args.dy0),
            W2L=float(bp["W2L"]),
            W2H=float(bp["W2H"]),
            eL=float(bp["eL"]),
            tL=float(bp["tL"]),
            tH=float(bp["tH"]),
        )
        best_cuts = snap_params(best_cuts, step)
        best_score, best_details = eval_one(best_cuts)

    if best_cuts is None or best_details is None or not np.isfinite(best_details["sigma_phys"]):
        raise RuntimeError("No valid solution found. Try loosening bounds or lowering min event requirements.")

    final_cuts = best_cuts
    final_details = best_details

    # -----------------------------
    # QA + 1D scans
    # -----------------------------
    ensure_dir(args.qa_dir)
    mD_best = build_cut_mask(dxD, dyD, W2D, eHD, ctD, final_cuts, apply_ct_cut=True)

    qa_dx_fit(final_details, args.qa_dir, (args.dx_win_low, args.dx_win_high))
    qa_dx_residual(final_details, args.qa_dir)
    qa_data_vars_before_after(dxD, dyD, W2D, ctD, eHD, mD_best, args.qa_dir)

    # W2 distribution + template fit QA (after best cuts)
    W2_bins = np.linspace(args.W2_plot_min, args.W2_plot_max, args.W2_plot_bins + 1)
    mN_best = build_cut_mask(dxN, dyN, W2N, eHN, ctN, final_cuts, apply_ct_cut=False)
    mP_best = build_cut_mask(dxP, dyP, W2P, eHP, ctP, final_cuts, apply_ct_cut=False)
    mB_best = None
    if bkg_tuple is not None:
        dxB, dyB, W2B, ctB, eHB, wB = bkg_tuple
        mB_best = build_cut_mask(dxB, dyB, W2B, eHB, ctB, final_cuts, apply_ct_cut=False)

    w2_fit = qa_W2_fit(
        args.qa_dir,
        W2_bins,
        W2D, wD, mD_best,
        W2N, wN, mN_best,
        W2P, wP, mP_best,
        bkg_tuple,
        mB=mB_best,
        title_extra=f"sigma_phys={final_details['sigma_phys']:.4g}, fn(dx)={final_details['fn']:.3f}"
    )
    if w2_fit is not None:
        print("[QA] W2 fit coefficients:", w2_fit["coeff"])
        print("[QA] W2 fit component totals:", w2_fit["counts"])

        if not args.no_1d_scan:
            scan_1d_sigma(eval_one, final_cuts, args, args.qa_dir)

        with open(os.path.join(args.qa_dir, "qa_best_fit_summary.json"), "w") as f:
            json.dump(
                {
                    "best_sigma_phys": float(final_details["sigma_phys"]),
                    "best_cuts": final_cuts.as_dict(),
                    "fn": final_details["fn"],
                    "Araw": final_details["Araw"],
                    "sigma_Araw": final_details["sigma_Araw"],
                    "Nplus": final_details["Nplus"],
                    "Nminus": final_details["Nminus"],
                    "counts_after_cuts": final_details["counts"],
                    "yields_in_dx_window": final_details["yields"],
                    "dx_window": [args.dx_win_low, args.dx_win_high],
                    "coin_time_cut_applied_to": "data_only",
                    "grid_step": args.grid_step,
                    "opt_center": bool(args.opt_center),
                },
                f,
                indent=2,
                sort_keys=True,
            )

        out = {
            "best_sigma_phys": float(final_details["sigma_phys"]),
            "best_score": float(final_details["score"]),
            "best_cuts": final_cuts.as_dict(),
            "meta": {
                "tree": args.tree,
                "data": args.data,
                "sim": args.sim,
                "bkg": args.bkg,
                "dx_hist": {"min": args.dx_min, "max": args.dx_max, "bins": args.dx_bins},
                "dx_window": [args.dx_win_low, args.dx_win_high],
                "method": args.method,
                "n_trials": args.n_trials,
                "seed": args.seed,
                "step_size": args.step_size,
                "dtype": args.dtype,
                "downsample_data": args.downsample_data,
                "downsample_sim": args.downsample_sim,
                "downsample_bkg": args.downsample_bkg,
                "qa_dir": args.qa_dir,
                "grid_step": args.grid_step,
                "opt_center": bool(args.opt_center),
            },
        }
        with open(args.out, "w") as f:
            json.dump(out, f, indent=2, sort_keys=True)

        print("\n=== DONE ===")
        print("Final Best sigma_phys:", out["best_sigma_phys"])
        print("Final Best cuts:", out["best_cuts"])
        print("Wrote:", os.path.abspath(args.out))
        print("QA plots saved in:", os.path.abspath(args.qa_dir))


if __name__ == "__main__":
    main()
