#!/usr/bin/env python3
"""
optimize_cuts_from_cpp_stream_final.py

Integrated "final" version based on your optimize_cuts_from_cpp_stream.py, with:
  ✅ Chunked reading (uproot.iterate) to avoid OOM
  ✅ Same C++-like preselection logic
  ✅ coin_time cut applied to DATA only (SIM/BKG assumed in-time; SIM/BKG may not have coin_time)
  ✅ Random or Optuna optimization + in-terminal progress bar
  ✅ 0.05-level cuts:
        - snap best solution to a grid (default step=0.05)
        - local grid refinement around snapped best (±step) on the grid (3^7 = 2187 evals)
  ✅ QA plots for the final best solution:
        - dx template fit plot (data points + NNLS fit components + total)
        - dx residual and pull
        - data variable distributions before/after (dy, W2, eHCAL, coin_time)
        - JSON summary of fit coefficients/yields

FOM:
  For each candidate cuts C:
    - apply dy/W2/eHCAL to all samples
    - apply coin_time cut to DATA only
    - histogram dx
    - NNLS fit: Data ≈ aN*N + aP*P + aB*B (a>=0)
    - integrate in dx window:
        FOM = Nn / sqrt(Nn + Nbkg)

Example:
  python3 -u optimize_cuts_from_cpp_stream_final.py \
    --data DATA.root --sim SIM_np.root --bkg BKG.root \
    --downsample-data 0.005 --downsample-bkg 0.003 --downsample-sim 0.5 \
    --dtype float32 --step-size "150 MB" \
    --tL-min 113 --tL-max 118 --tH-min 122 --tH-max 127 \
    --method random --n-trials 10000 \
    --cut-step 0.05 \
    --out best_cuts.json --qa-dir qa_best
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
class CutRanges:
    dyL: float
    dyH: float
    W2L: float
    W2H: float
    eL: float
    tL: float
    tH: float

    def as_dict(self) -> Dict[str, float]:
        return dict(dyL=self.dyL, dyH=self.dyH, W2L=self.W2L, W2H=self.W2H, eL=self.eL, tL=self.tL, tH=self.tH)


@dataclass
class ObjectiveConfig:
    dx_window: Tuple[float, float]
    min_data_events: int
    min_template_events: int
    dy_min_width: float
    W2_min_width: float
    t_min_width: float
    penalty_tight: float


# -----------------------------
# Progress bar
# -----------------------------
def print_progress(i: int, n: int, t0: float, best_fom: float, prefix: str = "", bar_width: int = 30, every: int = 50):
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
        f"best {best_fom:.6g}"
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


def preselect_cpp_like_data(vz, ePS, eSH, trP, eHCAL, helicity) -> np.ndarray:
    reject = (np.abs(vz) > 0.27) & (ePS < 0.2) & (np.abs((eSH + ePS) / trP - 1.0) > 0.2) & (eHCAL < 0.025) & (np.abs(helicity) != 1)
    return ~reject


def preselect_cpp_like_sim(vz, ePS, eSH, trP, eHCAL) -> np.ndarray:
    reject = (np.abs(vz) > 0.27) & (ePS < 0.2) & (np.abs((eSH + ePS) / trP - 1.0) > 0.2) & (eHCAL < 0.025)
    return ~reject


def build_cut_mask(dy, W2, eHCAL, ct, cuts: CutRanges, apply_ct_cut: bool) -> np.ndarray:
    mask = (
        (dy > cuts.dyL) & (dy < cuts.dyH) &
        (W2 > cuts.W2L) & (W2 < cuts.W2H) &
        (eHCAL > cuts.eL)
    )
    if apply_ct_cut:
        mask &= (ct > cuts.tL) & (ct < cuts.tH)
    return mask


def snap(x: float, step: float) -> float:
    return round(x / step) * step


def snap_cuts(c: CutRanges, step: float) -> CutRanges:
    return CutRanges(
        dyL=snap(c.dyL, step),
        dyH=snap(c.dyH, step),
        W2L=snap(c.W2L, step),
        W2H=snap(c.W2H, step),
        eL=snap(c.eL, step),
        tL=snap(c.tL, step),
        tH=snap(c.tH, step),
    )


def compute_fom_and_details(
    bins: np.ndarray,
    dxD, dyD, W2D, ctD, eD, wD,
    dxN, dyN, W2N, ctN, eN, wN,
    dxP, dyP, W2P, ctP, eP, wP,
    bkg: Optional[Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]],
    cuts: CutRanges,
    cfg: ObjectiveConfig,
    treat_proton_as_background_if_no_bkg: bool = True,
):
    """
    Fit model in the requested parameterization (dx-only):

        f(dx) = N [ p(dx) + R n(dx) + Nbg b(dx) ]

    where p,n,b are *unit-area* dx PDFs (from templates after cuts).

    We fit DATA histogram hD(dx) with NNLS on PDF bases:
        hD ≈ cP*p + cN*n + cB*b  (c>=0)
    and map:
        N   = cP
        R   = cN / cP
        Nbg = cB / cP

    FOM is kept as before:
      - Signal = neutron component in dx window
      - Background = bkg component in dx window (or proton if no bkg)
    """
    if not (cuts.dyL < cuts.dyH and cuts.W2L < cuts.W2H and cuts.tL < cuts.tH):
        return -np.inf, None
    if (cuts.dyH - cuts.dyL) < cfg.dy_min_width:
        return -np.inf, None
    if (cuts.W2H - cuts.W2L) < cfg.W2_min_width:
        return -np.inf, None
    if (cuts.tH - cuts.tL) < cfg.t_min_width:
        return -np.inf, None

    # DATA uses coin_time cut; templates DO NOT
    mD = build_cut_mask(dyD, W2D, eD, ctD, cuts, apply_ct_cut=True)
    mN = build_cut_mask(dyN, W2N, eN, ctN, cuts, apply_ct_cut=False)
    mP = build_cut_mask(dyP, W2P, eP, ctP, cuts, apply_ct_cut=False)

    nD = int(mD.sum())
    nN = int(mN.sum())
    nP = int(mP.sum())

    if nD < cfg.min_data_events:
        return -np.inf, None
    if nN < cfg.min_template_events or nP < cfg.min_template_events:
        return -np.inf, None

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
    hB = None
    b_pdf = None
    nB = 0
    if has_bkg:
        dxB, dyB, W2B, ctB, eB, wB = bkg
        mB = build_cut_mask(dyB, W2B, eB, ctB, cuts, apply_ct_cut=False)
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
    if np.all(T.sum(axis=0) <= 0):
        return -np.inf, None

    coeff = nnls_fit(hD, T)
    cP = float(coeff[0])  # N
    cN = float(coeff[1])  # N*R
    cB = float(coeff[2]) if has_bkg else 0.0  # N*Nbg

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

    Nn = float(modelN[win].sum())
    if has_bkg:
        Nbkg = float(modelB[win].sum())
    else:
        Nbkg = float(modelP[win].sum()) if treat_proton_as_background_if_no_bkg else 0.0

    if (Nn + Nbkg) <= 0:
        return -np.inf, None

    fom = Nn / math.sqrt(Nn + Nbkg)

    if cfg.penalty_tight > 0:
        inv_width = (
            1.0 / max(cuts.dyH - cuts.dyL, 1e-9) +
            1.0 / max(cuts.W2H - cuts.W2L, 1e-9) +
            1.0 / max(cuts.tH - cuts.tL, 1e-9)
        )
        fom -= cfg.penalty_tight * inv_width

    details = dict(
        fom=float(fom),
        cuts=cuts.as_dict(),
        params=dict(N=N, R=R, Nbg=Nbg),
        coeff=dict(cP=cP, cN=cN, cB=cB),
        counts=dict(nData=nD, nN=nN, nP=nP, nB=nB),
        yields=dict(Nn=Nn, Nbkg=Nbkg),
        hists=dict(
            hD=hD,
            p_pdf=p_pdf, n_pdf=n_pdf, b_pdf=b_pdf,
            modelP=modelP, modelN=modelN, modelB=modelB, hTot=modelTot
        ),
        bins=bins,
        centers=centers,
    )
    return float(fom), details


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
    plt.title(f"dx fit @ best cuts | FOM={details['fom']:.4g}")
    plt.legend()
    save_fig(os.path.join(qa_dir, "qa_dx_fit_best.png"),
             os.path.join(qa_dir, "qa_dx_fit_best.pdf"))


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


def qa_data_vars_before_after(dyD, W2D, ctD, eHD, mask_best, qa_dir: str):
    ensure_dir(qa_dir)
    vars_ = [
        ("dy", dyD),
        ("W2", W2D),
        ("eHCAL", eHD),
        ("coin_time", ctD),
    ]

    plt.figure(figsize=(10, 8))
    for i, (name, arr) in enumerate(vars_, 1):
        plt.subplot(2, 2, i)
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


# -----------------------------
# CLI
# -----------------------------
def parse_args():
    ap = argparse.ArgumentParser(description="Cut optimization with chunked I/O, data-only coin_time, 0.05 grid cuts, and QA plots.")

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
    ap.add_argument("--br-IHWP", default="IHWP")

    ap.add_argument("--br-weight", default="weight")
    ap.add_argument("--br-fnucl", default="fnucl")
    ap.add_argument("--fnucl-n", type=float, default=0.0)
    ap.add_argument("--fnucl-p", type=float, default=1.0)

    ap.add_argument("--dx-min", type=float, default=-4.0)
    ap.add_argument("--dx-max", type=float, default=3.0)
    ap.add_argument("--dx-bins", type=int, default=100)
    ap.add_argument("--dx-win-low", type=float, default=-0.4)
    ap.add_argument("--dx-win-high", type=float, default=0.4)

    ap.add_argument("--method", choices=["random", "optuna"], default="random")
    ap.add_argument("--n-trials", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=1)

    ap.add_argument("--dyL-min", type=float, default=-0.6)
    ap.add_argument("--dyL-max", type=float, default=-0.3)
    ap.add_argument("--dyH-min", type=float, default=0.3)
    ap.add_argument("--dyH-max", type=float, default=0.6)

    ap.add_argument("--W2L-min", type=float, default=-2.0)
    ap.add_argument("--W2L-max", type=float, default=0.0)
    ap.add_argument("--W2H-min", type=float, default=1.2)
    ap.add_argument("--W2H-max", type=float, default=2.0)

    ap.add_argument("--eL-min", type=float, default=0.025)
    ap.add_argument("--eL-max", type=float, default=0.5)

    # DATA only
    ap.add_argument("--tL-min", type=float, default=170.0)
    ap.add_argument("--tL-max", type=float, default=190.0)
    ap.add_argument("--tH-min", type=float, default=170.0)
    ap.add_argument("--tH-max", type=float, default=190.0)

    ap.add_argument("--min-data-events", type=int, default=200)
    ap.add_argument("--min-template-events", type=int, default=2000)
    ap.add_argument("--dy-min-width", type=float, default=0.05)
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

    ap.add_argument("--cut-step", type=float, default=0.05, help="Grid step for final cuts (snap + local refinement)")
    ap.add_argument("--no-local-refine", action="store_true", help="Disable local grid refinement; only snap to grid")

    ap.add_argument("--out", default="best_cuts.json")
    ap.add_argument("--qa-dir", default="qa_best")
    return ap.parse_args()


def main():
    args = parse_args()
    rng = np.random.default_rng(args.seed)
    dtype = np.float32 if args.dtype == "float32" else np.float64

    bins = np.linspace(args.dx_min, args.dx_max, args.dx_bins + 1)

    # -----------------------------
    # Load DATA
    # -----------------------------
    data_branches = (args.br_dx, args.br_dy, args.br_W2, args.br_ct, args.br_eHCAL,
                     args.br_vz, args.br_ePS, args.br_eSH, args.br_trP, args.br_helicity)

    def keep_data(chunk):
        return preselect_cpp_like_data(chunk[args.br_vz], chunk[args.br_ePS], chunk[args.br_eSH],
                                       chunk[args.br_trP], chunk[args.br_eHCAL], chunk[args.br_helicity])

    max_data = args.max_data if args.max_data > 0 else None
    D = iterate_filtered(args.data, args.tree, data_branches, keep_data, args.step_size, rng, args.downsample_data, max_data, dtype)

    dxD = D[args.br_dx]; dyD = D[args.br_dy]; W2D = D[args.br_W2]; ctD = D[args.br_ct]; eHD = D[args.br_eHCAL]
    wD = np.ones_like(dxD, dtype=np.float64)

    # -----------------------------
    # Load SIM (no coin_time needed)
    # -----------------------------
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

    # -----------------------------
    # Load BKG (no coin_time needed)
    # -----------------------------
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
        dy_min_width=args.dy_min_width,
        W2_min_width=args.W2_min_width,
        t_min_width=args.t_min_width,
        penalty_tight=args.penalty_tight,
    )

    def eval_one(cuts: CutRanges):
        return compute_fom_and_details(
            bins,
            dxD, dyD, W2D, ctD, eHD, wD,
            dxN, dyN, W2N, ctN, eHN, wN,
            dxP, dyP, W2P, ctP, eHP, wP,
            bkg_tuple,
            cuts,
            cfg,
            treat_proton_as_background_if_no_bkg=True,
        )

    best_fom = -np.inf
    best_cuts: Optional[CutRanges] = None
    best_details = None

    if args.method == "random":
        t0 = time.time()
        for i in range(args.n_trials):
            cuts = CutRanges(
                dyL=float(rng.uniform(args.dyL_min, args.dyL_max)),
                dyH=float(rng.uniform(args.dyH_min, args.dyH_max)),
                W2L=float(rng.uniform(args.W2L_min, args.W2L_max)),
                W2H=float(rng.uniform(args.W2H_min, args.W2H_max)),
                eL=float(rng.uniform(args.eL_min, args.eL_max)),
                tL=float(rng.uniform(args.tL_min, args.tL_max)),
                tH=float(rng.uniform(args.tH_min, args.tH_max)),
            )
            fom, det = eval_one(cuts)
            if fom > best_fom:
                best_fom = fom
                best_cuts = cuts
                best_details = det
                print(f"\n[trial {i+1}/{args.n_trials}] NEW BEST FOM={best_fom:.6f} cuts={best_cuts.as_dict()}")
            print_progress(i, args.n_trials, t0, best_fom, prefix="Optimize ", every=50)

    else:
        if optuna is None:
            raise RuntimeError("optuna not installed. pip install optuna")
        if scipy_nnls is None:
            print("Warning: scipy not found; NNLS uses fallback. Recommended: pip install scipy")

        def objective(trial):
            cuts = CutRanges(
                dyL=trial.suggest_float("dyL", args.dyL_min, args.dyL_max),
                dyH=trial.suggest_float("dyH", args.dyH_min, args.dyH_max),
                W2L=trial.suggest_float("W2L", args.W2L_min, args.W2L_max),
                W2H=trial.suggest_float("W2H", args.W2H_min, args.W2H_max),
                eL=trial.suggest_float("eL", args.eL_min, args.eL_max),
                tL=trial.suggest_float("tL", args.tL_min, args.tL_max),
                tH=trial.suggest_float("tH", args.tH_min, args.tH_max),
            )
            fom, _ = eval_one(cuts)
            return fom

        study = optuna.create_study(direction="maximize", sampler=optuna.samplers.TPESampler(seed=args.seed))
        study.optimize(objective, n_trials=args.n_trials)
        best_cuts = CutRanges(**study.best_params)
        best_fom, best_details = eval_one(best_cuts)

    # -----------------------------
    # Snap to grid + local refine
    # -----------------------------
    final_cuts = best_cuts
    final_fom = best_fom
    final_details = best_details

    if best_cuts is not None and np.isfinite(best_fom):
        snapped = snap_cuts(best_cuts, args.cut_step)
        snap_fom, snap_det = eval_one(snapped)
        final_cuts, final_fom, final_details = snapped, snap_fom, snap_det

        if not args.no_local_refine:
            step = args.cut_step
            base = snapped
            offsets = [-step, 0.0, step]
            best_grid_fom = -np.inf
            best_grid_cuts = None
            best_grid_det = None

            import itertools
            for d_dyL, d_dyH, d_W2L, d_W2H, d_eL, d_tL, d_tH in itertools.product(offsets, repeat=7):
                c = CutRanges(
                    dyL=base.dyL + d_dyL,
                    dyH=base.dyH + d_dyH,
                    W2L=base.W2L + d_W2L,
                    W2H=base.W2H + d_W2H,
                    eL=base.eL + d_eL,
                    tL=base.tL + d_tL,
                    tH=base.tH + d_tH,
                )

                # keep within original search bounds
                if not (args.dyL_min <= c.dyL <= args.dyL_max and args.dyH_min <= c.dyH <= args.dyH_max):
                    continue
                if not (args.W2L_min <= c.W2L <= args.W2L_max and args.W2H_min <= c.W2H <= args.W2H_max):
                    continue
                if not (args.eL_min <= c.eL <= args.eL_max):
                    continue
                if not (args.tL_min <= c.tL <= args.tL_max and args.tH_min <= c.tH <= args.tH_max):
                    continue

                c = snap_cuts(c, step)
                fom, det = eval_one(c)
                if fom > best_grid_fom:
                    best_grid_fom = fom
                    best_grid_cuts = c
                    best_grid_det = det

            if best_grid_cuts is not None and np.isfinite(best_grid_fom):
                final_cuts, final_fom, final_details = best_grid_cuts, best_grid_fom, best_grid_det

    # -----------------------------
    # QA
    # -----------------------------
    if final_cuts is not None and np.isfinite(final_fom) and final_details is not None:
        ensure_dir(args.qa_dir)
        mD_best = build_cut_mask(dyD, W2D, eHD, ctD, final_cuts, apply_ct_cut=True)

        qa_dx_fit(final_details, args.qa_dir, (args.dx_win_low, args.dx_win_high))
        qa_dx_residual(final_details, args.qa_dir)
        qa_data_vars_before_after(dyD, W2D, ctD, eHD, mD_best, args.qa_dir)

        with open(os.path.join(args.qa_dir, "qa_best_fit_summary.json"), "w") as f:
            json.dump(
                {
                    "best_fom": float(final_fom),
                    "best_cuts": final_cuts.as_dict(),
                    "params": final_details.get("params", {}),
                    "coeff_pdf_fit": final_details.get("coeff", {}),
                    "counts_after_cuts": final_details["counts"],
                    "yields_in_dx_window": final_details["yields"],
                    "dx_window": [args.dx_win_low, args.dx_win_high],
                    "coin_time_cut_applied_to": "data_only",
                    "cut_step": args.cut_step,
                    "local_refine": (not args.no_local_refine),
                },
                f,
                indent=2,
                sort_keys=True,
            )

    out = {
        "best_fom": float(final_fom) if np.isfinite(final_fom) else float("-inf"),
        "best_cuts": final_cuts.as_dict() if final_cuts else None,
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
            "max_data": args.max_data, "max_sim": args.max_sim, "max_bkg": args.max_bkg,
            "downsample_data": args.downsample_data, "downsample_sim": args.downsample_sim, "downsample_bkg": args.downsample_bkg,
            "qa_dir": args.qa_dir,
            "coin_time_cut_applied_to": "data_only",
            "cut_step": args.cut_step,
            "local_refine": (not args.no_local_refine),
        }
    }
    with open(args.out, "w") as f:
        json.dump(out, f, indent=2, sort_keys=True)

    print("\n=== DONE ===")
    print("Final Best FOM:", out["best_fom"])
    print("Final Best cuts:", out["best_cuts"])
    print("Wrote:", os.path.abspath(args.out))
    if final_cuts is not None and np.isfinite(final_fom):
        print("QA plots saved in:", os.path.abspath(args.qa_dir))


if __name__ == "__main__":
    main()
