#!/usr/bin/env python3
"""
optimize_cuts_from_cpp_stream.py

Memory-safe (OOM-proof) version of optimize_cuts_from_cpp.py:
- Reads ROOT TTrees in CHUNKS (uproot.iterate) instead of loading everything at once.
- Applies the SAME C++-like preselection you used in OptimizeCuts.cpp.
- Then runs random/optuna optimization with an in-terminal progress bar.

Default tree/branch names match your OptimizeCuts.cpp:
  Tree: Tout
  Branches: vz, ePS, eSH, trP, dx, dy, W2, coin_time, eHCAL, helicity, IHWP
  Sim label: fnucl (0=n, 1=p), weight: weight
  coin_time in sim/bkg is optional; if missing, it's treated as 0.

Typical use (random):
  python3 -u optimize_cuts_from_cpp_stream.py --data data.root --sim sim_np.root --bkg sim_inel.root \
      --method random --n-trials 10000 --out best.json

If you still want to reduce memory/time, cap events (highly recommended for large samples):
  --max-data 300000 --max-sim 400000 --max-bkg 300000

Chunk size:
  --step-size "200 MB"   (default)

Notes:
- The optimization loop still needs arrays in memory, but now you only keep *selected* rows.
- If even the selected rows are huge, cap them with --max-* or downsample (see --downsample-*).
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
    """Simple stdout progress bar with ETA."""
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
    """data_h (nbins,), templates (nbins, ntemp)"""
    if scipy_nnls is not None:
        coeff, _ = scipy_nnls(templates, data_h)
        return coeff
    coeff = np.linalg.lstsq(templates, data_h, rcond=None)[0]
    return np.maximum(coeff, 0.0)


def preselect_cpp_like_data(vz, ePS, eSH, trP, eHCAL, helicity) -> np.ndarray:
    # Mirrors OptimizeCuts.cpp (reject only if ALL are true):
    reject = (np.abs(vz) > 0.27) & (ePS < 0.2) & (np.abs((eSH + ePS) / trP - 1.0) > 0.2) & (eHCAL < 0.025) & (np.abs(helicity) != 1)
    return ~reject


def preselect_cpp_like_sim(vz, ePS, eSH, trP, eHCAL) -> np.ndarray:
    reject = (np.abs(vz) > 0.27) & (ePS < 0.2) & (np.abs((eSH + ePS) / trP - 1.0) > 0.2) & (eHCAL < 0.025)
    return ~reject


def build_cut_mask(dy, W2, eHCAL, ct, cuts: CutRanges, apply_ct_cut: bool = True) -> np.ndarray:
    mask = (
        (dy > cuts.dyL) & (dy < cuts.dyH) &
        (W2 > cuts.W2L) & (W2 < cuts.W2H) &
        (eHCAL > cuts.eL)
    )
    if apply_ct_cut:
        mask &= (ct > cuts.tL) & (ct < cuts.tH)
    return mask


def compute_fom(
    bins: np.ndarray,
    # data
    dxD, dyD, W2D, ctD, eD, wD,
    # neutron
    dxN, dyN, W2N, ctN, eN, wN,
    # proton
    dxP, dyP, W2P, ctP, eP, wP,
    # bkg (optional)
    bkg: Optional[Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]],
    cuts: CutRanges,
    cfg: ObjectiveConfig,
    treat_proton_as_background_if_no_bkg: bool = True,
) -> float:
    if not (cuts.dyL < cuts.dyH and cuts.W2L < cuts.W2H and cuts.tL < cuts.tH):
        return -np.inf
    if (cuts.dyH - cuts.dyL) < cfg.dy_min_width:
        return -np.inf
    if (cuts.W2H - cuts.W2L) < cfg.W2_min_width:
        return -np.inf
    if (cuts.tH - cuts.tL) < cfg.t_min_width:
        return -np.inf

    mD = build_cut_mask(dyD, W2D, eD, ctD, cuts, apply_ct_cut=True)
    mN = build_cut_mask(dyN, W2N, eN, ctN, cuts, apply_ct_cut=False)
    mP = build_cut_mask(dyP, W2P, eP, ctP, cuts, apply_ct_cut=False)

    if mD.sum() < cfg.min_data_events:
        return -np.inf
    if mN.sum() < cfg.min_template_events or mP.sum() < cfg.min_template_events:
        return -np.inf

    hD = hist1d(dxD[mD], wD[mD], bins)
    hN = hist1d(dxN[mN], wN[mN], bins)
    hP = hist1d(dxP[mP], wP[mP], bins)

    templates = [hN, hP]
    has_bkg = bkg is not None
    if has_bkg:
        dxB, dyB, W2B, ctB, eB, wB = bkg
        mB = build_cut_mask(dyB, W2B, eB, ctB, cuts, apply_ct_cut=False)
        if mB.sum() < cfg.min_template_events:
            return -np.inf
        hB = hist1d(dxB[mB], wB[mB], bins)
        templates.append(hB)

    T = np.vstack(templates).T
    if np.all(T.sum(axis=0) <= 0):
        return -np.inf

    coeff = nnls_fit(hD, T)
    aN = float(coeff[0])
    aP = float(coeff[1])
    aB = float(coeff[2]) if has_bkg else 0.0

    modelN = aN * hN
    modelP = aP * hP
    modelB = (aB * templates[2]) if has_bkg else (0.0 * hN)

    lo, hi = cfg.dx_window
    centers = 0.5 * (bins[:-1] + bins[1:])
    win = (centers >= lo) & (centers <= hi)

    Nn = modelN[win].sum()
    if has_bkg:
        Nbkg = modelB[win].sum()
    else:
        Nbkg = modelP[win].sum() if treat_proton_as_background_if_no_bkg else 0.0

    if (Nn + Nbkg) <= 0:
        return -np.inf

    fom = Nn / math.sqrt(Nn + Nbkg)

    if cfg.penalty_tight > 0:
        inv_width = (
            1.0 / max(cuts.dyH - cuts.dyL, 1e-9) +
            1.0 / max(cuts.W2H - cuts.W2L, 1e-9) +
            1.0 / max(cuts.tH - cuts.tL, 1e-9)
        )
        fom -= cfg.penalty_tight * inv_width

    return float(fom)


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
    """
    Stream over ROOT file, apply keep_fn(chunk)->bool mask, optionally downsample, and cap total kept rows.
    Returns concatenated arrays for each requested branch.
    """
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
            # take only what's needed to reach max_keep
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

    # Concatenate
    out2 = {}
    for br, pieces in out.items():
        out2[br] = np.concatenate(pieces) if len(pieces) else np.array([], dtype=dtype)
    return out2


# -----------------------------
# CLI
# -----------------------------
def parse_args():
    ap = argparse.ArgumentParser(description="Cut optimization with chunked I/O (prevents OOM).")

    ap.add_argument("--data", required=True)
    ap.add_argument("--sim", required=True)
    ap.add_argument("--bkg", default="")
    ap.add_argument("--tree", default="Tout")

    # Branch defaults from your OptimizeCuts.cpp
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

    # dx hist
    ap.add_argument("--dx-min", type=float, default=-4.0)
    ap.add_argument("--dx-max", type=float, default=3.0)
    ap.add_argument("--dx-bins", type=int, default=100)
    ap.add_argument("--dx-win-low", type=float, default=-0.7)
    ap.add_argument("--dx-win-high", type=float, default=0.7)

    # Optimization
    ap.add_argument("--method", choices=["random", "optuna"], default="random")
    ap.add_argument("--n-trials", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=1)

    # Search bounds
    ap.add_argument("--dyL-min", type=float, default=-0.5)
    ap.add_argument("--dyL-max", type=float, default=0.4)
    ap.add_argument("--dyH-min", type=float, default=-0.4)
    ap.add_argument("--dyH-max", type=float, default=0.5)

    ap.add_argument("--W2L-min", type=float, default=-2.0)
    ap.add_argument("--W2L-max", type=float, default=1.0)
    ap.add_argument("--W2H-min", type=float, default=0.8)
    ap.add_argument("--W2H-max", type=float, default=2.0)

    ap.add_argument("--eL-min", type=float, default=0.025)
    ap.add_argument("--eL-max", type=float, default=0.5)

    ap.add_argument("--tL-min", type=float, default=170.0)
    ap.add_argument("--tL-max", type=float, default=190.0)
    ap.add_argument("--tH-min", type=float, default=170.0)
    ap.add_argument("--tH-max", type=float, default=190.0)

    # Validity constraints
    ap.add_argument("--min-data-events", type=int, default=200)
    ap.add_argument("--min-template-events", type=int, default=200)
    ap.add_argument("--dy-min-width", type=float, default=0.05)
    ap.add_argument("--W2-min-width", type=float, default=0.30)
    ap.add_argument("--t-min-width", type=float, default=3.0)
    ap.add_argument("--penalty-tight", type=float, default=0.0)

    # Chunking / memory controls
    ap.add_argument("--step-size", default="200 MB", help='uproot.iterate step_size, e.g. "100 MB"')
    ap.add_argument("--dtype", choices=["float32", "float64"], default="float32", help="Store arrays in this dtype to save RAM")

    ap.add_argument("--max-data", type=int, default=0, help="Max kept DATA rows after preselection (0=all)")
    ap.add_argument("--max-sim", type=int, default=0, help="Max kept SIM rows after preselection (0=all)")
    ap.add_argument("--max-bkg", type=int, default=0, help="Max kept BKG rows after preselection (0=all)")

    ap.add_argument("--downsample-data", type=float, default=1.0, help="Keep fraction of DATA rows (after preselection)")
    ap.add_argument("--downsample-sim", type=float, default=1.0, help="Keep fraction of SIM rows (after preselection)")
    ap.add_argument("--downsample-bkg", type=float, default=1.0, help="Keep fraction of BKG rows (after preselection)")

    ap.add_argument("--out", default="best_cuts.json")
    return ap.parse_args()


def main():
    args = parse_args()
    rng = np.random.default_rng(args.seed)
    dtype = np.float32 if args.dtype == "float32" else np.float64

    bins = np.linspace(args.dx_min, args.dx_max, args.dx_bins + 1)

    # -----------------------------
    # Load DATA (chunked)
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
    # Load SIM (chunked)
    # -----------------------------
    sim_has_ct = False
    with uproot.open(args.sim) as f:
        sim_has_ct = args.br_ct in f[args.tree].keys()

    sim_branches = (args.br_dx, args.br_dy, args.br_W2, args.br_eHCAL,
                    args.br_vz, args.br_ePS, args.br_eSH, args.br_trP,
                    args.br_weight, args.br_fnucl) + ((args.br_ct,) if sim_has_ct else ())

    def keep_sim(chunk):
        return preselect_cpp_like_sim(chunk[args.br_vz], chunk[args.br_ePS], chunk[args.br_eSH],
                                      chunk[args.br_trP], chunk[args.br_eHCAL])

    max_sim = args.max_sim if args.max_sim > 0 else None
    S = iterate_filtered(args.sim, args.tree, sim_branches, keep_sim, args.step_size, rng, args.downsample_sim, max_sim, dtype)

    dxS = S[args.br_dx]; dyS = S[args.br_dy]; W2S = S[args.br_W2]
    eHS = S[args.br_eHCAL]; wS = S[args.br_weight]; fnucl = S[args.br_fnucl]
    ctS = S[args.br_ct] if sim_has_ct else np.zeros_like(dxS, dtype=dtype)

    ok_np = (fnucl == args.fnucl_n) | (fnucl == args.fnucl_p)
    dxS, dyS, W2S, ctS, eHS, wS, fnucl = dxS[ok_np], dyS[ok_np], W2S[ok_np], ctS[ok_np], eHS[ok_np], wS[ok_np], fnucl[ok_np]

    isN = (fnucl == args.fnucl_n)
    isP = (fnucl == args.fnucl_p)
    if isN.sum() == 0 or isP.sum() == 0:
        raise RuntimeError(f"Sim split produced empty sample(s): n={isN.sum()} p={isP.sum()} (check fnucl values).")

    dxN, dyN, W2N, ctN, eHN, wN = dxS[isN], dyS[isN], W2S[isN], ctS[isN], eHS[isN], wS[isN]
    dxP, dyP, W2P, ctP, eHP, wP = dxS[isP], dyS[isP], W2S[isP], ctS[isP], eHS[isP], wS[isP]

    # -----------------------------
    # Load BKG (chunked, optional)
    # -----------------------------
    bkg_tuple = None
    if args.bkg:
        bkg_has_ct = False
        with uproot.open(args.bkg) as f:
            bkg_has_ct = args.br_ct in f[args.tree].keys()

        bkg_branches = (args.br_dx, args.br_dy, args.br_W2, args.br_eHCAL,
                        args.br_vz, args.br_ePS, args.br_eSH, args.br_trP,
                        args.br_weight) + ((args.br_ct,) if bkg_has_ct else ())

        def keep_bkg(chunk):
            return preselect_cpp_like_sim(chunk[args.br_vz], chunk[args.br_ePS], chunk[args.br_eSH],
                                          chunk[args.br_trP], chunk[args.br_eHCAL])

        max_bkg = args.max_bkg if args.max_bkg > 0 else None
        B = iterate_filtered(args.bkg, args.tree, bkg_branches, keep_bkg, args.step_size, rng, args.downsample_bkg, max_bkg, dtype)

        dxB = B[args.br_dx]; dyB = B[args.br_dy]; W2B = B[args.br_W2]; eHB = B[args.br_eHCAL]
        wB = B[args.br_weight]
        ctB = B[args.br_ct] if bkg_has_ct else np.zeros_like(dxB, dtype=dtype)
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

    def eval_one(cuts: CutRanges) -> float:
        return compute_fom(
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
            fom = eval_one(cuts)
            if fom > best_fom:
                best_fom = fom
                best_cuts = cuts
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
            return eval_one(cuts)

        study = optuna.create_study(direction="maximize", sampler=optuna.samplers.TPESampler(seed=args.seed))
        study.optimize(objective, n_trials=args.n_trials)
        best_fom = float(study.best_value)
        best_cuts = CutRanges(**study.best_params)

    out = {
        "best_fom": best_fom,
        "best_cuts": best_cuts.as_dict() if best_cuts else None,
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
        }
    }
    with open(args.out, "w") as f:
        json.dump(out, f, indent=2, sort_keys=True)

    print("\n=== DONE ===")
    print("Best FOM:", best_fom)
    print("Best cuts:", best_cuts.as_dict() if best_cuts else None)
    print("Wrote:", os.path.abspath(args.out))


if __name__ == "__main__":
    main()
