#!/usr/bin/env python3
"""
optimize_cuts_from_cpp.py

Python cut-optimizer matching OptimizeCuts.cpp branch/tree defaults.

From OptimizeCuts.cpp:
  Tree: "Tout"
  Data branches:
    vz, ePS, eSH, trP, dx, dy, W2, coin_time, eHCAL, helicity, IHWP
  Sim branches (n/p):
    vz, ePS, eSH, trP, dx, dy, W2, coin_time (optional), eHCAL, weight, fnucl
  Bkg sim branches:
    vz, ePS, eSH, trP, dx, dy, W2, coin_time (optional), eHCAL, weight

Method:
  For each trial set of cut edges:
    - Apply cuts to Data and templates (n, p, bkg)
    - Histogram dx for data/templates
    - Fit data dx with non-negative least squares:
          Data ≈ aN*N + aP*P + aB*B
    - FOM in dx window:  Nn / sqrt(Nn + Nbkg)
        where Nbkg defaults to fitted background template (aB*B).
        (You can switch to treat protons as "background" if desired.)

Optimization:
  - Random search (default)
  - Optuna (Bayesian/TPE) if installed

Requirements:
  pip install uproot numpy
Optional (recommended):
  pip install scipy optuna

Examples:
  python optimize_cuts_from_cpp.py \
    --data data.root --sim sim_np.root --bkg sim_bkg.root \
    --n-trials 5000 --method random --out best_cuts.json

  python optimize_cuts_from_cpp.py \
    --data data.root --sim sim_np.root --bkg sim_bkg.root \
    --n-trials 2000 --method optuna --out best_cuts.json
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
# ROOT loading helpers
# -----------------------------
def load_tree_arrays(rootfile: str, treename: str, branches: Tuple[str, ...]) -> Dict[str, np.ndarray]:
    with uproot.open(rootfile) as f:
        t = f[treename]
        arr = t.arrays(list(branches), library="np")
    return arr


def has_branch(rootfile: str, treename: str, branch: str) -> bool:
    with uproot.open(rootfile) as f:
        t = f[treename]
        return branch in t.keys()


def safe_get(arr: Dict[str, np.ndarray], key: str) -> np.ndarray:
    if key not in arr:
        raise KeyError(f"Branch '{key}' not found. Available (first ~30): {list(arr.keys())[:30]}")
    return arr[key]


# -----------------------------
# Numerics
# -----------------------------
def hist1d(x: np.ndarray, w: np.ndarray, bins: np.ndarray) -> np.ndarray:
    h, _ = np.histogram(x, bins=bins, weights=w)
    return h.astype(np.float64)


def nnls_fit(data_h: np.ndarray, templates: np.ndarray) -> np.ndarray:
    """
    data_h: (nbins,)
    templates: (nbins, ntemp)
    returns coeff: (ntemp,)
    """
    if scipy_nnls is not None:
        coeff, _ = scipy_nnls(templates, data_h)
        return coeff

    # Fallback: clipped least squares (not as good as NNLS)
    coeff = np.linalg.lstsq(templates, data_h, rcond=None)[0]
    return np.maximum(coeff, 0.0)


def preselect_cpp_like(vz, ePS, eSH, trP, eHCAL, helicity: Optional[np.ndarray], is_data: bool) -> np.ndarray:
    """
    Mirrors your OptimizeCuts.cpp logic:

    Data reject condition:
      if (abs(vz)>0.27 && ePS<0.2 &&
          abs((eSH+ePS)/trP - 1.0) > 0.2 && eHCAL<0.025 &&
          abs(helicity)!=1) continue;

    Sim/bkg reject condition:
      if (abs(vz)>0.27 && ePS<0.2 &&
          abs((eSH+ePS)/trP - 1.0) > 0.2 && eHCAL<0.025) continue;

    NOTE: This rejects only when ALL those are simultaneously true.
    """
    cond_core = (np.abs(vz) > 0.27) & (ePS < 0.2) & (np.abs((eSH + ePS) / trP - 1.0) > 0.2) & (eHCAL < 0.025)
    if is_data:
        if helicity is None:
            raise ValueError("Data preselection needs helicity branch.")
        cond = cond_core & (np.abs(helicity) != 1)
    else:
        cond = cond_core

    # keep mask is NOT rejected
    return ~cond


def build_cut_mask(dy, W2, eHCAL, ct, cuts: CutRanges) -> np.ndarray:
    return (
        (dy > cuts.dyL) & (dy < cuts.dyH) &
        (W2 > cuts.W2L) & (W2 < cuts.W2H) &
        (eHCAL > cuts.eL) &
        (ct > cuts.tL) & (ct < cuts.tH)
    )


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
    # Sanity constraints
    if not (cuts.dyL < cuts.dyH and cuts.W2L < cuts.W2H and cuts.tL < cuts.tH):
        return -np.inf
    if (cuts.dyH - cuts.dyL) < cfg.dy_min_width:
        return -np.inf
    if (cuts.W2H - cuts.W2L) < cfg.W2_min_width:
        return -np.inf
    if (cuts.tH - cuts.tL) < cfg.t_min_width:
        return -np.inf

    mD = build_cut_mask(dyD, W2D, eD, ctD, cuts)
    mN = build_cut_mask(dyN, W2N, eN, ctN, cuts)
    mP = build_cut_mask(dyP, W2P, eP, ctP, cuts)

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
        mB = build_cut_mask(dyB, W2B, eB, ctB, cuts)
        if mB.sum() < cfg.min_template_events:
            return -np.inf
        hB = hist1d(dxB[mB], wB[mB], bins)
        templates.append(hB)

    T = np.vstack(templates).T  # (nbins, ntemp)
    if np.all(T.sum(axis=0) <= 0):
        return -np.inf

    coeff = nnls_fit(hD, T)
    aN = float(coeff[0])
    aP = float(coeff[1])
    aB = float(coeff[2]) if has_bkg else 0.0

    modelN = aN * hN
    modelP = aP * hP
    modelB = (aB * templates[2]) if has_bkg else (0.0 * hN)

    # Integrate in dx window
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

    # Optional mild penalty for very tight cuts
    if cfg.penalty_tight > 0:
        inv_width = (
            1.0 / max(cuts.dyH - cuts.dyL, 1e-9) +
            1.0 / max(cuts.W2H - cuts.W2L, 1e-9) +
            1.0 / max(cuts.tH - cuts.tL, 1e-9)
        )
        fom -= cfg.penalty_tight * inv_width

    return float(fom)


# -----------------------------
# CLI
# -----------------------------

def print_progress(i: int, n: int, t0: float, best_fom: float, bar_width: int = 30, every: int = 50):
    """
    Simple stdout progress bar with ETA. Updates every `every` iterations.
    Call with i = current iteration index (0-based).
    """
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
        f"\r[{bar}] {100*frac:6.2f}% "
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


def parse_args():
    ap = argparse.ArgumentParser(description="Optimize dy/W2/eHCAL/coin_time cuts using dx template fits (matching OptimizeCuts.cpp defaults).")

    ap.add_argument("--data", required=True, help="Data ROOT file")
    ap.add_argument("--sim", required=True, help="Sim ROOT file containing n+p (fnucl==0/1)")
    ap.add_argument("--bkg", default="", help="Optional bkg sim ROOT file (template background)")

    ap.add_argument("--tree", default="Tout", help='Tree name (default: "Tout" from OptimizeCuts.cpp)')

    # Branch defaults from your C++ macro
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

    ap.add_argument("--br-weight", default="weight", help='Sim/bkg weight branch (default: "weight")')
    ap.add_argument("--br-fnucl", default="fnucl", help='Sim label branch (default: "fnucl")')
    ap.add_argument("--fnucl-n", type=float, default=0.0, help="Value for neutron in fnucl (default 0)")
    ap.add_argument("--fnucl-p", type=float, default=1.0, help="Value for proton in fnucl (default 1)")

    # Histogram settings (matches your macro defaults: nbinsDx=100, dx range [-4,3])
    ap.add_argument("--dx-min", type=float, default=-4.0)
    ap.add_argument("--dx-max", type=float, default=3.0)
    ap.add_argument("--dx-bins", type=int, default=100)

    # FOM dx window (set to what you use in C++ if different)
    ap.add_argument("--dx-win-low", type=float, default=-0.7)
    ap.add_argument("--dx-win-high", type=float, default=0.7)

    # Search bounds (tune as you like)
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

    # Optimization controls
    ap.add_argument("--method", choices=["random", "optuna"], default="random")
    ap.add_argument("--n-trials", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=1)

    # Validity constraints
    ap.add_argument("--min-data-events", type=int, default=200)
    ap.add_argument("--min-template-events", type=int, default=200)
    ap.add_argument("--dy-min-width", type=float, default=0.05)
    ap.add_argument("--W2-min-width", type=float, default=0.30)
    ap.add_argument("--t-min-width", type=float, default=3.0)
    ap.add_argument("--penalty-tight", type=float, default=0.0)

    ap.add_argument("--out", default="best_cuts.json")
    return ap.parse_args()


# -----------------------------
# Main
# -----------------------------
def main():
    args = parse_args()
    rng = np.random.default_rng(args.seed)

    bins = np.linspace(args.dx_min, args.dx_max, args.dx_bins + 1)

    # coin_time may be missing in sim/bkg (your C++ checks this)
    sim_has_ct = has_branch(args.sim, args.tree, args.br_ct)
    bkg_has_ct = has_branch(args.bkg, args.tree, args.br_ct) if args.bkg else False

    # Load DATA branches (exactly as in C++ data section)
    data_br = (
        args.br_vz, args.br_ePS, args.br_eSH, args.br_trP,
        args.br_dx, args.br_dy, args.br_W2, args.br_ct, args.br_eHCAL,
        args.br_helicity, args.br_IHWP
    )
    D = load_tree_arrays(args.data, args.tree, data_br)

    # Load SIM branches (as in C++ sim section)
    sim_br = (
        args.br_vz, args.br_ePS, args.br_eSH, args.br_trP,
        args.br_dx, args.br_dy, args.br_W2, args.br_eHCAL,
        args.br_weight, args.br_fnucl
    )
    if sim_has_ct:
        sim_br = sim_br + (args.br_ct,)
    S = load_tree_arrays(args.sim, args.tree, sim_br)

    # Load BKG branches (as in C++ bkg section)
    B = None
    if args.bkg:
        bkg_br = (
            args.br_vz, args.br_ePS, args.br_eSH, args.br_trP,
            args.br_dx, args.br_dy, args.br_W2, args.br_eHCAL,
            args.br_weight
        )
        if bkg_has_ct:
            bkg_br = bkg_br + (args.br_ct,)
        B = load_tree_arrays(args.bkg, args.tree, bkg_br)

    # -----------------------------
    # Extract + preselection
    # -----------------------------
    # DATA
    vzD = safe_get(D, args.br_vz)
    ePSD = safe_get(D, args.br_ePS)
    eSHD = safe_get(D, args.br_eSH)
    trPD = safe_get(D, args.br_trP)
    dxD = safe_get(D, args.br_dx)
    dyD = safe_get(D, args.br_dy)
    W2D = safe_get(D, args.br_W2)
    ctD = safe_get(D, args.br_ct)
    eHD = safe_get(D, args.br_eHCAL)
    helD = safe_get(D, args.br_helicity)

    keepD = preselect_cpp_like(vzD, ePSD, eSHD, trPD, eHD, helD, is_data=True)

    dxD, dyD, W2D, ctD, eHD = dxD[keepD], dyD[keepD], W2D[keepD], ctD[keepD], eHD[keepD]
    wD = np.ones_like(dxD, dtype=np.float64)  # C++ sets data weight = 1

    # SIM
    vzS = safe_get(S, args.br_vz)
    ePSS = safe_get(S, args.br_ePS)
    eSHS = safe_get(S, args.br_eSH)
    trPS = safe_get(S, args.br_trP)
    dxS = safe_get(S, args.br_dx)
    dyS = safe_get(S, args.br_dy)
    W2S = safe_get(S, args.br_W2)
    eHS = safe_get(S, args.br_eHCAL)
    wS = safe_get(S, args.br_weight)
    fnucl = safe_get(S, args.br_fnucl)
    ctS = safe_get(S, args.br_ct) if sim_has_ct else np.zeros_like(dxS, dtype=np.float64)

    keepS = preselect_cpp_like(vzS, ePSS, eSHS, trPS, eHS, helicity=None, is_data=False)
    dxS, dyS, W2S, ctS, eHS, wS, fnucl = dxS[keepS], dyS[keepS], W2S[keepS], ctS[keepS], eHS[keepS], wS[keepS], fnucl[keepS]

    # Accept only fnucl 0 or 1, like C++
    ok_np = (fnucl == args.fnucl_n) | (fnucl == args.fnucl_p)
    dxS, dyS, W2S, ctS, eHS, wS, fnucl = dxS[ok_np], dyS[ok_np], W2S[ok_np], ctS[ok_np], eHS[ok_np], wS[ok_np], fnucl[ok_np]

    isN = (fnucl == args.fnucl_n)
    isP = (fnucl == args.fnucl_p)
    if isN.sum() == 0 or isP.sum() == 0:
        raise RuntimeError(f"After preselection/split, sim is empty: n={isN.sum()} p={isP.sum()} (check fnucl values).")

    dxN, dyN, W2N, ctN, eHN, wN = dxS[isN], dyS[isN], W2S[isN], ctS[isN], eHS[isN], wS[isN]
    dxP, dyP, W2P, ctP, eHP, wP = dxS[isP], dyS[isP], W2S[isP], ctS[isP], eHS[isP], wS[isP]

    # BKG
    bkg_tuple = None
    if B is not None:
        vzB = safe_get(B, args.br_vz)
        ePSB = safe_get(B, args.br_ePS)
        eSHB = safe_get(B, args.br_eSH)
        trPB = safe_get(B, args.br_trP)
        dxB = safe_get(B, args.br_dx)
        dyB = safe_get(B, args.br_dy)
        W2B = safe_get(B, args.br_W2)
        eHB = safe_get(B, args.br_eHCAL)
        wB = safe_get(B, args.br_weight)
        ctB = safe_get(B, args.br_ct) if bkg_has_ct else np.zeros_like(dxB, dtype=np.float64)

        keepB = preselect_cpp_like(vzB, ePSB, eSHB, trPB, eHB, helicity=None, is_data=False)
        dxB, dyB, W2B, ctB, eHB, wB = dxB[keepB], dyB[keepB], W2B[keepB], ctB[keepB], eHB[keepB], wB[keepB]
        bkg_tuple = (dxB, dyB, W2B, ctB, eHB, wB)

    # Objective config
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

    # -----------------------------
    # Optimization
    # -----------------------------
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

            # Progress bar update (prints every N trials; last trial always prints)
            print_progress(i, args.n_trials, t0, best_fom, every=50)


    elif args.method == "optuna":
        if optuna is None:
            raise RuntimeError("optuna not installed. Install with: pip install optuna")
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

    else:
        raise RuntimeError("Unknown method")

    # -----------------------------
    # Output
    # -----------------------------
    out = {
        "best_fom": best_fom,
        "best_cuts": best_cuts.as_dict() if best_cuts else None,
        "meta": {
            "tree": args.tree,
            "data": args.data,
            "sim": args.sim,
            "bkg": args.bkg,
            "sim_has_coin_time": sim_has_ct,
            "bkg_has_coin_time": bkg_has_ct,
            "dx_hist": {"min": args.dx_min, "max": args.dx_max, "bins": args.dx_bins},
            "dx_window": [args.dx_win_low, args.dx_win_high],
            "method": args.method,
            "n_trials": args.n_trials,
            "seed": args.seed,
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

