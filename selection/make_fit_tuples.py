from pathlib import Path
import numpy as np
import pandas as pd
import ROOT

from bdt_run3 import (
    INPUT_FEATURES,
    AUX_FEATURES,
    rdf_to_frame,
    skim_dir,
    MODEL_STEM,
    N_SPLITS,
    OUT_DIR,
    SCORED_TUPLE,
)
import bdt_util as bu
from sb_fit import x_min2, x_max2, BDT_CUT, store_count, _predict_one, ensemble_score, NEW_CUTS
from mass_vetoes import run_veto_pipeline, mass_vetoes

ROOT.gROOT.SetBatch(True)

OS_REGION = "sideband"  # [sideband or full]

MC_3PI = skim_dir / "taumu_sim_OS.root"
MC_3PIPI0 = skim_dir / "taumu_pi0_sim_OS.root"
DAT_OS = skim_dir / "taumu_dat_OS.root"

GEN_3PI = skim_dir / "taumu_gen_counts.npz"
GEN_3PIPI0 = skim_dir / "taumu_pi_gen_counts.npz"

BR_3PI = 0.0931
BR_3PIPI0 = 0.0462
BR_TOT = BR_3PI + BR_3PIPI0

CACHE_FILE = OUT_DIR / "cache" / f"sideband_fit_os.npz"

OUT_MC = OUT_DIR / "cache" / "scored_mc_combined.root"
OUT_OS = OUT_DIR / "cache" / "scored_data_os.root"

FEAT = [n for n, _ in INPUT_FEATURES]

END_VZ_CUT_R = "tau_END_VZ - Lb_END_VZ > -1"
HLT1_CUT_R = "Lb_Hlt1TrackMVADecision_TOS || Lb_Hlt1TwoTrackMVADecision_TOS"
PID_CUT_R = "(proton_PROBNN_P > 0.4) && (proton_PROBNN_K < 0.6) && (proton_PROBNN_PI < 0.4) && (kaon_PROBNN_PI < 0.5) && (kaon_PROBNN_K > 0.4) && (mu_PROBNN_MU > 0.5)"

END_VZ_CUT = "tau_ENDVERTEX_Z - Lb_ENDVERTEX_Z > -1"
HLT1_CUT = "Lb_HLT1Track_TOS || Lb_HLT1TwoTrack_TOS"
PID_CUT = "(proton_ProbNNp > 0.4) && (proton_ProbNNk < 0.6) && (proton_ProbNNpi < 0.4) && (kaon_ProbNNpi < 0.5) && (kaon_ProbNNk > 0.4) && (mu_ProbNNmu > 0.5)"


def _read_window():
    if not CACHE_FILE.exists():
        raise FileNotFoundError(
            f"{CACHE_FILE} not found -- run sb_fit.py first")
    d = np.load(CACHE_FILE)
    return float(d["window_lo"]), float(d["window_hi"])


def _3pi_from_oof():
    rdf = ROOT.RDataFrame("tree", str(SCORED_TUPLE)).Filter("label==1")
    n_all = int(rdf.Count().GetValue())
    print(f"  SCORED_TUPLE label==1 entries: {n_all}")

    rdf = rdf.Filter(END_VZ_CUT)
    n_end_vz = int(rdf.Count().GetValue())
    store_count(n_end_vz, "n_end_vz", GEN_3PI)

    rdf = rdf.Filter(HLT1_CUT)
    n_hlt1 = int(rdf.Count().GetValue())
    store_count(n_hlt1, "n_hlt1", GEN_3PI)

    rdf = rdf.Filter(PID_CUT)
    n_pid = int(rdf.Count().GetValue())
    store_count(n_pid, "n_pid", GEN_3PI)

    rdf = run_veto_pipeline(rdf)
    rdf = rdf.Filter("passVetoes")
    n_veto = int(rdf.Count().GetValue())
    store_count(n_veto, "n_veto", GEN_3PI)

    cols = FEAT + ["Lb_DTF_M_aux", "bdt"]
    arr = rdf.AsNumpy(cols)
    df = pd.DataFrame({c: np.asarray(arr[c])
                       for c in cols}).dropna(subset=FEAT)
    df = df.rename(columns={"Lb_DTF_M_aux": "Lb_M"})
    return df[FEAT + ["Lb_M", "bdt"]]


def _load_gen(path):
    if not path.exists():
        raise FileNotFoundError(f"{path} not found (generator count)")
    return float(np.load(path)["n_gen"])


def score_frame(models, path, mass_cut, mass_cut2=None, gen=None):
    if not Path(path).exists():
        raise FileNotFoundError(f"skim not found: {path}")
    rdf = ROOT.RDataFrame("DecayTree", str(path)).Filter(mass_cut)

    if gen is not None:
        prdf = rdf.Filter(mass_cut2)

        n_vtx = int(prdf.Count().GetValue())
        store_count(n_vtx, "n_vtx", gen)

        prdf = prdf.Filter(END_VZ_CUT_R)
        n_end_vz = int(prdf.Count().GetValue())
        store_count(n_end_vz, "n_end_vz", gen)

        prdf = prdf.Filter(HLT1_CUT_R)
        n_hlt1 = int(prdf.Count().GetValue())
        store_count(n_hlt1, "n_hlt1", gen)

        prdf = prdf.Filter(PID_CUT_R)
        n_pid = int(prdf.Count().GetValue())
        store_count(n_pid, "n_pid", gen)

        prdf = run_veto_pipeline(prdf)
        prdf = prdf.Filter("passVetoes")
        n_veto = int(prdf.Count().GetValue())
        store_count(n_veto, "n_veto", gen)

        rdf = rdf.Filter(f"{END_VZ_CUT_R} && {HLT1_CUT_R} && {PID_CUT_R}")

    else:
        rdf = rdf.Filter(NEW_CUTS)

    rdf = run_veto_pipeline(rdf)
    rdf = rdf.Filter("passVetoes")

    df = rdf_to_frame(rdf, INPUT_FEATURES, AUX_FEATURES).dropna(subset=FEAT)
    df = df.rename(columns={"Lb_DTF_M_aux": "Lb_M"})
    df["bdt"] = ensemble_score(models, df[FEAT].to_numpy())
    return df[FEAT + ["Lb_M", "bdt"]]


def _snapshot(df, cols, out_path):
    out = {c: np.asarray(df[c].to_numpy(), dtype=np.float64) for c in cols}
    ROOT.RDF.FromNumpy(out).Snapshot("tree", str(out_path))
    print(f"  wrote {len(df):>7d} rows -> {out_path}")


# ---------------------------------------------------------------- main
def main():
    OUT_MC.parent.mkdir(parents=True, exist_ok=True)
    print(f"=== Build fit tuples ===")

    window_lo, window_hi = _read_window()
    print(f"signal window [{window_lo:.1f}, {window_hi:.1f}] MeV")
    print(f"full window (OS) [{x_min2:.1f}, {x_max2:.1f}] MeV")

    models = bu.load_models(str(MODEL_STEM), N_SPLITS)
    print(f"loaded {len(models)} fold models from {MODEL_STEM}_fold*.json")

    df1_cut1 = f"Lb_DTF_M > {x_min2} && Lb_DTF_M < {x_max2}"
    df1_cut2 = f"Lb_DTF_M > {window_lo} && Lb_DTF_M < {window_hi}"
    df0_cut = (f"(Lb_DTF_M > {window_hi} && Lb_DTF_M < {x_max2}) || "
               f"(Lb_DTF_M > {x_min2} && Lb_DTF_M < {window_lo})")

    # ---- signal MC, mode 0: 3pi nu ----
    print("\n[mode 0] tau -> 3pi nu  (OOF scores from SCORED_TUPLE)")
    df0 = _3pi_from_oof()
    n_gen0 = _load_gen(GEN_3PI)
    df0["mc_weight"] = BR_3PI / n_gen0
    df0["mode"] = 0
    print(f"  {len(df0)} events, n_gen={n_gen0:.0f}, "
          f"weight={BR_3PI / n_gen0:.4e}")

    # ---- signal MC, mode 0(b): 3pi nu ----
    print("\n[mode 0] tau -> 3pi nu  (Scored rest)")
    df0_b = score_frame(models, MC_3PI, df0_cut)
    df0_b["mc_weight"] = BR_3PI / n_gen0
    df0_b["mode"] = 0
    print(f"  {len(df0_b)} events")

    # ---- signal MC, mode 1: 3pi pi0 nu ----
    print("\n[mode 1] tau -> 3pi pi0 nu")
    df1 = score_frame(models, MC_3PIPI0, df1_cut1, df1_cut2, gen=GEN_3PIPI0)
    n_gen1 = _load_gen(GEN_3PIPI0)
    df1["mc_weight"] = BR_3PIPI0 / n_gen1
    df1["mode"] = 1
    print(f"  {len(df1)} events, n_gen={n_gen1:.0f}, "
          f"weight={BR_3PIPI0 / n_gen1:.4e}")

    # ---- combine + write ----
    mc = pd.concat([df0, df0_b, df1], ignore_index=True)
    print("\ncombined signal MC:")
    _snapshot(mc, FEAT + ["bdt", "Lb_M", "mc_weight", "mode"], OUT_MC)

    # ---- OS data (sidebands only, unless unblinding) ----
    print(f"\n[OS data] region = {OS_REGION}")
    if OS_REGION == "sideband":
        os_cut = (f"(Lb_DTF_M > {window_hi} && Lb_DTF_M < {x_max2}) || "
                  f"(Lb_DTF_M > {x_min2} && Lb_DTF_M < {window_lo})")
    elif OS_REGION == "full":
        os_cut = f"Lb_DTF_M > {x_min2} && Lb_DTF_M < {x_max2}"
        print("  WARNING: signal box is UNBLINDED")
    dos = score_frame(models, DAT_OS, os_cut)
    _snapshot(dos, FEAT + ["bdt", "Lb_M"], OUT_OS)

    print("\ndone.")


if __name__ == "__main__":
    main()
