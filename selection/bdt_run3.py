from pathlib import Path

import numpy as np
import pandas as pd
import ROOT

import bdt_util as bu

skim_dir = Path(
    "/afs/cern.ch/user/a/adalmia/workingDirectory/vrd-lb2pktaumu/run3/skims")

LB_PDG_M = 5619.5
NSIGMA = 3.0
x_min, x_max = 4500, 8000
BKG_MODE = "ss"  #os, combined, ss

INPUT_FEATURES = [
    ("Lb_FD_CHI2", "Lb_BPVFDCHI2"),
    ("Lb_CHI2_DOF", "Lb_CHI2DOF"),
    ("Lb_pT", "Lb_PT"),
    ("Lb_DIRA", "Lb_BPVDIRA"),
    ("L1520_FD_CHI2", "L1520_BPVFDCHI2"),
    ("L1520_pT", "L1520_PT"),
    ("proton_IP_CHI2", "proton_BPVIPCHI2"),
    ("proton_pT", "proton_PT"),
    ("kaon_IP_CHI2", "kaon_BPVIPCHI2"),
    ("Jpsi_pT", "Jpsi_PT"),
    ("Lb_CONE_C_0p5_MaxP", "Lb_HEAD_CC_B_ChargedIso_DR2_0p5_Max_P"),
    ("Lb_CONE_C_0p5_PTASY", "Lb_HEAD_CC_B_ChargedIso_DR2_0p5_PTASY"),
    ("Lb_CONE_N_0p5_PTASY", "Lb_HEAD_NC_B_NeutralIso_DR2_0p5_PTASY"),
    ("Lb_CONE_C_0p5_CMULT", "Lb_HEAD_CC_B_ChargedIso_DR2_0p5_CMULT"),
]

AUX_FEATURES = [("Lb_DTF_M_aux", "Lb_DTF_M"), ("tau_ENDVERTEX_Z",
                                               "tau_END_VZ"),
                ("Lb_ENDVERTEX_Z", "Lb_END_VZ"),
                ("Lb_HLT1Track_TOS", "Lb_Hlt1TrackMVADecision_TOS"),
                ("Lb_HLT1TwoTrack_TOS", "Lb_Hlt1TwoTrackMVADecision_TOS"),
                ("proton_ProbNNp", "proton_PROBNN_P"),
                ("proton_ProbNNk", "proton_PROBNN_K"),
                ("proton_ProbNNpi", "proton_PROBNN_PI"),
                ("kaon_ProbNNk", "kaon_PROBNN_K"),
                ("kaon_ProbNNpi", "kaon_PROBNN_PI"),
                ("mu_ProbNNmu", "mu_PROBNN_MU")]

VETO_INPUTS = [
    f"{p}_{c}" for p in ["proton", "kaon", "pi1", "pi2", "pi3", "mu", "tau"]
    for c in ["PX", "PY", "PZ"]
] + ["proton_CHARGE", "mu_CHARGE"]
VETO_AUX = [(f"{b}__raw", b) for b in VETO_INPUTS]

N_SPLITS = 5
SEED = 42
RUN_OPTUNA = False

OUT_DIR = Path("./bdt_work")
MODEL_STEM = OUT_DIR / "models" / f"lb_pktaumu_xgb_{BKG_MODE}"
PLOT_DIR = OUT_DIR / "plots" / BKG_MODE
SCORED_TUPLE = OUT_DIR / "cache" / f"lb_pktaumu_scored_{BKG_MODE}.root"


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------
def rms90(histo, imean=None):
    """RMS of the central 90% -- same definition as the study scripts."""
    import math
    nbins = histo.GetNbinsX()
    if imean is None:
        imean = histo.GetMaximumBin()
    entries = 0.9 * histo.GetEntries()
    sumw = histo.GetBinContent(imean)
    sumwx = sumw * histo.GetBinCenter(imean)
    sumwx2 = sumw * histo.GetBinCenter(imean)**2
    for i in range(1, nbins):
        if sumw >= entries:
            break
        if imean - i > 0:
            w = histo.GetBinContent(imean - i)
            x = histo.GetBinCenter(imean - i)
            sumw += w
            sumwx += w * x
            sumwx2 += w * x * x
        if imean + i <= nbins:
            w = histo.GetBinContent(imean + i)
            x = histo.GetBinCenter(imean + i)
            sumw += w
            sumwx += w * x
            sumwx2 += w * x * x
    mean90 = sumwx / sumw
    return math.sqrt(abs(sumwx2 / sumw - mean90 * mean90))


def rdf_to_frame(rdf, feature_defs, aux_defs):
    """Define every feature/aux expression, then dump to a pandas frame.
    Everything is cast to double so XGBoost sees clean floats."""
    node = rdf
    for name, expr in feature_defs + aux_defs:
        node = node.Define(name, f"(double)({expr})")
    cols = [n for n, _ in feature_defs] + [n for n, _ in aux_defs]
    arr = node.AsNumpy(cols)
    df = pd.DataFrame({c: np.asarray(arr[c]) for c in cols})
    df["Lb_CONE_C_0p5_MaxP"] = df["Lb_CONE_C_0p5_MaxP"].fillna(0)
    return df


# ----------------------------------------------------------------------
def main():
    for d in [OUT_DIR / "models", PLOT_DIR, OUT_DIR / "cache"]:
        d.mkdir(parents=True, exist_ok=True)

    sample = "taumu"

    # --- read the cached skims (selection already applied) --------------
    sim_base = ROOT.RDataFrame("DecayTree",
                               str(skim_dir / f"{sample}_sim_OS.root"))
    ss_base = ROOT.RDataFrame("DecayTree",
                              str(skim_dir / f"{sample}_dat_SS.root"))
    os_base = ROOT.RDataFrame("DecayTree",
                              str(skim_dir / f"{sample}_dat_OS.root"))

    # --- find the DTF peak + RMS90 from sim, define the windows ---------
    h = sim_base.Filter(f"Lb_DTF_M > {x_min} && Lb_DTF_M < {x_max}").Histo1D(
        ROOT.RDF.TH1DModel("h_fine", "", 3500, x_min, x_max), "Lb_DTF_M")
    m_peak = h.GetBinCenter(h.GetMaximumBin())
    sigma = rms90(h.GetValue())
    window_lo = m_peak - NSIGMA * sigma
    window_hi = m_peak + NSIGMA * sigma
    print(
        f"\nDTF peak {m_peak:.1f} MeV (PDG {LB_PDG_M}), RMS90 {sigma:.1f} MeV")
    print(f"signal window [{window_lo:.1f}, {window_hi:.1f}] MeV")

    signal_window_cut = f"Lb_DTF_M > {window_lo} && Lb_DTF_M < {window_hi}"
    upper_sideband_cut = f"Lb_DTF_M > {window_hi} && Lb_DTF_M < {x_max}"
    full_range_cut = f"Lb_DTF_M > {x_min} && Lb_DTF_M < {x_max}"

    # --- signal -----------------------------------------------------------
    sig_frame = rdf_to_frame(
        sim_base.Filter(signal_window_cut), INPUT_FEATURES,
        AUX_FEATURES + VETO_AUX)
    sig_frame["label"] = 1
    sig_frame["source"] = -1  # -1 = signal
    print(f"signal candidates: {len(sig_frame)}")

    # --- background: SS (full range) + OS upper sideband ----------------
    # NB: os_base is the full OS spectrum incl. the blinded region. We ONLY
    # ever touch it through upper_sideband_cut -- never loosen this.
    ss_frame = rdf_to_frame(
        ss_base.Filter(full_range_cut), INPUT_FEATURES,
        AUX_FEATURES + VETO_AUX)
    ss_frame["label"] = 0
    ss_frame["source"] = 0  # 0 = SS data
    #print(f"SS background candidates: {len(ss_frame)}")

    os_frame = rdf_to_frame(
        os_base.Filter(upper_sideband_cut), INPUT_FEATURES,
        AUX_FEATURES + VETO_AUX)
    os_frame["label"] = 0
    os_frame["source"] = 1  # 1 = OS sideband
    #print(f"OS sideband background candidates: {len(os_frame)}")
    '''
    # --- combine + weights ----------------------------------------------
    # Balance classes (total bkg weight == total sig weight) and balance the
    # two background proxies (SS and OS each get half the background budget).
    bkg_frames = {"ss": [ss_frame], "os": [os_frame],
              "combined": [ss_frame, os_frame]}[BKG_MODE]
    data = pd.concat([sig_frame, *bkg_frames],
                     ignore_index=True).dropna()
    n_sig = (data.label == 1).sum()

    present = [s for s in (0, 1) if (data.source == s).any()]   # 0 = SS, 1 = OS
    per_source = n_sig / len(present)

    #n_ss = ((data.label == 0) & (data.source == 0)).sum()
    #n_os = ((data.label == 0) & (data.source == 1)).sum()

    w = np.ones(len(data), dtype=float)

    for src in present:
        n_src = ((data.label == 0) & (data.source == src)).sum()
        if n_src > 0:
            w[(data.label == 0) & (data.source == src)] = per_source / n_src

    
   # if n_ss > 0:
     #   w[(data.label == 0) & (data.source == 0)] = (n_sig / 2.0) / n_ss
   # if n_os > 0:
     #   w[(data.label == 0) & (data.source == 1)] = (n_sig / 2.0) / n_os
    

    data["weight"] = w

    feat_names = [n for n, _ in INPUT_FEATURES]
    X = data[feat_names].to_numpy()
    y = data["label"].to_numpy()
    weights = data["weight"].to_numpy()

    if BKG_MODE == "combined":
        strat = (data.label.to_numpy() * 10 + (data.source.to_numpy() + 1))
    else:
        strat = data.label.to_numpy()
    '''

    N_BKG_PER_SIG = 2.0  # background events drawn per signal event
    BKG_WEIGHT = 0.5  # weight given to every drawn background event

    feat_names = [n for n, _ in INPUT_FEATURES]

    bkg_frames = {
        "ss": [ss_frame],
        "os": [os_frame],
        "combined": [ss_frame, os_frame]
    }[BKG_MODE]
    bkg_pool = pd.concat(
        bkg_frames, ignore_index=True).dropna(subset=feat_names)
    sig_frame_raw = sig_frame
    sig_frame = sig_frame.dropna(subset=feat_names)

    for name in feat_names:
        n_nan = sig_frame_raw[name].isna().sum()  # the frame BEFORE dropna
        print(f"{name:<20} {n_nan:6d}  ({100*n_nan/len(sig_frame_raw):.1f}%)")

    n_sig = len(sig_frame)
    n_draw = int(N_BKG_PER_SIG * n_sig)
    if len(bkg_pool) < n_draw:
        print(f"WARNING: only {len(bkg_pool)} bkg available, wanted {n_draw}; "
              f"using all -> classes will not be exactly balanced.")
        n_draw = len(bkg_pool)
    bkg_draw = bkg_pool.sample(n=n_draw, random_state=SEED)

    data = pd.concat([sig_frame, bkg_draw], ignore_index=True)
    data["weight"] = np.where(data.label == 1, 1.0, BKG_WEIGHT)
    print(f"training set: {n_sig} signal (w=1.0) + {n_draw} background "
          f"(w={BKG_WEIGHT}) [{BKG_MODE}]")

    X = data[feat_names].to_numpy()
    y = data["label"].to_numpy()
    weights = data["weight"].to_numpy()

    # stratify on source only when both backgrounds are present
    if BKG_MODE == "combined":
        strat = data.label.to_numpy() * 10 + (data.source.to_numpy() + 1)
    else:
        strat = data.label.to_numpy()

    params = None

    # --- train ----------------------------------------------------------
    models, oof, fold_id = bu.train_kfold(
        X,
        y,
        weights,
        n_splits=N_SPLITS,
        params=params,
        strat_key=strat,
        seed=SEED)
    data["bdt"] = oof
    data["fold"] = fold_id

    # --- diagnostics ----------------------------------------------------
    bu.plot_roc(y, oof, weights, PLOT_DIR / "roc.pdf", fold_id=fold_id)
    bu.plot_importance(models, feat_names, PLOT_DIR / "importance.pdf")
    bu.correlation_matrix(
        data[data.label == 1],
        feat_names + ["bdt"],
        PLOT_DIR / "corr_sig.pdf",
        title="signal")
    bu.correlation_matrix(
        data[data.label == 0],
        feat_names + ["bdt"],
        PLOT_DIR / "corr_bkg.pdf",
        title="background")
    bu.response_by_source(data["bdt"], data["label"], data["source"],
                          PLOT_DIR / "response_by_source.pdf")

    bkg = data[data.label == 0]
    #bu.sculpting_check(bkg["Lb_DTF_M_aux"], bkg["bdt"],
    #                  PLOT_DIR / "sculpting.pdf", weights_bkg=bkg["weight"])

    # --- save models + scored tuple -------------------------------------
    bu.save_models(models, str(MODEL_STEM))
    out_cols = {
        c: np.asarray(data[c].to_numpy(), dtype=np.float64)
        for c in feat_names + [
            "Lb_DTF_M_aux", "tau_ENDVERTEX_Z", "Lb_ENDVERTEX_Z",
            "Lb_HLT1Track_TOS", "Lb_HLT1TwoTrack_TOS", "proton_ProbNNp",
            "proton_ProbNNk", "proton_ProbNNpi", "kaon_ProbNNk",
            "kaon_ProbNNpi", "mu_ProbNNmu", "label", "source", "weight",
            "fold", "bdt"
        ]
    }
    for b in VETO_INPUTS:
        out_cols[b] = np.asarray(
            data[f"{b}__raw"].to_numpy(), dtype=np.float64)
    ROOT.RDF.FromNumpy(out_cols).Snapshot("tree", str(SCORED_TUPLE))
    print(f"saved scored tuple -> {SCORED_TUPLE}")
    bu.bdt_mass_correlation(SCORED_TUPLE, PLOT_DIR / "bdt_corr.pdf")

    print("\nDone. Diagnostics in", PLOT_DIR)


if __name__ == "__main__":
    main()
