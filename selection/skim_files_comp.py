import ROOT
import numpy as np
from pathlib import Path

from simplePlots import get_all_input_files
from punzi_scan import load_cache

# Multi-core: the skim is a pure read+filter+write, parallelises perfectly
ROOT.EnableImplicitMT()

friend_dir = Path("/eos/lhcb/user/f/fabudine/vrd-lb2pktaumu/run3/friends")

# Local output dir for the skims - put this somewhere on your work area,
# NOT on EOS, so reruns read locally.
skim_dir = Path(
    "/afs/cern.ch/user/a/adalmia/workingDirectory/vrd-lb2pktaumu/run3/skims")
skim_dir.mkdir(exist_ok=True)

chain_prefix = "DTF3"

gen_tree_paths = ["MCDT_tuple1/MCDecayTree", "MCDT_tuple2/MCDecayTree"]

BR_TAU = 0.0931
BR_TAU_PI = 0.0462

sample_configs = {
    "taumu": {
        "sim_sample_name": "Lb_pKtaumu_3pi",
        "sim_pi_sample_name": "Lb_pKtaumu_3pipi0",
        "dat_sample_name": "turbo_rd",
        "sim_chain_name": "Hlt2RD_LbToPKTauMu_TauTo3Pi_OS",
        "dat_chain_name": "Hlt2RD_LbToPKTauMu_TauTo3Pi_SS",
        "lepton": "mu",
        "intpart": "L1520",
        "leptonID": "13",
    }
}

# Columns to keep in the skim: everything used for cuts, the mass window,
# and every candidate BDT variable. Keep this in sync with the variables
# dict in the comparison script.
keep_columns = [
    "Lb_DTF_M",
    "Lb_DTF_CHI2DOF",
    "Lb_Hlt1TrackMVADecision_TOS",
    "Lb_Hlt1TwoTrackMVADecision_TOS",
    "Lb_Hlt1TrackMuonMVADecision_TOS",
    "proton_PROBNN_P",
    "proton_PROBNN_K",
    "proton_PROBNN_PI",
    "kaon_PROBNN_K",
    "kaon_PROBNN_PI",
    "mu_PROBNN_MU",
    "tau_END_VZ",
    "Lb_END_VZ",
    #candidates
    "proton_P",
    "proton_PT",
    "proton_PID_P",
    "proton_BPVIPCHI2",
    "proton_PID_K",
    "kaon_PT",
    "kaon_P",
    "kaon_PID_K",
    "kaon_BPVIPCHI2",
    "L1520_M",
    "L1520_PT",
    "L1520_CHI2DOF",
    "L1520_BPVFDCHI2",
    "L1520_MAXSDOCACHI2",
    "L1520_P",
    "L1520_BPVIPCHI2",
    "L1520_BPVDIRA",
    "mu_PT",
    "mu_P",
    "mu_PID_MU",
    "mu_BPVIPCHI2",
    "mu_CHI2DOF",
    "tau_M",
    "tau_PT",
    "tau_MAXSDOCACHI2",
    "tau_BPVFDCHI2",
    "tau_BPVDIRA",
    "tau_P",
    "tau_CHI2DOF",
    "tau_DOCACHI2_13",
    "tau_DOCACHI2_23",
    "tau_MAX_PT",
    "tau_MAX_BPVIPCHI2",
    "pi1_PT",
    "pi1_P",
    "pi1_PID_K",
    "pi1_PID_E",
    "pi1_PID_P",
    "pi1_PID_MU",
    "pi1_BPVIPCHI2",
    "pi1_CHI2DOF",
    "pi2_PT",
    "pi2_P",
    "pi2_PID_K",
    "pi2_PID_E",
    "pi2_PID_P",
    "pi2_PID_MU",
    "pi2_BPVIPCHI2",
    "pi2_CHI2DOF",
    "pi3_PT",
    "pi3_P",
    "pi3_PID_K",
    "pi3_PID_E",
    "pi3_PID_P",
    "pi3_PID_MU",
    "pi3_BPVIPCHI2",
    "pi3_CHI2DOF",
    "Jpsi_M",
    "Jpsi_CHI2DOF",
    "Jpsi_BPVIPCHI2",
    "Jpsi_P",
    "Jpsi_PT",
    "Jpsi_BPVFDCHI2",
    "Jpsi_BPVDIRA",
    "Jpsi_MAXSDOCACHI2",
    "Lb_M",
    "Lb_CHI2DOF",
    "Lb_P",
    "Lb_PT",
    "Lb_BPVIPCHI2",
    "Lb_BPVFDCHI2",
    "Lb_BPVDIRA",
    "Lb_MAXSDOCACHI2",
    "Lb_HEAD_CC_B_ChargedIso_DR2_0p5_PTASY",
    "Lb_HEAD_CC_B_ChargedIso_DR2_0p5_PASY",
    "Lb_HEAD_CC_B_ChargedIso_DR2_0p5_Max_P",
    "Lb_HEAD_CC_B_ChargedIso_DR2_0p5_CMULT",
    "Lb_HEAD_CC_B_ChargedIso_DR2_0p5_Max_PT",
    "Lb_HEAD_NC_B_NeutralIso_DR2_0p5_CMULT",
    "Lb_HEAD_NC_B_NeutralIso_DR2_0p5_PTASY",
    "Lb_HEAD_NC_B_NeutralIso_DR2_0p5_PASY",
    "Lb_HEAD_NC_B_NeutralIso_DR2_0p5_Max_P",
    "Lb_HEAD_NC_B_NeutralIso_DR2_0p5_Max_PT",
    #--vetoes--
    "proton_PX",
    "proton_PY",
    "proton_PZ",
    "kaon_PX",
    "kaon_PY",
    "kaon_PZ",
    "pi1_PX",
    "pi1_PY",
    "pi1_PZ",
    "pi2_PX",
    "pi2_PY",
    "pi2_PZ",
    "pi3_PX",
    "pi3_PY",
    "pi3_PZ",
    "mu_PX",
    "mu_PY",
    "mu_PZ",
    "tau_PX",
    "tau_PY",
    "tau_PZ",
    "proton_CHARGE",
    "mu_CHARGE",
]
'''
    #Lb_CISO
    "Lb_HEAD_CC_B_ChargedIso_DR2_0p05_CMULT", "Lb_HEAD_CC_B_ChargedIso_DR2_0p05_PTASY", "Lb_HEAD_CC_B_ChargedIso_DR2_0p05_PASY", "Lb_HEAD_CC_B_ChargedIso_DR2_0p05_Max_P", "Lb_HEAD_CC_B_ChargedIso_DR2_0p05_Max_PT",
    "Lb_HEAD_NC_B_NeutralIso_DR2_0p05_CMULT", "Lb_HEAD_NC_B_NeutralIso_DR2_0p05_PTASY", "Lb_HEAD_NC_B_NeutralIso_DR2_0p05_PASY", "Lb_HEAD_NC_B_NeutralIso_DR2_0p05_Max_P", "Lb_HEAD_NC_B_NeutralIso_DR2_0p05_Max_PT",
    "Lb_HEAD_CC_B_ChargedIso_DR2_0p1_CMULT", "Lb_HEAD_CC_B_ChargedIso_DR2_0p1_PTASY", "Lb_HEAD_CC_B_ChargedIso_DR2_0p1_PASY", "Lb_HEAD_CC_B_ChargedIso_DR2_0p1_Max_P", "Lb_HEAD_CC_B_ChargedIso_DR2_0p1_Max_PT",
    "Lb_HEAD_NC_B_NeutralIso_DR2_0p1_CMULT", "Lb_HEAD_NC_B_NeutralIso_DR2_0p1_PTASY", "Lb_HEAD_NC_B_NeutralIso_DR2_0p1_PASY", "Lb_HEAD_NC_B_NeutralIso_DR2_0p1_Max_P", "Lb_HEAD_NC_B_NeutralIso_DR2_0p1_Max_PT",
    "Lb_HEAD_CC_B_ChargedIso_DR2_0p15_CMULT", "Lb_HEAD_CC_B_ChargedIso_DR2_0p15_PTASY", "Lb_HEAD_CC_B_ChargedIso_DR2_0p15_PASY", "Lb_HEAD_CC_B_ChargedIso_DR2_0p15_Max_P", "Lb_HEAD_CC_B_ChargedIso_DR2_0p15_Max_PT",
    "Lb_HEAD_NC_B_NeutralIso_DR2_0p15_CMULT", "Lb_HEAD_NC_B_NeutralIso_DR2_0p15_PTASY", "Lb_HEAD_NC_B_NeutralIso_DR2_0p15_PASY", "Lb_HEAD_NC_B_NeutralIso_DR2_0p15_Max_P", "Lb_HEAD_NC_B_NeutralIso_DR2_0p15_Max_PT",
    "Lb_HEAD_CC_B_ChargedIso_DR2_0p2_CMULT", "Lb_HEAD_CC_B_ChargedIso_DR2_0p2_PTASY", "Lb_HEAD_CC_B_ChargedIso_DR2_0p2_PASY", "Lb_HEAD_CC_B_ChargedIso_DR2_0p2_Max_P", "Lb_HEAD_CC_B_ChargedIso_DR2_0p2_Max_PT",
    "Lb_HEAD_NC_B_NeutralIso_DR2_0p2_CMULT", "Lb_HEAD_NC_B_NeutralIso_DR2_0p2_PTASY", "Lb_HEAD_NC_B_NeutralIso_DR2_0p2_PASY", "Lb_HEAD_NC_B_NeutralIso_DR2_0p2_Max_P", "Lb_HEAD_NC_B_NeutralIso_DR2_0p2_Max_PT",
    "Lb_HEAD_CC_B_ChargedIso_DR2_0p5_CMULT", "Lb_HEAD_CC_B_ChargedIso_DR2_0p5_PTASY", "Lb_HEAD_CC_B_ChargedIso_DR2_0p5_PASY", "Lb_HEAD_CC_B_ChargedIso_DR2_0p5_Max_P", "Lb_HEAD_CC_B_ChargedIso_DR2_0p5_Max_PT",
    "Lb_HEAD_NC_B_NeutralIso_DR2_0p5_CMULT", "Lb_HEAD_NC_B_NeutralIso_DR2_0p5_PTASY", "Lb_HEAD_NC_B_NeutralIso_DR2_0p5_PASY", "Lb_HEAD_NC_B_NeutralIso_DR2_0p5_Max_P", "Lb_HEAD_NC_B_NeutralIso_DR2_0p5_Max_PT"]


    #Lb_ISO
    "tau_VTXISO_tau_OneTrack_DChi2Vtx_3p0_NParts", "tau_VTXISO_tau_OneTrack_DChi2Vtx_3p0_Smallest_DELTACHI2", "tau_VTXISO_tau_OneTrack_DChi2Vtx_3p0_Smallest_CHI2", "tau_VTXISO_tau_OneTrack_DChi2Vtx_3p0_Smallest_DELTACHI2_MASS",
    "tau_VTXISO_tau_TwoTracks_DChi2Vtx_3p0_NParts", "tau_VTXISO_tau_TwoTracks_DChi2Vtx_3p0_Smallest_DELTACHI2", "tau_VTXISO_tau_TwoTracks_DChi2Vtx_3p0_Smallest_CHI2", "tau_VTXISO_tau_TwoTracks_DChi2Vtx_3p0_Smallest_DELTACHI2_MASS",
    "tau_VTXISO_tau_OneTrack_DChi2Vtx_9p0_NParts", "tau_VTXISO_tau_OneTrack_DChi2Vtx_9p0_Smallest_DELTACHI2", "tau_VTXISO_tau_OneTrack_DChi2Vtx_9p0_Smallest_CHI2", "tau_VTXISO_tau_OneTrack_DChi2Vtx_9p0_Smallest_DELTACHI2_MASS",
    "tau_VTXISO_tau_TwoTracks_DChi2Vtx_9p0_NParts", "tau_VTXISO_tau_TwoTracks_DChi2Vtx_9p0_Smallest_DELTACHI2", "tau_VTXISO_tau_TwoTracks_DChi2Vtx_9p0_Smallest_CHI2", "tau_VTXISO_tau_TwoTracks_DChi2Vtx_9p0_Smallest_DELTACHI2_MASS", 
    "tau_VTXISO_tau_OneTrack_DChi2Vtx_25p0_NParts", "tau_VTXISO_tau_OneTrack_DChi2Vtx_25p0_Smallest_DELTACHI2", "tau_VTXISO_tau_OneTrack_DChi2Vtx_25p0_Smallest_CHI2", "tau_VTXISO_tau_OneTrack_DChi2Vtx_25p0_Smallest_DELTACHI2_MASS",
    "tau_VTXISO_tau_TwoTracks_DChi2Vtx_25p0_NParts", "tau_VTXISO_tau_TwoTracks_DChi2Vtx_25p0_Smallest_DELTACHI2", "tau_VTXISO_tau_TwoTracks_DChi2Vtx_25p0_Smallest_CHI2", "tau_VTXISO_tau_TwoTracks_DChi2Vtx_25p0_Smallest_DELTACHI2_MASS",
    "tau_VTXISO_tau_OneTrack_DChi2Vtx_100p0_NParts", "tau_VTXISO_tau_OneTrack_DChi2Vtx_100p0_Smallest_DELTACHI2", "tau_VTXISO_tau_OneTrack_DChi2Vtx_100p0_Smallest_CHI2", "tau_VTXISO_tau_OneTrack_DChi2Vtx_100p0_Smallest_DELTACHI2_MASS",
    "tau_VTXISO_tau_TwoTracks_DChi2Vtx_100p0_NParts", "tau_VTXISO_tau_TwoTracks_DChi2Vtx_100p0_Smallest_DELTACHI2", "tau_VTXISO_tau_TwoTracks_DChi2Vtx_100p0_Smallest_CHI2", "tau_VTXISO_tau_TwoTracks_DChi2Vtx_100p0_Smallest_DELTACHI2_MASS",
    "tau_VTXISO_tau_OneTrack_DChi2Vtx_NO_NParts", "tau_VTXISO_tau_OneTrack_DChi2Vtx_NO_Smallest_DELTACHI2", "tau_VTXISO_tau_OneTrack_DChi2Vtx_NO_Smallest_CHI2", "tau_VTXISO_tau_OneTrack_DChi2Vtx_NO_Smallest_DELTACHI2_MASS",
    "tau_VTXISO_tau_TwoTracks_DChi2Vtx_NO_NParts", "tau_VTXISO_tau_TwoTracks_DChi2Vtx_NO_Smallest_DELTACHI2", "tau_VTXISO_tau_TwoTracks_DChi2Vtx_NO_Smallest_CHI2", "tau_VTXISO_tau_TwoTracks_DChi2Vtx_NO_Smallest_DELTACHI2_MASS"
]
'''

# Truth-match columns are only present in simulation
truth_columns = [
    "pi1_TRUEID",
    "pi2_TRUEID",
    "pi3_TRUEID",
    "pi1_MC_MOTHER_ID",
    "pi1_MC_MOTHER_KEY",
    "pi2_MC_MOTHER_KEY",
    "pi3_MC_MOTHER_KEY",
    "proton_TRUEID",
    "kaon_TRUEID",
    "mu_TRUEID",
    "proton_MC_MOTHER_ID",
    "proton_MC_MOTHER_KEY",
    "kaon_MC_MOTHER_KEY",
    "mu_MC_MOTHER_KEY",
    "pi1_MC_GD_MOTHER_KEY",
    "Lb_BKGCAT",
]


def build_chain(sample_name, chain_name):
    chain = ROOT.TChain(f"{chain_name}_{chain_prefix}/DecayTree")
    friend_chain = ROOT.TChain("recomass")
    for sname, in_file in get_all_input_files(sample_name):
        in_file_stem = Path(in_file).name.split('.')[0]
        chain.Add(in_file)
        friend_file = friend_dir / f"recomass_{chain_name}_{chain_prefix}_{in_file_stem}.root"
        friend_chain.Add(str(friend_file))
    chain.AddFriend(friend_chain)
    return chain, friend_chain


def count_gen_entries(sample_name, tree_paths):
    total = 0
    for tp in tree_paths:
        chain = ROOT.TChain(tp)
        for sname, in_file in get_all_input_files(sample_name):
            chain.AddFile(in_file)
        total += int(chain.GetEntries())
    return total


def skim(sample_name, chain_name, cut, out_path, columns, weight=None):
    """Apply cut, materialize the kept columns, write a small local file."""
    chain, friend = build_chain(sample_name, chain_name)
    rdf = ROOT.RDataFrame(chain).Filter(cut)

    if weight is not None:
        rdf = rdf.Define("Lb_tau_BR", f"{weight!r}")
        columns.append("Lb_tau_BR")

    # Only keep columns that actually exist in this chain
    available = set(str(c) for c in rdf.GetColumnNames())
    cols = [c for c in columns if c in available]
    missing = [c for c in columns if c not in available]
    if missing:
        print(f"  (skipping {len(missing)} absent columns: {missing})")

    print(f"Writing skim -> {out_path}")
    rdf.Snapshot("DecayTree", str(out_path), cols)


if __name__ == '__main__':

    for sample in ["taumu"]:

        cfg = sample_configs[sample]
        lepton = cfg["lepton"]
        intpart = cfg["intpart"]
        leptonID = cfg["leptonID"]

        # --- cuts (same definitions as the comparison script) ---
        tm_cut = f"abs(pi1_TRUEID)==211 && abs(pi2_TRUEID)==211 && abs(pi3_TRUEID)==211 && abs(pi1_MC_MOTHER_ID)==15 && pi1_MC_MOTHER_KEY==pi2_MC_MOTHER_KEY && pi1_MC_MOTHER_KEY== pi3_MC_MOTHER_KEY && abs(proton_TRUEID)==2212 && abs(kaon_TRUEID)==321 && abs({lepton}_TRUEID)=={leptonID} && abs(proton_MC_MOTHER_ID)==5122 && proton_MC_MOTHER_KEY==kaon_MC_MOTHER_KEY && proton_MC_MOTHER_KEY=={lepton}_MC_MOTHER_KEY && proton_MC_MOTHER_KEY==pi1_MC_GD_MOTHER_KEY"

        convergence_cut = "Lb_DTF_CHI2DOF < 50"
        Lb_cut = "Lb_CHI2DOF < 20 && Lb_HEAD_CC_B_ChargedIso_DR2_0p5_PTASY > -0.5 && Lb_BPVDIRA > 0.9995 && Lb_HEAD_NC_B_NeutralIso_DR2_0p5_PTASY > -0.5"

        sim_cut = " && ".join([
            tm_cut,
            Lb_cut,  # initial_cuts,
            convergence_cut
        ])
        dat_cut = " && ".join([
            Lb_cut,  # pid_cut, initial_cuts,
            convergence_cut
        ])

        sim_sample_name = cfg["sim_sample_name"]
        sim_pi_sample_name = cfg["sim_pi_sample_name"]
        dat_sample_name = cfg["dat_sample_name"]
        sim_chain_name = cfg["sim_chain_name"]  # OS line
        dat_chain_name = cfg["dat_chain_name"]  # SS line

        # Make skim files (run first 1 out of 2)

        n_gen = count_gen_entries(sim_sample_name, gen_tree_paths)
        print(f"Total generated signal decays ({sample}): {n_gen}")

        n_gen_pi = count_gen_entries(sim_pi_sample_name, gen_tree_paths)
        print(
            f"Total generated signal decays ({sample}) (with pi0): {n_gen_pi}")

        w_mc = BR_TAU / n_gen
        print(f"MC weight = {BR_TAU} / {n_gen} = {w_mc:.4e}")

        w_mc_pi = BR_TAU_PI / n_gen_pi
        print(
            f"MC weight (with pi0) = {BR_TAU_PI} / {n_gen_pi} = {w_mc_pi:.4e}")

        # 1a. Signal: simulation, OS line, truth-matched
        skim(
            sim_sample_name,
            sim_chain_name,
            sim_cut,
            skim_dir / f"{sample}_sim_OS.root",
            keep_columns + truth_columns,
            weight=w_mc)

        # 1b. Signal: simulation, OS line, truth-matched, pi0
        skim(
            sim_pi_sample_name,
            sim_chain_name,
            sim_cut,
            skim_dir / f"{sample}_pi0_sim_OS.root",
            keep_columns + truth_columns,
            weight=w_mc_pi)

        # 2. Background proxy A: data, SS line (full spectrum)
        skim(dat_sample_name, dat_chain_name, dat_cut,
             skim_dir / f"{sample}_dat_SS.root", keep_columns)

        # 3. Background proxy B: data, OS line (sideband applied later)
        skim(dat_sample_name, sim_chain_name, dat_cut,
             skim_dir / f"{sample}_dat_OS.root", keep_columns)

        # Make skim files (run second 2 out of 2)

        print(f"\n Counting entries for {sim_sample_name}...")

        sim_chain_c, sim_friend_c = build_chain(sim_sample_name,
                                                sim_chain_name)
        n_reco = sim_chain_c.GetEntries()
        print(f"Reconstructed candidates before sim_cut ({sample}): {n_reco}")

        _, _, _, window_lo, window_hi = load_cache()
        n_win = sim_chain_c.GetEntries(
            f"Lb_DTF_M > {window_lo} && Lb_DTF_M < {window_hi}")
        print(f"Candidates in signal window ({sample}): {n_win}")

        gen_out = skim_dir / f"{sample}_gen_counts.npz"
        np.savez(gen_out, n_gen=n_gen, n_reco=n_reco, n_win=n_win)

        print(f"Saved counts -> {gen_out}")

        print(f"\n Counting entries for {sim_pi_sample_name}...")

        sim_chain_p, sim_friend_p = build_chain(sim_pi_sample_name,
                                                sim_chain_name)
        n_reco_pi = sim_chain_p.GetEntries()
        print(
            f"Reconstructed candidates before sim_cut ({sample}) (with pi0): {n_reco_pi}"
        )

        n_win_pi = sim_chain_p.GetEntries(
            f"Lb_DTF_M > {window_lo} && Lb_DTF_M < {window_hi}")
        print(f"Candidates in signal window ({sample}) (with pi0): {n_win_pi}")

        gen_out_pi = skim_dir / f"{sample}_pi_gen_counts.npz"
        np.savez(gen_out_pi, n_gen=n_gen_pi, n_reco=n_reco_pi, n_win=n_win_pi)

        print(f"Saved counts -> {gen_out_pi}")

    print("\nSkims done.")
