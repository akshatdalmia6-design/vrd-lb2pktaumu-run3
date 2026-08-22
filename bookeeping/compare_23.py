import ROOT
import sys
import math
import numpy as np
from run3_samples_v1r5791 import samples
from pathlib import Path
import glob

from simplePlots import get_all_input_files
sys.path.append("../utils")
from lhcbPlotStyle import setLHCbPlotStyle

setLHCbPlotStyle()

friend_dir = Path("  ")

run2_inputDir = Path(
    " ")

HLT1_cut_mu = "(Lb_Hlt1TrackMVADecision_TOS || Lb_Hlt1TwoTrackMVADecision_TOS || Lb_Hlt1TrackMuonMVADecision_TOS)"

# Common cuts

belowD_cut = "tau_M>800 && tau_M<1500"

dira_cut = "Lb_BPVDIRA>0.9998"

taudir_cut = "tau_END_VZ-Lb_END_VZ>-1"

protonID_cut = "proton_PROBNN_P>0.4 && proton_PROBNN_K<0.6 &&proton_PROBNN_PI<0.4"

kID_cut = "kaon_PROBNN_K>0.4 && kaon_PROBNN_PI<0.5"

sample_configs = {
    "taumu": {
        "sim_sample_name": "Lb_pKtaumu_3pi",
        "dat_sample_name": "turbo_rd",
        "sim_chain_name": "Hlt2RD_LbToPKTauMu_TauTo3Pi_OS",
        "dat_chain_name": "Hlt2RD_LbToPKTauMu_TauTo3Pi_SS",
        "lepton": "mu",
        "intpart": "L1520",
        "HLT1_cut": HLT1_cut_mu,
        "leptonID_cut": "mu_PROBNN_MU>0.2",
        "leptonID": "13",
        "dtfsuffix": "pKmudtf",
        "run2_sim_file": "mcLb_pKtaumu_3pi_2018_MagUp_mass1.root",
        "run2_dat_file": "Lb_pKtaumuSS_3pi_2018_MagUp_filtered_mass1.root",
        "run2_intpart": "Lc",
    }
}

run2_decayTree = "LbTuple1/DecayTree"

c1 = ROOT.TCanvas('c1', 'c1', 600, 500)
c1.SetTopMargin(0.06)
c1.SetRightMargin(0.05)
c1.SetBottomMargin(0.14)
c1.SetLeftMargin(0.16)
c1.SetTicks(1, 1)
c1.cd()


def plotInCanvas(histograms=None,
                 fileName=None,
                 unit=None,
                 legends=True,
                 var_units=None):
    h0 = histograms[0]
    if var_units is not None:
        h0.GetXaxis().SetTitle(var_units)
    binwidth = h0.GetXaxis().GetBinWidth(1)
    binwidth = round(binwidth, -int(math.floor(math.log10(abs(binwidth)) - 2)))
    h0.GetYaxis().SetTitle(f"Fraction of candidates / {binwidth} {unit}")

    for h in histograms:
        h.GetXaxis().SetTitleSize(0.055)
        h.GetXaxis().SetLabelSize(0.045)
        h.GetXaxis().SetTitleOffset(1.10)
        h.GetXaxis().SetNdivisions(505)
        h.GetYaxis().SetTitleSize(0.055)
        h.GetYaxis().SetLabelSize(0.045)
        h.GetYaxis().SetTitleOffset(1.30)
        h.GetYaxis().SetMaxDigits(3)

    ymax = max(h.GetMaximum() for h in histograms)
    h0.GetYaxis().SetRangeUser(0.0, 1.45 * ymax)

    h0.Draw("Hist E")
    for h in histograms[1:]:
        h.Draw("Hist E SAME")

    leg = None
    if legends:
        leg = ROOT.TLegend(0.62, 0.72, 0.92, 0.90)
        for h in histograms:
            leg.AddEntry(h.GetPtr(), h.GetTitle(), "l")
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextFont(42)
        leg.SetTextSize(0.050)
        leg.Draw()

    c1.RedrawAxis()
    c1.SaveAs(fileName)
    return leg


histograms = dict()

variables = {
    "proton_p": [
        "proton_P", "proton_P", 80, 0, 200000, "#it{p} (p) [MeV/#it{c}]",
        "[MeV/#it{c}]"
    ],
    # "L1520_Mass": ["L1520_M", "L1520_M", 80, 1300, 5700, "#it{m}(#it{#Lambda}(1520)) [MeV/#it{c}^{2}]", "[MeV/#it{c}^{2}]"],
    "L1520_pT": [
        "L1520_PT", "L1520_PT", 80, 0, 30000,
        "#it{p}_{T} (p + K^{-}) [MeV/#it{c}]", "[MeV/#it{c}]"
    ],
    "L1520_CHI2_DOF": [
        "L1520_ENDVERTEX_CHI2/1.0", "L1520_CHI2DOF", 80, -0.5, 8,
        "#it{#chi}^{2}_{vtx}/#it{n}_{dof} (p + K^{-})", ""
    ],
    # "L1520_FD_CHI2": ["L1520_FDCHI2_OWNPV", "L1520_BPVFDCHI2", 80, 0, 5000, "#it{#Lambda}(1520) FD #it{#chi}^{2}", ""],
    "mu_pT": [
        "mu_PT", "mu_PT", 80, 0, 15000, "#it{p}_{T} (#mu) [MeV/#it{c}]",
        "[MeV/#it{c}]"
    ],
    # "mu_p": ["mu_P","mu_P", 80, 1000, 30000, "#it{p}(#mu) [MeV/#it{c}]", "[MeV/#it{c}]"],
    "tau_Mass": [
        "tau_M", "tau_M", 80, 500, 2000, "#it{m} (#tau) [MeV/#it{c}^{2}]",
        "[MeV/#it{c}^{2}]"
    ],
    # "tau_FD_CHI2": ["tau_FDCHI2_OWNPV", "tau_BPVFDCHI2", 80, 0, 2000, "#tau FD #it{#chi}^{2}", ""],
    # "tau_CHI2_DOF": ["tau_ENDVERTEX_CHI2/3.0", "tau_CHI2DOF",  80, 0, 100, "#tau #it{#chi}^{2}_{vtx}/#it{n}_{dof}]", ""],
    # "Lb_Mass": ["Lb_M", "Lb_M", 80, 2000, 10000, "#it{m}(#it{#Lambda}_#it{b}) [MeV/#it{c}^{2}]", "[MeV/#it{c}^{2}]"],
    # "Lb_pT": ["Lb_PT","Lb_PT", 80, 500, 20000, "#it{p}_{T}(#it{#Lambda}_#it{b}) [MeV/#it{c}]", "[MeV/#it{c}]"],
    # "Lb_p": ["Lb_P","Lb_P", 80, 2000, 300000, "#it{p}(#it{#Lambda}_#it{b}) [MeV/#it{c}]", "[MeV/#it{c}]"],
    # "Lb_CHI2_DOF": ["Lb_ENDVERTEX_CHI2/5.0", "Lb_CHI2DOF",  80, 0, 80, "#it{#Lambda}_{#it{b}} #it{#chi}^{2}_{vtx}/#it{n}_{dof}]", ""],
    # "Lb_IP_CHI2": ["Lb_IPCHI2_OWNPV", "Lb_BPVIPCHI2", 80, 0, 40, "#it{#Lambda}_{#it{b}} IP #it{#chi}^{2}", ""],
    # "Lb_FD_CHI2": ["Lb_FDCHI2_OWNPV", "Lb_BPVFDCHI2", 80, 0, 5000, "#it{#Lambda}_{#it{b}} FD #it{#chi}^{2}", ""],
    "proton_pID_p": [
        "proton_PIDp", "proton_PID_P", 80, -10, 250,
        "#Delta log #it{L(p)} (p)", " "
    ],
    "proton_pID_diff": [
        "proton_PIDp-proton_PIDK", "proton_PID_P-proton_PID_K", 80, -50, 250,
        "#Delta log #it{L(p)} - #Delta log #it{L(K)} (p)", " "
    ],
    "kaon_pID_K": [
        "kaon_PIDK", "kaon_PID_K", 80, -20, 200,
        "#Delta log #it{L(K)} (K^{-})", " "
    ],
    "mu_pID_mu": [
        "mu_PIDmu", "mu_PID_MU", 80, -20, 40, "#Delta log #it{L(#mu)} (#mu)",
        " "
    ],
    "pi1_pID_K": [
        "pi1_PIDK", "pi1_PID_K", 80, -250, 200, "#Delta log #it{L(K)} (#pi 1)",
        " "
    ],
    "pi1_pID_p": [
        "pi1_PIDp", "pi1_PID_P", 80, -250, 200, "#Delta log #it{L(p)} (#pi 1)",
        " "
    ],
    "pi2_pID_K": [
        "pi2_PIDK", "pi2_PID_K", 80, -250, 200, "#Delta log #it{L(K)} (#pi 2)",
        " "
    ],
    "pi2_pID_p": [
        "pi2_PIDp", "pi2_PID_P", 80, -250, 200, "#Delta log #it{L(p)} (#pi 2)",
        " "
    ],
    "pi3_pID_K": [
        "pi3_PIDK", "pi3_PID_K", 80, -250, 200, "#Delta log #it{L(K)} (#pi 3)",
        " "
    ],
    "pi3_pID_p": [
        "pi3_PIDp", "pi3_PID_P", 80, -250, 200, "#Delta log #it{L(p)} (#pi 3)",
        " "
    ]
}

run2_line_color = ROOT.kBlue
run2_line_style = ROOT.kDashed

run3_line_color = ROOT.kRed
run3_line_style = ROOT.kSolid


def scale_histo(rdf_run2, rdf_run3, var, histo_run2, histo_run3, n_run2,
                n_run3):

    min2 = rdf_run2.Min(var).GetValue()
    max2 = rdf_run2.Max(var).GetValue()
    min3 = rdf_run3.Min(var).GetValue()
    max3 = rdf_run3.Max(var).GetValue()

    print(
        f"Range->  Run2: [{min2:.2f}, {max2:.2f}] and Run3: [{min3:.2f}, {max3:.2f}]"
    )

    tol = 0.01 * (max(max2, max3) - min(min2, min3))

    if max2 - max3 > tol:
        frac = rdf_run2.Filter(f"{var} < {max3}").Count().GetValue() / n_run2
        histo_run2.Scale(frac)

    elif max3 - max2 > tol:
        frac = rdf_run3.Filter(f"{var} < {max2}").Count().GetValue() / n_run3
        histo_run3.Scale(frac)

    if min2 - min3 > tol:
        frac = rdf_run3.Filter(f"{var} > {min2}").Count().GetValue() / n_run3
        histo_run3.Scale(frac)

    elif min3 - min2 > tol:
        frac = rdf_run2.Filter(f"{var} > {min3}").Count().GetValue() / n_run2
        histo_run2.Scale(frac)


chain_prefix = "DTF3"

for sample in ["taumu"]:

    lepton = sample_configs[sample]["lepton"]
    intpart = sample_configs[sample]["intpart"]
    leptonID = sample_configs[sample]["leptonID"]

    tm_cut = f"abs(pi1_TRUEID)==211 && abs(pi2_TRUEID)==211 && abs(pi3_TRUEID)==211 && abs(pi1_MC_MOTHER_ID)==15 && pi1_MC_MOTHER_KEY==pi2_MC_MOTHER_KEY && pi1_MC_MOTHER_KEY== pi3_MC_MOTHER_KEY && abs(proton_TRUEID)==2212 && abs(kaon_TRUEID)==321 && abs({lepton}_TRUEID)=={leptonID} && abs(proton_MC_MOTHER_ID)==5122 && proton_MC_MOTHER_KEY==kaon_MC_MOTHER_KEY && proton_MC_MOTHER_KEY=={lepton}_MC_MOTHER_KEY && proton_MC_MOTHER_KEY==pi1_MC_GD_MOTHER_KEY"

    sim_sample_name = sample_configs[sample]["sim_sample_name"]
    sim_chain_name = sample_configs[sample]["sim_chain_name"]
    run2_sim_file = sample_configs[sample]["run2_sim_file"]

    run3_chain = ROOT.TChain(f"{sim_chain_name}_{chain_prefix}/DecayTree")
    run3_friend_chain = ROOT.TChain("recomass")

    for sample_name, in_file in get_all_input_files(sim_sample_name):
        in_file_stem = Path(in_file).name.split('.')[0]
        print(f"Adding Run3 file {in_file}")
        run3_chain.Add(in_file)
        friend_file = friend_dir / f"recomass_{sim_chain_name}_{chain_prefix}_{in_file_stem}.root"
        print(f"Adding friend {friend_file}")
        run3_friend_chain.Add(str(friend_file))

    run3_chain.AddFriend(run3_friend_chain)

    run3_sample_rdf = ROOT.RDataFrame(run3_chain)
    run3_sample_rdf = run3_sample_rdf.Filter(tm_cut)

    run2_chain = ROOT.TChain(run2_decayTree)
    run2_chain.Add(str(run2_inputDir / run2_sim_file))

    run2_sample_rdf = ROOT.RDataFrame(run2_chain)
    run2_sample_rdf = run2_sample_rdf.Filter(tm_cut)

    histograms["sim"] = dict()

    plotfile_name = f"./plots/{sample}_sim_variables.pdf"

    c1.SaveAs(f"{plotfile_name}[")

    for var, cons in variables.items():

        name = "sim"

        histograms[name][var] = dict()

        run2_sample_rdf = run2_sample_rdf.Define(var, cons[0])
        run3_sample_rdf = run3_sample_rdf.Define(var, cons[1])

        var_cut = f"{var} > {cons[3]} && {var} < {cons[4]}"

        run2_sample_rdf = run2_sample_rdf.Filter(var_cut)
        run3_sample_rdf = run3_sample_rdf.Filter(var_cut)

        histograms[name][var]["histo_run2"] = run2_sample_rdf.Histo1D(
            ROOT.RDF.TH1DModel(f"h_{var}_run2", "", cons[2], cons[3], cons[4]),
            var)
        histograms[name][var]["histo_run2"].SetLineColor(run2_line_color)
        histograms[name][var]["histo_run2"].SetLineStyle(run2_line_style)
        histograms[name][var]["histo_run2"].SetLineWidth(2)
        histograms[name][var]["histo_run2"].SetMarkerSize(0)
        histograms[name][var]["histo_run2"].SetTitle("Run 2")

        run2_hist_entries = histograms[name][var]["histo_run2"].GetEntries()
        print("Run2 entries: ", run2_hist_entries)
        if run2_hist_entries > 0:
            histograms[name][var]["histo_run2"].Scale(1.0 / run2_hist_entries)
        else:
            print(
                f"  WARNING: Run2 histogram for '{var}' has 0 entries -- skipping scale."
            )

        histograms[name][var]["histo_run3"] = run3_sample_rdf.Histo1D(
            ROOT.RDF.TH1DModel(f"h_{var}_run3", "", cons[2], cons[3], cons[4]),
            var)
        histograms[name][var]["histo_run3"].SetLineColor(run3_line_color)
        histograms[name][var]["histo_run3"].SetLineStyle(run3_line_style)
        histograms[name][var]["histo_run3"].SetLineWidth(2)
        histograms[name][var]["histo_run3"].SetMarkerSize(0)
        histograms[name][var]["histo_run3"].SetTitle("Run 3")

        run3_hist_entries = histograms[name][var]["histo_run3"].GetEntries()
        print("Run3 entries: ", run3_hist_entries)
        if run3_hist_entries > 0:
            histograms[name][var]["histo_run3"].Scale(1.0 / run3_hist_entries)
        else:
            print(
                f"  WARNING: Run3 histogram for '{var}' has 0 entries -- skipping scale."
            )

        n_run2 = run2_sample_rdf.Count().GetValue()
        n_run3 = run3_sample_rdf.Count().GetValue()

        if n_run2 > 0 and n_run3 > 0:
            scale_histo(run2_sample_rdf, run3_sample_rdf, var,
                        histograms[name][var]["histo_run2"],
                        histograms[name][var]["histo_run3"], n_run2, n_run3)
        else:
            print(
                f"  WARNING: '{var}' has empty sample (run2={n_run2}, run3={n_run3}) -- skipping overlap scale."
            )

        maxima = list()
        maxima.append(histograms[name][var]["histo_run2"].GetMaximum())
        maxima.append(histograms[name][var]["histo_run3"].GetMaximum())
        ''' 
        for histogram in [histograms[name][var][f"histo_{run}"] for run in ["run2", "run3"]]:
            histogram.GetYaxis().SetRangeUser(0, 1.05 * max(maxima))
        '''
        c1.cd()

        plotInCanvas([
            histograms[name][var]["histo_run2"],
            histograms[name][var]["histo_run3"]
        ], plotfile_name, cons[6], True, cons[5])

    c1.SaveAs(f"{plotfile_name}]")

c1.Close()
