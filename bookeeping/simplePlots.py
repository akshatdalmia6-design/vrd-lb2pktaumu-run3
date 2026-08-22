import ROOT
import sys
import math
import numpy as np
from run3_samples_v1r5791 import samples
from pathlib import Path
import glob

sys.path.append("../utils")
from lhcbPlotStyle import setLHCbPlotStyle

setLHCbPlotStyle()

friend_dir = Path(" ")

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
    }
}

decayTree = "LbTuple1/DecayTree"

x_min = 4500
x_max = 8000


def get_all_input_files(sample_name):
    """Generator to yield sample context and file path for a specific sample name."""
    # Check if the sample name exists in your samples configuration
    if sample_name not in samples:
        print(
            f"Warning: Sample '{sample_name}' not found in samples dictionary."
        )
        return

    # Loop only over the dictionary keys (block paths) for this specific sample
    for block, block_path in samples[sample_name].items():

        # if "turbo" in sample_name and block not in ["blk8"]:
        #     continue

        for file_path in glob.glob(block_path):
            yield sample_name, file_path


def plotInCanvas(histograms=None, fileName=None, legends=True):

    histograms[0].GetXaxis().SetTitle(
        "#it{#Lambda}_{#it{b}} mass [MeV/#it{c}^{2}]")
    binwidth = histograms[0].GetXaxis().GetBinWidth(1)
    binwidth = round(binwidth, -int(math.floor(math.log10(abs(binwidth)) - 2)))
    histograms[0].GetYaxis().SetTitle("Candidates per {0} {1}".format(
        binwidth, "MeV/#it{c}^{2}"))

    histograms[0].Draw("Hist E")

    for histogram in histograms[1:]:
        histogram.Draw("Hist E SAME")

    if legends:
        leg = c1.BuildLegend(0.62, 0.60, 0.90, 0.90, '', '')
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextFont(42)
        leg.SetTextSize(0.055)

    c1.SaveAs(fileName)


def rms90(histo=None, imean=None):

    # 1. Determine the central 90%
    nbins = histo.GetNbinsX()

    # Bin of the mean (in this case maximum)
    if imean is None:
        imean = histo.GetMaximumBin()

    entries = 0.9 * histo.GetEntries()  # Target entries

    sumw = histo.GetBinContent(imean)
    sumwx = sumw * histo.GetBinCenter(imean)
    sumwx2 = sumw * pow(histo.GetBinCenter(imean),
                        2)  # Symmetrically expand outward

    for i in range(1, nbins):

        if (sumw >= entries):
            break

        # Check bin to the left
        if (imean - i > 0):
            w = histo.GetBinContent(imean - i)
            x = histo.GetBinCenter(imean - i)
            sumw += w
            sumwx += w * x
            sumwx2 += w * x * x

        # Check bin to the right
        if (imean + i <= nbins):
            w = histo.GetBinContent(imean + i)
            x = histo.GetBinCenter(imean + i)
            sumw += w
            sumwx += w * x
            sumwx2 += w * x * x

    # 2. Calculate the RMS for the central 90%
    mean90 = sumwx / sumw
    # Variance (RMS^2)
    rms2 = abs(sumwx2 / sumw - mean90 * mean90)
    result = pow(rms2, 0.5)

    print(f"RMS of central 90% = {result}, Total RMS = {histo.GetRMS()}")

    return result


histograms = dict()

# Comparison plots
if __name__ == '__main__':

    c1 = ROOT.TCanvas('c1', 'c1', 600, 400)
    c1.SetTopMargin(0.07)
    c1.SetLeftMargin(0.17)
    c1.cd()

    for sample in ["taumu"]:

        lepton = sample_configs[sample]["lepton"]
        intpart = sample_configs[sample]["intpart"]

        HLT1_cut = sample_configs[sample]["HLT1_cut"]
        leptonID_cut = sample_configs[sample]["leptonID_cut"]
        leptonID = sample_configs[sample]["leptonID"]

        dtfsuffix = sample_configs[sample]["dtfsuffix"]

        masses = {
            "Lb_M_corr":
            ["DTF3", "Lb_M_corr > 0", ROOT.kGray + 3, ROOT.kDashed],
            "Lb_M_reco":
            ["DTF3", "Delta_sol > -1.e13", ROOT.kBlue + 1, ROOT.kSolid],
            "Lb_M_DTF0":
            ["DTF0", f"Lb_DTF_CHI2DOF < 50", ROOT.kMagenta + 2, ROOT.kDashed],
            "Lb_M_DTF1": [
                "DTF1", f"Lb_DTF_CHI2DOF < 50", ROOT.kOrange + 7,
                ROOT.kDashDotted
            ],
            "Lb_M_DTF2": [
                "DTF2", f"Lb_DTF_CHI2DOF < 50", ROOT.kRed + 1, ROOT.kDotted
            ],
            "Lb_M_DTF3": ["DTF3", f"Lb_DTF_CHI2DOF < 50", ROOT.kGreen + 3, 5],
            "Lb_M_hybr": [
                "DTF3", f"(Delta_sol > 0 || Lb_DTF_CHI2DOF < 50)",
                ROOT.kBlack + 3, 4
            ],
        }

        tm_cut = f"abs(pi1_TRUEID)==211 && abs(pi2_TRUEID)==211 && abs(pi3_TRUEID)==211 && abs(pi1_MC_MOTHER_ID)==15 && pi1_MC_MOTHER_KEY==pi2_MC_MOTHER_KEY && pi1_MC_MOTHER_KEY== pi3_MC_MOTHER_KEY && abs(proton_TRUEID)==2212 && abs(kaon_TRUEID)==321 && abs({lepton}_TRUEID)=={leptonID} && abs(proton_MC_MOTHER_ID)==5122 && proton_MC_MOTHER_KEY==kaon_MC_MOTHER_KEY && proton_MC_MOTHER_KEY=={lepton}_MC_MOTHER_KEY && proton_MC_MOTHER_KEY==pi1_MC_GD_MOTHER_KEY"

        # Initial

        Lc_cut = f"{intpart}_M > 2300"

        displ_cut = f"proton_BPVIPCHI2 > 25 && kaon_BPVIPCHI2 > 25 && {lepton}_BPVIPCHI2 > 25 && pi1_BPVIPCHI2 > 25 && pi2_BPVIPCHI2 > 25 && pi3_BPVIPCHI2 > 25"

        FD_cut = f"{intpart}_BPVFDCHI2 > 25"

        initial_cuts = " && ".join([Lc_cut, displ_cut, HLT1_cut, FD_cut])

        # Tighter

        tauvtx_cut = f"tau_END_VZ-{intpart}_END_VZ > -1"

        fullnopid_cut = " && ".join(
            [belowD_cut, dira_cut, taudir_cut, tauvtx_cut, initial_cuts])

        if lepton == "electron":
            fullnopid_cut += " && Lb_PT>4500 && tau_PT>1500"

        full_cut = " && ".join(
            [protonID_cut, leptonID_cut, kID_cut, fullnopid_cut])

        print(full_cut)

        sim_sample_name = sample_configs[sample]["sim_sample_name"]
        sim_chain_name = sample_configs[sample]["sim_chain_name"]
        dat_sample_name = sample_configs[sample]["dat_sample_name"]
        dat_chain_name = sample_configs[sample]["dat_chain_name"]

        for sample_name, chain_name, name, cut in [
            (sim_sample_name, sim_chain_name, "sim", tm_cut),  # tm_cut),
                #                         (dat_chain, "dat", fullnopid_cut)]:
            (dat_sample_name, dat_chain_name, "dat", full_cut)
        ]:

            histograms[name] = dict()

            # 2. Optimised Hybrid Logic
            # We compare squared momenta to avoid two sqrt() calls
            hybrid_logic = f"""
            if (Delta_sol > 0) {{
                double p2_max = tau_PX_reco_max*tau_PX_reco_max + tau_PY_reco_max*tau_PY_reco_max + tau_PZ_reco_max*tau_PZ_reco_max;
                double p2_min = tau_PX_reco_min*tau_PX_reco_min + tau_PY_reco_min*tau_PY_reco_min + tau_PZ_reco_min*tau_PZ_reco_min;
                return (p2_max> p2_min) ? Lb_M_reco_max : Lb_M_reco_min;
            }} else {{
                return (Lb_DTF_CHI2DOF < 50) ? Lb_DTF_M : -1.0;
            }}
            """

            maxima = list()

            for mass in masses.keys():

                print(f"\n{mass}")

                mass_var = "Lb_DTF_M" if "DTF" in mass else mass

                chain_prefix = masses[mass][0]

                chain = ROOT.TChain(f"{chain_name}_{chain_prefix}/DecayTree")
                friend_chain = ROOT.TChain("recomass")

                for sample_name, in_file in get_all_input_files(sample_name):
                    in_file_stem = Path(in_file).name.split('.')[0]
                    print(f"Adding {in_file}")

                    chain.Add(in_file)

                    if chain_prefix == "DTF3":

                        friend_file = friend_dir / f"recomass_{chain_name}_{chain_prefix}_{in_file_stem}.root"
                        print(f"Adding friend {friend_file}")
                        friend_chain.Add(str(friend_file))

                if chain_prefix == "DTF3":
                    chain.AddFriend(friend_chain)

                sample_rdf = ROOT.RDataFrame(chain)

                if mass == "Lb_M_hybr":
                    sample_rdf = sample_rdf.Define("Lb_M_hybr", hybrid_logic)

                entries_0 = sample_rdf.Count().GetValue()

                if cut is not None:
                    sample_rdf = sample_rdf.Filter(cut)

                entries_1 = sample_rdf.Count().GetValue()

                print(
                    f"Entries in {name}: {entries_0}, {entries_1}, ratio = {float(entries_1/entries_0):.4f}"
                )

                histograms[name][mass] = dict()

                convergence_cut = masses[mass][1]
                line_color = masses[mass][2]
                line_style = masses[mass][3]

                mass_cut = f"{mass_var} > {x_min} && {mass_var} < {x_max}"

                sample_entries_inrange = sample_rdf.Filter(
                    mass_cut).Count().GetValue()

                histograms[name][mass]["histo_fine"] = sample_rdf.Filter(
                    f"{mass_cut} && {convergence_cut}").Histo1D(
                        ROOT.RDF.TH1DModel(f"h_{mass}", "", 3500, x_min,
                                           x_max), mass_var)
                histograms[name][mass]["histo"] = sample_rdf.Filter(
                    f"{mass_cut} && {convergence_cut}").Histo1D(
                        ROOT.RDF.TH1DModel(f"h_{mass}", "", 100, x_min, x_max),
                        mass_var)
                histograms[name][mass]["histo"].SetLineColor(line_color)
                histograms[name][mass]["histo"].SetLineStyle(line_style)
                histograms[name][mass]["histo"].SetLineWidth(2)
                histograms[name][mass]["histo"].SetMarkerSize(0)
                histograms[name][mass]["histo"].SetTitle(mass)
                maxima.append(histograms[name][mass]["histo"].GetMaximum())

                sample_entries_converged = histograms[name][mass][
                    "histo"].GetEntries()
                print("Entries after selection: ", sample_entries_converged)

                histograms[name][mass][
                    "convrate"] = sample_entries_converged / sample_entries_inrange
                histograms[name][mass]["convrate_err"] = np.sqrt(
                    histograms[name][mass]["convrate"] *
                    (1 - histograms[name][mass]["convrate"]) /
                    sample_entries_inrange)

                histograms[name][mass][
                    "effcy"] = sample_entries_converged / entries_1
                histograms[name][mass]["effcy_err"] = np.sqrt(
                    histograms[name][mass]["effcy"] *
                    (1 - histograms[name][mass]["effcy"]) / entries_1)

                histograms[name][mass]["bias"] = histograms[name][mass][
                    "histo_fine"].GetBinCenter(
                        histograms[name][mass]["histo_fine"].
                        GetMaximumBin()) - 5619.5

                histograms[name][mass]["rms90"] = rms90(
                    histograms[name][mass]["histo_fine"])

                print(histograms[name][mass]["histo"].GetRMS())

            c1.cd()

            plotfile_name = f"./plots/{sample}_{name}_masses.pdf"

            c1.SaveAs(f"{plotfile_name}[")

            # All y axes in the same range
            for histogram in [
                    histograms[name][mass]["histo"] for mass in masses.keys()
            ]:
                histogram.GetYaxis().SetRangeUser(0, 1.05 * max(maxima))

            plotInCanvas([
                histograms[name]["Lb_M_reco"]["histo"],
                histograms[name]["Lb_M_corr"]["histo"]
            ], plotfile_name)

            plotInCanvas([
                histograms[name]["Lb_M_hybr"]["histo"],
                histograms[name]["Lb_M_reco"]["histo"]
            ], plotfile_name)

            plotInCanvas([histograms[name]["Lb_M_hybr"]["histo"]] + [
                histograms[name][f"Lb_M_DTF{i}"]["histo"] for i in range(0, 4)
            ], plotfile_name)

            c1.SaveAs(f"{plotfile_name}]")

        header = r'''
        \begin{table}[]
            \centering
            \begin{tabular}{c c c r c c}
            \toprule
                Method & Convergence rate & Efficiency & Bias & RMS90 & Bkg. in fit range \\\midrule
        '''

        bottom = r'''
            \\\bottomrule
            \end{tabular}
        '''

        print(header)

        for mass in masses.keys():

            method = rf'\texttt{{{mass}}}'
            convrate = histograms["sim"][mass]["convrate"]
            convrate_err = histograms["sim"][mass]["convrate_err"]

            effcy = histograms["sim"][mass]["effcy"]
            effcy_err = histograms["sim"][mass]["effcy_err"]

            bias = histograms["sim"][mass]["bias"]
            rms90 = histograms["sim"][mass]["rms90"]

            bkgfrac = histograms["dat"][mass]["effcy"]
            bkgfrac_err = histograms["dat"][mass]["effcy_err"]

            print(
                rf"{method} & {convrate:.4f} & {effcy:.2f} & {bias} & {rms90:.0f} & {bkgfrac:.2f} \\"
            )

        print(bottom)
