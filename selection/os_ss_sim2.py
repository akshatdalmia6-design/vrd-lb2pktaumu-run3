import ROOT
import sys
import math
from pathlib import Path

sys.path.append("../utils")
from lhcbPlotStyle import setLHCbPlotStyle
setLHCbPlotStyle()

ROOT.EnableImplicitMT()

skim_dir = Path(
    " ")

LB_PDG_M = 5619.5
NSIGMA = 3.0
x_min = 4500
x_max = 8000

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
        leg = ROOT.TLegend(0.50, 0.72, 0.80, 0.90)
        for h in histograms:
            leg.AddEntry(h, h.GetTitle(), "l")
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextFont(42)
        leg.SetTextSize(0.050)
        leg.Draw()

    c1.RedrawAxis()
    c1.SaveAs(fileName)
    return leg


def rms90(histo=None, imean=None):

    nbins = histo.GetNbinsX()
    if imean is None:
        imean = histo.GetMaximumBin()

    entries = 0.9 * histo.GetEntries()

    sumw = histo.GetBinContent(imean)
    sumwx = sumw * histo.GetBinCenter(imean)
    sumwx2 = sumw * pow(histo.GetBinCenter(imean), 2)

    for i in range(1, nbins):
        if (sumw >= entries):
            break
        if (imean - i > 0):
            w = histo.GetBinContent(imean - i)
            x = histo.GetBinCenter(imean - i)
            sumw += w
            sumwx += w * x
            sumwx2 += w * x * x
        if (imean + i <= nbins):
            w = histo.GetBinContent(imean + i)
            x = histo.GetBinCenter(imean + i)
            sumw += w
            sumwx += w * x
            sumwx2 += w * x * x

    mean90 = sumwx / sumw
    rms2 = abs(sumwx2 / sumw - mean90 * mean90)
    result = pow(rms2, 0.5)
    print(f"RMS of central 90% = {result}, Total RMS = {histo.GetRMS()}")
    return result


def separation(h_sig, h_bkg):
    """TMVA-style separation <S^2> = 1/2 sum (s-b)^2/(s+b) for unit-area histograms."""
    sep = 0.0
    for i in range(1, h_sig.GetNbinsX() + 1):
        s = h_sig.GetBinContent(i)
        b = h_bkg.GetBinContent(i)
        if s + b > 0:
            sep += 0.5 * (s - b) * (s - b) / (s + b)
    return sep


# Candidate BDT input variables.
# Format: [expression, nbins, xmin, xmax, axis title, unit]
variables = {
    # --- proton ---
    "proton_p":
    ["proton_P", 80, 0, 100000, "#it{p}(p) [MeV/#it{c}]", "[MeV/#it{c}]"],
    "proton_pT": [
        "proton_PT", 80, -100, 14000, "#it{p}_{T} (p) [MeV/#it{c}]",
        "[MeV/#it{c}]"
    ],
    "proton_pID_P": ["proton_PID_P", 80, 0, 500, "PID(p) of p", ""],
    "proton_IP_CHI2":
    ["proton_BPVIPCHI2", 80, -50, 1200, "#{IP} #it{#chi}^{2} (p)", ""],

    # --- kaon ---
    "kaon_pT": [
        "kaon_PT", 80, 0, 20000, "#it{p}_{T}(#it{K}) [MeV/#it{c}]",
        "[MeV/#it{c}]"
    ],
    "kaon_p":
    ["kaon_P", 80, 0, 100000, "#it{p}(#it{K}) [MeV/#it{c}]", "[MeV/#it{c}]"],
    "kaon_pID_K": ["kaon_PID_K", 80, 0, 500, "PID(#it{K}) of #it{K}", ""],
    "kaon_IP_CHI2": [
        "kaon_BPVIPCHI2", 80, -100, 1500, "#it{IP} #it{#chi}^{2} (K^{-})", ""
    ],

    # --- L(1520) ---
    "L1520_Mass": [
        "L1520_M", 80, 1300, 5700,
        "#it{m}(#it{#Lambda}(1520)) [MeV/#it{c}^{2}]", "[MeV/#it{c}^{2}]"
    ],
    "L1520_pT": [
        "L1520_PT", 80, 0, 18000, "#it{p}_{T}(p + K^{-}) [MeV/#it{c}]",
        "[MeV/#it{c}]"
    ],
    "L1520_CHI2_DOF": [
        "L1520_CHI2DOF", 80, 0, 20,
        "#it{#Lambda}(1520) #it{#chi}^{2}_{vtx}/#it{n}_{dof}", ""
    ],
    "L1520_FD_CHI2": [
        "L1520_BPVFDCHI2", 80, -50, 1800, "#it{FD} #it{#chi}^{2} (p + K^{-})",
        ""
    ],
    "L1520_MAXS_DOCA": [
        "L1520_MAXSDOCACHI2", 80, 0, 50,
        "#it{#Lambda}(1520) max SDOCA #it{#chi}^{2}", ""
    ],
    "L1520_p": [
        "L1520_P", 80, 0, 100000, "#it{p}(#it{#Lambda}(1520)) [MeV/#it{c}]",
        "[MeV/#it{c}]"
    ],
    "L1520_DIRA": [
        "L1520_BPVDIRA", 80, 0.99, 1.0, "#it{#Lambda}(1520) DIRA", ""
    ],
    "L1520_IP_CHI2": [
        "L1520_BPVIPCHI2", 80, 0, 2000, "#Lambda(1520) IP #it{#chi}^{2}", ""
    ],

    # --- muon ---
    "mu_pT": [
        "mu_PT", 80, 0, 20000, "#it{p}_{T}(#mu) [MeV/#it{c}]", "[MeV/#it{c}]"
    ],
    "mu_p": [
        "mu_P", 80, 0, 100000, "#it{p}(#mu) [MeV/#it{c}]", "[MeV/#it{c}]"
    ],
    "mu_pID_MU": ["mu_PID_MU", 80, 0, 500, "PID(#mu) of #mu", ""],
    "mu_IP_CHI2": ["mu_BPVIPCHI2", 80, 0, 2000, "#mu IP #it{#chi}^{2}", ""],
    "mu_CHI2_DOF":
    ["mu_CHI2DOF", 80, 0, 20, "#it{#mu} #it{#chi}^{2}_{vtx}/#it{n}_{dof}", ""],

    # --- tau (3pi) ---
    "tau_Mass": [
        "tau_M", 80, 500, 2000, "#it{m}(3#pi) [MeV/#it{c}^{2}]",
        "[MeV/#it{c}^{2}]"
    ],
    "tau_pT": [
        "tau_PT", 80, 0, 20000, "#it{p}_{T}(3#pi) [MeV/#it{c}]", "[MeV/#it{c}]"
    ],
    "tau_MAXS_DOCA": [
        "tau_MAXSDOCACHI2", 80, 0, 50, "3#pi max SDOCA #it{#chi}^{2}", ""
    ],
    "tau_FD_CHI2": ["tau_BPVFDCHI2", 80, 0, 2000, "#tau FD #it{#chi}^{2}", ""],
    "tau_DIRA": ["tau_BPVDIRA", 80, 0.99, 1.0, "#tau DIRA", ""],
    "tau_CHI2_DOF": [
        "tau_CHI2DOF", 80, 0, 20, "#tau #it{#chi}^{2}_{vtx}/#it{n}_{dof}", ""
    ],
    "tau_DOCA_CHI2_13": [
        "tau_DOCACHI2_13", 80, 0, 50, "3#pi DOCA #it{#chi}^{2} (1,3)", ""
    ],
    "tau_DOCA_CHI2_23": [
        "tau_DOCACHI2_23", 80, 0, 50, "3#pi DOCA #it{#chi}^{2} (2,3)", ""
    ],
    "tau_MAX_pT": [
        "tau_MAX_PT", 80, 0, 15000, "max #it{p}_{T} of 3#pi [MeV/#it{c}]",
        "[MeV/#it{c}]"
    ],
    "tau_MAX_IPCHI2": [
        "tau_MAX_BPVIPCHI2", 80, 0, 5000, "max IP #it{#chi}^{2} of 3#pi", ""
    ],
    "tau_p": [
        "tau_P", 80, 0, 80000, "#it{p}(3#pi) [MeV/#it{c}]", "[MeV/#it{c}]"
    ],

    # --- pi1 ---
    "pi1_pT": [
        "pi1_PT", 80, 0, 15000, "#it{p}_{T}(#pi_{1}) [MeV/#it{c}]",
        "[MeV/#it{c}]"
    ],
    "pi1_p": [
        "pi1_P", 80, 0, 80000, "#it{p}(#pi_{1}) [MeV/#it{c}]", "[MeV/#it{c}]"
    ],
    "pi1_pID_K": ["pi1_PID_K", 80, 0, 5, "PID(#it{K}) of #pi_{1}", ""],
    "pi1_pID_E": ["pi1_PID_E", 80, 0, 5, "PID(#it{e}) of #pi_{1}", ""],
    "pi1_pID_P": ["pi1_PID_P", 80, 0, 10, "PID(p) of #pi_{1}", ""],
    "pi1_pID_MU": ["pi1_PID_MU", 80, 0, 10, "PID(#mu) of #pi_{1}", ""],
    "pi1_IP_CHI2": [
        "pi1_BPVIPCHI2", 80, 0, 2000, "#pi_{1} IP #it{#chi}^{2}", ""
    ],
    "pi1_CHI2_DOF":
    ["pi1_CHI2DOF", 80, 0, 20, "#pi_{1} #it{#chi}^{2}_{vtx}/#it{n}_{dof}", ""],

    # --- pi2 ---
    "pi2_pT": [
        "pi2_PT", 80, 0, 15000, "#it{p}_{T}(#pi_{2}) [MeV/#it{c}]",
        "[MeV/#it{c}]"
    ],
    "pi2_p": [
        "pi2_P", 80, 0, 80000, "#it{p}(#pi_{2}) [MeV/#it{c}]", "[MeV/#it{c}]"
    ],
    "pi2_pID_K": ["pi2_PID_K", 80, 0, 5, "PID(#it{K}) of #pi_{2}", ""],
    "pi2_pID_E": ["pi2_PID_E", 80, 0, 5, "PID(#it{e}) of #pi_{2}", ""],
    "pi2_pID_P": ["pi2_PID_P", 80, 0, 10, "PID(p) of #pi_{2}", ""],
    "pi2_pID_MU": ["pi2_PID_MU", 80, 0, 10, "PID(#mu) of #pi_{2}", ""],
    "pi2_IP_CHI2": [
        "pi2_BPVIPCHI2", 80, 0, 2000, "#pi_{2} IP #it{#chi}^{2}", ""
    ],
    "pi2_CHI2_DOF": [
        "pi2_CHI2DOF", 80, 0, 20, "#pi_{2} #it{#chi}^{2}_{vtx}/#it{n}_{dof}",
        ""
    ],

    # --- pi3 ---
    "pi3_pT": [
        "pi3_PT", 80, 0, 15000, "#it{p}_{T}(#pi_{3}) [MeV/#it{c}]",
        "[MeV/#it{c}]"
    ],
    "pi3_p": [
        "pi3_P", 80, 0, 80000, "#it{p}(#pi_{3}) [MeV/#it{c}]", "[MeV/#it{c}]"
    ],
    "pi3_pID_K": ["pi3_PID_K", 80, 0, 5, "PID(#it{K}) of #pi_{3}", ""],
    "pi3_pID_E": ["pi3_PID_E", 80, 0, 5, "PID(#it{e}) of #pi_{3}", ""],
    "pi3_pID_P": ["pi3_PID_P", 80, 0, 10, "PID(p) of #pi_{3}", ""],
    "pi3_pID_MU": ["pi3_PID_MU", 80, 0, 10, "PID(#mu) of #pi_{3}", ""],
    "pi3_IP_CHI2": [
        "pi3_BPVIPCHI2", 80, 0, 2000, "#pi_{3} IP #it{#chi}^{2}", ""
    ],
    "pi3_CHI2_DOF": [
        "pi3_CHI2DOF", 80, 0, 20, "#pi_{3} #it{#chi}^{2}_{vtx}/#it{n}_{dof}",
        ""
    ],

    # --- Jpsi  ---
    "Jpsi_Mass": [
        "Jpsi_M", 80, 100, 5100, "#it{m}(#it{J}/#psi cand.) [MeV/#it{c}^{2}]",
        "[MeV/#it{c}^{2}]"
    ],
    "Jpsi_CHI2_DOF": [
        "Jpsi_CHI2DOF", 80, 0, 20,
        "#it{J}/#psi #it{#chi}^{2}_{vtx}/#it{n}_{dof}", ""
    ],
    "Jpsi_pT": [
        "Jpsi_PT", 80, 0, 16000, "#it{p}_{T} (#tau + #mu) [MeV/#it{c}]",
        "[MeV/#it{c}]"
    ],
    "Jpsi_p": [
        "Jpsi_P", 80, 0, 80000, "#it{p} #it{J}/#psi [MeV/#it{c}]",
        "[MeV/#it{c}]"
    ],
    "Jpsi_IP_CHI2": [
        "Jpsi_BPVIPCHI2", 80, 0, 40, "#it{J}/#psi IP #it{#chi}^{2}", ""
    ],
    "Jpsi_FD_CHI2": [
        "Jpsi_BPVFDCHI2", 80, 0, 5000, "#it{J}/#psi FD #it{#chi}^{2}", ""
    ],
    "Jpsi_DIRA": ["Jpsi_BPVDIRA", 80, 0.995, 1.0, "#it{J}/#psi DIRA", ""],
    "Jpsi_MAXS_DOCA": [
        "Jpsi_MAXSDOCACHI2", 80, 0, 50, "#it{J}/#psi max SDOCA #it{#chi}^{2}",
        ""
    ],

    # --- Lb ---
    "Lb_Mass": [
        "Lb_M", 80, 2000, 10000,
        "#it{m}(#it{#Lambda}_{#it{b}}) [MeV/#it{c}^{2}]", "[MeV/#it{c}^{2}]"
    ],
    "Lb_CHI2_DOF": [
        "Lb_CHI2DOF", 80, -5, 30,
        "#it{#chi}^{2}_{vtx}/#it{n}_{dof} (#Lambda_{b})", ""
    ],
    "Lb_p": [
        "Lb_P", 80, 0, 300000, "#it{p}(#it{#Lambda}_{#it{b}}) [MeV/#it{c}]",
        "[MeV/#it{c}]"
    ],
    "Lb_pT": [
        "Lb_PT", 80, 0, 30000, "#it{p}_{T}(#Lambda_{b}) [MeV/#it{c}]",
        "[MeV/#it{c}]"
    ],
    "Lb_IP_CHI2": [
        "Lb_BPVIPCHI2", 80, 0, 40, "#it{#Lambda}_{#it{b}} IP #it{#chi}^{2}", ""
    ],
    "Lb_FD_CHI2": [
        "Lb_BPVFDCHI2", 80, 0, 4500, "#it{FD} #it{#chi}^{2} (#Lambda_{b})", ""
    ],
    "Lb_DIRA": ["Lb_BPVDIRA", 80, 0.999, 1.001, "#it{DIRA} (#Lambda_{b})", ""],
    "Lb_MAXS_DOCA": [
        "Lb_MAXSDOCACHI2", 80, 0, 50,
        "#it{#Lambda}_{#it{b}} max SDOCA #it{#chi}^{2}", ""
    ],
}

sig_line_color = ROOT.kRed
sig_line_style = ROOT.kSolid
bkg_line_color = ROOT.kBlue
bkg_line_style = ROOT.kDashed

histograms = dict()

if __name__ == '__main__':

    for sample in ["taumu"]:

        # Skims are already preselected + DTF-converged; read them directly.
        sim_rdf_base = ROOT.RDataFrame("DecayTree",
                                       str(skim_dir / f"{sample}_sim_OS.root"))
        ss_rdf_base = ROOT.RDataFrame("DecayTree",
                                      str(skim_dir / f"{sample}_dat_SS.root"))
        os_rdf_base = ROOT.RDataFrame("DecayTree",
                                      str(skim_dir / f"{sample}_dat_OS.root"))

        # ----- peak + RMS90 from the simulation DTF mass -----
        histo_fine = sim_rdf_base.Filter(
            f"Lb_DTF_M > {x_min} && Lb_DTF_M < {x_max}").Histo1D(
                ROOT.RDF.TH1DModel("h_Lb_DTF_M_fine", "", 3500, x_min, x_max),
                "Lb_DTF_M")

        m_peak = histo_fine.GetBinCenter(histo_fine.GetMaximumBin())
        sigma_peak = rms90(histo_fine.GetValue())
        print(
            f"\nDTF mass peak at {m_peak:.1f} MeV (PDG: {LB_PDG_M}), RMS90 = {sigma_peak:.1f} MeV"
        )

        window_lo = m_peak - NSIGMA * sigma_peak
        window_hi = m_peak + NSIGMA * sigma_peak
        print(f"Signal window: [{window_lo:.1f}, {window_hi:.1f}] MeV")

        signal_window_cut = f"Lb_DTF_M > {window_lo} && Lb_DTF_M < {window_hi}"
        upper_sideband_cut = f"Lb_DTF_M > {window_hi} && Lb_DTF_M < {x_max}"
        full_range_cut = f"Lb_DTF_M > {x_min} && Lb_DTF_M < {x_max}"

        sim_rdf = sim_rdf_base.Filter(signal_window_cut)

        bkg_proxies = {
            "ssdata": {
                "rdf": ss_rdf_base.Filter(full_range_cut),
                "legend": "Run 3 SS data (bkg)"
            },
            "osdata_sideband": {
                "rdf": os_rdf_base.Filter(upper_sideband_cut),
                "legend": "Run 3 OS sideband (bkg)"
            },
        }

        # Define every variable once on each sample.
        sig_def = sim_rdf
        bkg_defs = {
            bkg_name: proxy["rdf"]
            for bkg_name, proxy in bkg_proxies.items()
        }
        for var, cons in variables.items():
            sig_def = sig_def.Define(var, cons[0])
            for bkg_name in bkg_defs:
                bkg_defs[bkg_name] = bkg_defs[bkg_name].Define(var, cons[0])
        '''

        all_defs = {"sig": sig_def, **bkg_defs}

        minmax = {s: {var: {"min": d.Min(var), "max": d.Max(var)} for var in variables}
                  for s, d in all_defs.items()}

        results = [r for s in minmax.values() for v in s.values() for r in v.values()]
        print(f"\nFinding min/max for {len(variables)} variables...")
        ROOT.RDF.RunGraphs(results)

        # Common range per variable = union of all samples' [min, max].
        ranges = {}
        for var in variables:
            xlo = min(minmax[s][var]["min"].GetValue() for s in all_defs)
            xhi = max(minmax[s][var]["max"].GetValue() for s in all_defs)
            if not (xhi > xlo):            # constant variable
                xlo, xhi = xlo - 0.5, xhi + 0.5
            ranges[var] = (xlo, xhi)
            print(f"  {var:<20} range=[{xlo:.3g}, {xhi:.3g}]")
        '''

        # ==============================================================
        # PHASE 1: book everything (no data read yet)
        # ==============================================================
        booked = {bkg_name: dict() for bkg_name in bkg_proxies}

        for bkg_name, proxy in bkg_proxies.items():
            for var, cons in variables.items():
                var_cut = f"{cons[0]} > {cons[2]} && {cons[0]} < {cons[3]}"
                #xlo, xhi = ranges[var]

                sig_node = sim_rdf.Define(var, cons[0]).Filter(var_cut)
                bkg_node = proxy["rdf"].Define(var, cons[0]).Filter(var_cut)

                booked[bkg_name][var] = {
                    "h_sig":
                    sig_node.Histo1D(
                        ROOT.RDF.TH1DModel(f"h_{var}_sig_{bkg_name}", "",
                                           cons[1], cons[2], cons[3]), var),
                    "h_bkg":
                    bkg_node.Histo1D(
                        ROOT.RDF.TH1DModel(f"h_{var}_bkg_{bkg_name}", "",
                                           cons[1], cons[2], cons[3]), var),
                }

        # ==============================================================
        # PHASE 2: trigger ALL of it in one go
        # ==============================================================
        all_results = [
            h for d in booked.values() for v in d.values() for h in v.values()
        ]
        print(
            f"\nRunning {len(all_results)} histograms in a single pass per sample..."
        )
        ROOT.RDF.RunGraphs(all_results)

        # ==============================================================
        # PHASE 3: style, normalize, separation, plot (no data reads)
        # ==============================================================
        separations = {bkg_name: dict() for bkg_name in bkg_proxies}

        for bkg_name, proxy in bkg_proxies.items():

            plotfile_name = f"./plots/{sample}_bdtvars_sim_vs_{bkg_name}.pdf"
            c1.SaveAs(f"{plotfile_name}[")

            for var, cons in variables.items():

                h_sig = booked[bkg_name][var]["h_sig"].GetValue()
                h_bkg = booked[bkg_name][var]["h_bkg"].GetValue()

                h_sig.SetLineColor(sig_line_color)
                h_sig.SetLineStyle(sig_line_style)
                h_sig.SetLineWidth(2)
                h_sig.SetMarkerSize(0)
                h_sig.SetTitle("Run 3 sim (sig)")

                h_bkg.SetLineColor(bkg_line_color)
                h_bkg.SetLineStyle(bkg_line_style)
                h_bkg.SetLineWidth(2)
                h_bkg.SetMarkerSize(0)
                h_bkg.SetTitle(proxy["legend"])

                if h_sig.GetEntries() > 0:
                    h_sig.Scale(1.0 / h_sig.GetEntries())
                if h_bkg.GetEntries() > 0:
                    h_bkg.Scale(1.0 / h_bkg.GetEntries())

                separations[bkg_name][var] = separation(h_sig, h_bkg)
                print(
                    f"[{bkg_name}] {var} - separation <S^2> = {separations[bkg_name][var]:.4f}"
                )

                ymax = 1.05 * max(h_sig.GetMaximum(), h_bkg.GetMaximum())
                h_sig.GetYaxis().SetRangeUser(0, ymax)
                h_bkg.GetYaxis().SetRangeUser(0, ymax)

                c1.cd()
                plotInCanvas([h_sig, h_bkg], plotfile_name, cons[5], True,
                             cons[4])

            c1.SaveAs(f"{plotfile_name}]")

        # ==============================================================
        # PHASE 4: two-column table, sorted by the SS separation
        # ==============================================================
        def in_table(var):
            expr = variables[var][0]
            if expr.endswith("_M"):
                return False
            if "_PID_" in expr:
                return False
            return True

        sorted_vars = [(var, sep_ss) for var, sep_ss in sorted(
            separations["ssdata"].items(), key=lambda x: x[1], reverse=True)
                       if in_table(var)]

        top10 = sorted_vars[:10]
        rest = sorted_vars[10:]

        print("\nTop 10 variables")
        print("=" * 45)
        print(f"{'Variable':<20} {'sep vs SS':>10} {'sep vs OS sb':>12}")
        print("-" * 45)
        for var, sep_ss in top10:
            sep_os = separations["osdata_sideband"][var]
            print(f"{var:<20} {sep_ss:10.4f} {sep_os:12.4f}")
        print("=" * 45)

        print("\nRemaining variables")
        print("=" * 45)
        print(f"{'Variable':<20} {'sep vs SS':>10} {'sep vs OS sb':>12}")
        print("-" * 45)
        for var, sep_ss in rest:
            sep_os = separations["osdata_sideband"][var]
            print(f"{var:<20} {sep_ss:10.4f} {sep_os:12.4f}")
        print("=" * 45)

    c1.Close
