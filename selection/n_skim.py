import ROOT
import numpy as np
import math
import sys
from pathlib import Path

from simplePlots import get_all_input_files

sys.path.append("../utils")
from lhcbPlotStyle import setLHCbPlotStyle
setLHCbPlotStyle()

# Multi-core: the skim is a pure read+filter+write, parallelises perfectly
ROOT.EnableImplicitMT()
ROOT.gROOT.SetBatch(True)

# Local output dir for the skims - put this somewhere on your work area,
# NOT on EOS, so reruns read locally.
skim_dir = Path("./skims")
skim_dir.mkdir(exist_ok=True)

plot_dir = Path("./plots")
plot_dir.mkdir(exist_ok=True)

gen_tree_paths = ["MCDT/MCDecayTree"]

HLT1_cut = "(Lb_Hlt1TrackMVADecision_TOS || Lb_Hlt1TwoTrackMVADecision_TOS)"

PID_cut = "(H1_PROBNN_P > 0.4) && (H1_PROBNN_K < 0.6) && (H1_PROBNN_PI < 0.4) && (H2_PROBNN_PI < 0.5) && (H2_PROBNN_K > 0.4) && (L1_PROBNN_MU > 0.5) && (L2_PROBNN_MU > 0.5)"

tm_cut = "Lb_BKGCAT < 51"

iso_cut = "Lb_HEAD_CC_Lb_ChargedIso_DR2_0p5_PTASY > -0.5 && Lb_HEAD_NC_Lb_NeutralIso_DR2_0p5_PTASY > -0.5"

x_min = 5400
x_max = 5800

sample_configs = {
    "mumu": {
        "sim_sample_name": "Lb_JpsipK_mm",
        "dat_sample_name": "turbo_rdlow",
        "chain_name": "Hlt2RD_LambdabToPKMuMu",
    }
}

# Columns to keep in the skim: everything used for cuts, the mass window,
# and every candidate BDT variable. Keep this in sync with the variables
# dict in the comparison script.
keep_columns = [
    "Lb_Hlt1TrackMVADecision_TOS",
    "Lb_Hlt1TwoTrackMVADecision_TOS",
    "Lb_Hlt1TrackMuonMVADecision_TOS",
    "H1_PROBNN_P",
    "H1_PROBNN_K",
    "H1_PROBNN_PI",
    "H2_PROBNN_K",
    "H2_PROBNN_PI",
    "L1_PROBNN_MU",
    "L2_PROBNN_MU",
    "Lb_END_VZ",
    #candidates
    "H1_P",
    "H1_PT",
    "H1_PID_P",
    "H1_BPVIPCHI2",
    "H1_PID_K",
    "H2_PT",
    "H2_P",
    "H2_PID_K",
    "H2_BPVIPCHI2",
    "L1_PT",
    "L1_P",
    "L1_PID_MU",
    "L1_BPVIPCHI2",
    "L1_CHI2DOF",
    "L2_PT",
    "L2_P",
    "L2_PID_MU",
    "L2_BPVIPCHI2",
    "L2_CHI2DOF",
    "Lb_M",
    "Lb_CHI2DOF",
    "Lb_P",
    "Lb_PT",
    "Lb_BPVIPCHI2",
    "Lb_BPVFDCHI2",
    "Lb_BPVDIRA",
    "Lb_MAXSDOCACHI2",
    "Lb_HEAD_CC_Lb_ChargedIso_DR2_0p5_PTASY",
    "Lb_HEAD_CC_Lb_ChargedIso_DR2_0p5_PASY",
    "Lb_HEAD_CC_Lb_ChargedIso_DR2_0p5_Max_P",
    "Lb_HEAD_CC_Lb_ChargedIso_DR2_0p5_CMULT",
    "Lb_HEAD_CC_Lb_ChargedIso_DR2_0p5_Max_PT",
    "Lb_HEAD_NC_Lb_NeutralIso_DR2_0p5_CMULT",
    "Lb_HEAD_NC_Lb_NeutralIso_DR2_0p5_PTASY",
    "Lb_HEAD_NC_Lb_NeutralIso_DR2_0p5_PASY",
    "Lb_HEAD_NC_Lb_NeutralIso_DR2_0p5_Max_P",
    "Lb_HEAD_NC_Lb_NeutralIso_DR2_0p5_Max_PT"
]

# Truth-match columns are only present in simulation
truth_columns = ["Lb_BKGCAT"]


def build_chain(sample_name, chain_name):
    chain = ROOT.TChain(f"{chain_name}/DecayTree")
    for sname, in_file in get_all_input_files(sample_name):
        in_file_stem = Path(in_file).name.split('.')[0]
        chain.Add(in_file)
    return chain


def count_gen_entries(sample_name, tree_paths):
    total = 0
    for tp in tree_paths:
        chain = ROOT.TChain(tp)
        for sname, in_file in get_all_input_files(sample_name):
            chain.AddFile(in_file)
        total += int(chain.GetEntries())
    return total


def skim(sample_name, chain_name, cut, out_path, columns):
    """Apply cut, materialize the kept columns, write a small local file."""
    chain = build_chain(sample_name, chain_name)
    rdf = ROOT.RDataFrame(chain).Filter(cut)

    # Only keep columns that actually exist in this chain
    available = set(str(c) for c in rdf.GetColumnNames())
    cols = [c for c in columns if c in available]
    missing = [c for c in columns if c not in available]
    if missing:
        print(f"  (skipping {len(missing)} absent columns: {missing})")

    print(f"Writing skim -> {out_path}")
    rdf.Snapshot("DecayTree", str(out_path), cols)


TABLE_ROWS = [
    ("Generated", "n_gen"),
    ("HLT2 trigger", "n_reco"),
    ("Signal window", "n_win"),
    ("Fiducial cuts", "n_vtx"),
    ("delZ cuts", "n_end_vz"),
    ("HLT1 cuts", "n_hlt1"),
    ("PID cuts", "n_pid"),
]

GEN_COUNTS = skim_dir / f"mumu_norm_gen_counts.npz"


def efficiency_table():
    """2-column table: selection name, successive efficiency N_i / N_{i-1}.
    Final row is the total efficiency sig_bdt / total."""
    d = np.load(GEN_COUNTS)
    rows = [(lbl, key) for lbl, key in TABLE_ROWS if key in d.files]

    def binom_err(k, n):
        if n <= 0:
            return float("nan")
        e = k / n
        return (e * (1.0 - e) / n)**0.5

    N = [float(d[key]) for _, key in rows]
    N0 = N[0]
    header = f"{'Selection':<22}{'Efficiency':>24}"
    sep = "-" * len(header)
    lines = [header, sep]

    for i in range(1, len(rows)):
        lbl = rows[i][0]
        prev, n = N[i - 1], N[i]
        e = n / prev if prev > 0 else float("nan")
        err = binom_err(n, prev)
        lines.append(f"{lbl:<22}{e:>12.5f} +/- {err:<8.5f}")
    lines.append(sep)
    total = N[-1] / N0 if N0 > 0 else float("nan")
    total_err = binom_err(N[-1], N0)
    lines.append(f"{'TOTAL':<22}{total:>12.3e} +/- {total_err:<8.3e}")

    table = "\n".join(lines)
    print("\n" + table)


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


variables = {
    "Lb_Mass": [
        "Lb_M", 100, 5350, 6050, "#it{m} (#Lambda_{b}) [MeV/#it{c}^{2}]",
        "[MeV/#it{c}^{2}]"
    ],
}

WINDOW_CUT = f"Lb_M > {x_min} && Lb_M < {x_max}"
sig_line_color = ROOT.kRed
sig_line_style = ROOT.kSolid
bkg_line_color = ROOT.kBlue
bkg_line_style = ROOT.kDashed


def plot_mass(sample):

    sim_rdf = ROOT.RDataFrame(
        "DecayTree", str(
            skim_dir / f"{sample}_norm_sim.root")).Filter(WINDOW_CUT)
    dat_rdf = ROOT.RDataFrame(
        "DecayTree", str(
            skim_dir / f"{sample}_norm_dat.root")).Filter(WINDOW_CUT)

    booked = {}
    for var, cons in variables.items():
        expr, nbins, xlo, xhi = cons[0], cons[1], cons[2], cons[3]
        sim_def = sim_rdf.Define(var, expr)
        dat_def = dat_rdf.Define(var, expr)
        booked[var] = {
            "h_sig":
            sim_def.Histo1D(
                ROOT.RDF.TH1DModel(f"h_{var}_sim", "", nbins, xlo, xhi), var),
            "h_bkg":
            dat_def.Histo1D(
                ROOT.RDF.TH1DModel(f"h_{var}_dat", "", nbins, xlo, xhi), var),
        }

    all_results = [h for v in booked.values() for h in v.values()]
    ROOT.RDF.RunGraphs(all_results)  # one event loop per sample

    for var, cons in variables.items():
        h_sig = booked[var]["h_sig"].GetValue()
        h_bkg = booked[var]["h_bkg"].GetValue()

        h_sig.SetLineColor(sig_line_color)
        h_sig.SetLineStyle(sig_line_style)
        h_sig.SetLineWidth(2)
        h_sig.SetMarkerSize(0)
        h_sig.SetTitle("Run 3 sim (norm)")  # legend reads GetTitle()

        h_bkg.SetLineColor(bkg_line_color)
        h_bkg.SetLineStyle(bkg_line_style)
        h_bkg.SetLineWidth(2)
        h_bkg.SetMarkerSize(0)
        h_bkg.SetTitle("Run 3 data (norm)")

        print(f"\n{var}: sim = {h_sig.GetEntries():.0f}, "
              f"data = {h_bkg.GetEntries():.0f}")

        #c unit-area normalisation -- exactly as in compare_fast_skimmed.py
        if h_sig.GetEntries() > 0:
            h_sig.Scale(1.0 / h_sig.GetEntries())
        if h_bkg.GetEntries() > 0:
            h_bkg.Scale(1.0 / h_bkg.GetEntries())

        out = plot_dir / f"{sample}_{var}_sim_vs_dat.pdf"
        c1.cd()
        plotInCanvas([h_sig, h_bkg], str(out), cons[5], True, cons[4])
        print(f"saved {out}")


if __name__ == '__main__':

    for sample in ["mumu"]:

        cfg = sample_configs[sample]

        sim_cut = " && ".join([tm_cut, iso_cut, HLT1_cut, PID_cut])
        dat_cut = " && ".join([iso_cut, HLT1_cut, PID_cut])

        sim_sample_name = cfg["sim_sample_name"]
        dat_sample_name = cfg["dat_sample_name"]
        chain_name = cfg["chain_name"]

        # Make skim files (run first 1 out of 2)

        # 1. Signal: simulation, OS line, truth-matched
        skim(sim_sample_name, chain_name, sim_cut,
             skim_dir / f"{sample}_norm_sim.root",
             keep_columns + truth_columns)

        # 2. Background proxy: data, (full spectrum)
        skim(dat_sample_name, chain_name, dat_cut,
             skim_dir / f"{sample}_norm_dat.root", keep_columns)

        # Make skim files (run second 2 out of 2)

        print(f"\n Counting entries for {sim_sample_name}...")

        sim_chain_c = build_chain(sim_sample_name, chain_name)
        n_reco = sim_chain_c.GetEntries()
        print(f"Reconstructed candidates before sim_cut ({sample}): {n_reco}")

        n_gen = count_gen_entries(sim_sample_name, gen_tree_paths)
        print(f"Total generated signal decays ({sample}): {n_gen}")

        n_win = sim_chain_c.GetEntries(f"Lb_M > {x_min} && Lb_M < {x_max}")
        print(f"Candidates in signal window ({sample}): {n_win}")

        n_vtx = sim_chain_c.GetEntries(
            f"Lb_M > {x_min} && Lb_M < {x_max} && {iso_cut} && {tm_cut}")
        print(f"Candidates after fiducial cuts ({sample}): {n_vtx}")

        n_hlt1 = sim_chain_c.GetEntries(
            f"Lb_M > {x_min} && Lb_M < {x_max} && {iso_cut} && {tm_cut} && {HLT1_cut}"
        )
        print(f"Candidates after HLT1 cuts ({sample}): {n_hlt1}")

        n_pid = sim_chain_c.GetEntries(
            f"Lb_M > {x_min} && Lb_M < {x_max} && {iso_cut} && {tm_cut} && {HLT1_cut} && {PID_cut}"
        )
        print(f"Candidates after PID cuts ({sample}): {n_pid}")

        gen_out = skim_dir / f"{sample}_norm_gen_counts.npz"
        np.savez(
            gen_out,
            n_gen=n_gen,
            n_reco=n_reco,
            n_win=n_win,
            n_vtx=n_vtx,
            n_hlt1=n_hlt1,
            n_pid=n_pid)
        print(f"Saved counts -> {gen_out}")

        plot_mass(sample)

    c1.Close()

    print("\nSkims done.")

    efficiency_table()
