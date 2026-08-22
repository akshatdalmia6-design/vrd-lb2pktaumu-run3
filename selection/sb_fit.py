from pathlib import Path
import numpy as np
import ROOT
import sys

sys.path.append("../utils")
from lhcbPlotStyle import setLHCbPlotStyle

from bdt_run3 import (
    INPUT_FEATURES,
    AUX_FEATURES,
    rdf_to_frame,
    rms90,
    skim_dir,
    NSIGMA,
    BKG_MODE,
    N_SPLITS,
    MODEL_STEM,
    OUT_DIR,
    x_min,
    x_max,
)
import bdt_util as bu
from mass_vetoes import mass_vetoes, run_veto_pipeline

ROOT.gROOT.SetBatch(True)

SAMPLE = "taumu"
FIT_NBINS = 80
x_min2 = 4800
x_max2 = 7000
BDT_CUT = 0.87

PLOT_DIR = OUT_DIR / "plots" / BKG_MODE
CACHE_DIR = OUT_DIR / "cache"
CACHE_FILE = CACHE_DIR / f"sideband_fit_{BKG_MODE}.npz"
GEN_COUNTS = skim_dir / f"{SAMPLE}_gen_counts.npz"

NEW_CUTS = "(tau_END_VZ - Lb_END_VZ > -1) && (Lb_Hlt1TrackMVADecision_TOS || Lb_Hlt1TwoTrackMVADecision_TOS) && (proton_PROBNN_P > 0.4) && (proton_PROBNN_K < 0.6) && (proton_PROBNN_PI < 0.4) && (kaon_PROBNN_PI < 0.5) && (kaon_PROBNN_K > 0.4) && (mu_PROBNN_MU > 0.5)"


def store_count(n_vtx, cut_name, gen_counts):
    store = {}
    if gen_counts.exists():
        d = np.load(gen_counts)
        store = {k: d[k] for k in d.files}
    store[cut_name] = n_vtx
    np.savez(gen_counts, **store)
    print(f"Calculated and stored signal count {cut_name}: {n_vtx}")


def get_windows():
    sim = ROOT.RDataFrame("DecayTree", str(skim_dir / f"{SAMPLE}_sim_OS.root"))
    h = sim.Filter(f"Lb_DTF_M > {x_min} && Lb_DTF_M < {x_max}").Histo1D(
        ROOT.RDF.TH1DModel("h_fine", "", 3500, x_min, x_max), "Lb_DTF_M")
    m_peak = h.GetBinCenter(h.GetMaximumBin())
    sigma = rms90(h.GetValue())
    lo, hi = m_peak - NSIGMA * sigma, m_peak + NSIGMA * sigma
    print(f"DTF peak {m_peak:.1f} MeV, RMS90 {sigma:.1f} MeV")
    print(f"signal window [{lo:.1f}, {hi:.1f}] MeV | "
          f"upper sideband [{hi:.1f}, {x_max2}] MeV | "
          f"lower sideband [{x_min2}, {lo:.1f}] MeV")
    sim = sim.Filter(f"Lb_DTF_M > {lo} && Lb_DTF_M < {hi}")
    n_vtx = int(sim.Count().GetValue())
    store_count(n_vtx, "n_vtx", GEN_COUNTS)
    return lo, hi


def _predict_one(model, X):
    bi = getattr(model, "best_iteration", None)
    if bi is not None:
        return model.predict_proba(X, iteration_range=(0, bi + 1))[:, 1]
    return model.predict_proba(X)[:, 1]


def ensemble_score(models, X):
    X = np.asarray(X, dtype=np.float32)
    return np.mean([_predict_one(m, X) for m in models], axis=0)


def score_full_sideband(models, window_hi, window_lo):
    cut = f"(Lb_DTF_M > {window_hi} && Lb_DTF_M < {x_max2}) || (Lb_DTF_M > {x_min2} && Lb_DTF_M < {window_lo})"
    rdf = ROOT.RDataFrame(
        "DecayTree",
        str(skim_dir / f"{SAMPLE}_dat_{BKG_MODE.upper()}.root")).Filter(cut)
    rdf = rdf.Filter(NEW_CUTS)
    rdf = run_veto_pipeline(rdf)
    rdf = rdf.Filter("passVetoes")
    df = rdf_to_frame(rdf, INPUT_FEATURES, AUX_FEATURES)
    feat = [n for n, _ in INPUT_FEATURES]
    df = df.dropna(subset=feat)
    mass = df["Lb_DTF_M_aux"].to_numpy(dtype=float)
    score = ensemble_score(models, df[feat].to_numpy())
    return mass, score


def score_middle(models, window_hi, window_lo):
    cut = f"(Lb_DTF_M < {window_hi} && Lb_DTF_M > {window_lo})"
    rdf = ROOT.RDataFrame(
        "DecayTree",
        str(skim_dir / f"{SAMPLE}_dat_{BKG_MODE.upper()}.root")).Filter(cut)
    rdf = rdf.Filter(NEW_CUTS)
    rdf = run_veto_pipeline(rdf)
    rdf = rdf.Filter("passVetoes")
    df = rdf_to_frame(rdf, INPUT_FEATURES, AUX_FEATURES)
    feat = [n for n, _ in INPUT_FEATURES]
    df = df.dropna(subset=feat)
    mass = df["Lb_DTF_M_aux"].to_numpy(dtype=float)
    score = ensemble_score(models, df[feat].to_numpy())
    return mass, score


def fit_sideband(mass, window_lo, window_hi):
    h = ROOT.TH1D("h_sb", "", FIT_NBINS, x_min2, x_max2)
    h.Sumw2()
    for m_ in mass:
        h.Fill(float(m_))
    for b in range(1, h.GetNbinsX() + 1):
        lo, hi = h.GetBinLowEdge(b), h.GetBinLowEdge(b) + h.GetBinWidth(b)
        fully_lower = hi <= window_lo
        fully_upper = lo >= window_hi
        if not (fully_lower or fully_upper):
            h.SetBinContent(b, 0.0)
            h.SetBinError(b, 0.0)
    f = ROOT.TF1("f_expo", "expo", x_min2, x_max2)  # exp(p0 + p1*x)
    res = h.Fit(f, "QNRS")
    chi2 = res.Chi2()
    ndf = res.Ndf()
    chi2_ndf = chi2 / ndf if ndf > 0 else float("nan")
    lam = -f.GetParameter(1)
    R = f.Integral(window_lo, window_hi) / (
        f.Integral(window_hi, x_max2) + f.Integral(x_min2, window_lo))
    print(f"slope lambda = {lam:.3e} MeV   extrapolation factor R = {R:.4f}")
    print(f"chi2/ndf = {chi2:.1f} / {ndf} = {chi2_ndf:.2f}")
    return h, f, R, lam, chi2_ndf


def plot_sideband(h_fit, f, window_lo, window_hi, out_path):
    c = ROOT.TCanvas("c1", "c1", 600, 500)
    c.SetTopMargin(0.06)
    c.SetRightMargin(0.05)
    c.SetBottomMargin(0.14)
    c.SetLeftMargin(0.16)
    c.SetTicks(1, 1)
    c.cd()

    def _slice(name, m_lo, m_hi, lower):
        hn_bins = []
        for b in range(1, h_fit.GetNbinsX() + 1):
            lo = h_fit.GetBinLowEdge(b)
            hi = lo + h_fit.GetBinWidth(b)
            keep = (hi <= window_lo) if lower else (lo >= window_hi)
            if keep:
                hn_bins.append(b)
        if not hn_bins:
            return ROOT.TH1D(name, "", 1, m_lo, m_hi)
        b1, b2 = hn_bins[0], hn_bins[-1]
        hn = ROOT.TH1D(name, "", b2 - b1 + 1, h_fit.GetBinLowEdge(b1),
                       h_fit.GetBinLowEdge(b2) + h_fit.GetBinWidth(b2))
        hn.SetDirectory(0)
        for i, b in enumerate(range(b1, b2 + 1), 1):
            hn.SetBinContent(i, h_fit.GetBinContent(b))
            hn.SetBinError(i, h_fit.GetBinError(b))
        return hn

    h_low = _slice("h_low", x_min2, window_lo, True)
    h_up = _slice("h_up", window_hi, x_max2, False)

    ymax = 1.45 * max(h_low.GetMaximum(), h_up.GetMaximum())
    frame = c.DrawFrame(x_min2, 0.0, x_max2, ymax)
    frame.GetXaxis().SetTitle("m_{DTF} [MeV]")
    frame.GetYaxis().SetTitle("Candidates")
    for ax in (frame.GetXaxis(), frame.GetYaxis()):
        ax.SetTitleSize(0.055)
        ax.SetLabelSize(0.045)
    frame.GetXaxis().SetTitleOffset(1.10)
    frame.GetXaxis().SetNdivisions(505)
    frame.GetYaxis().SetTitleOffset(1.30)
    frame.GetYaxis().SetMaxDigits(3)

    for hh in (h_low, h_up):
        hh.SetMarkerStyle(20)
        hh.SetMarkerColor(ROOT.kBlack)
        hh.SetLineColor(ROOT.kBlack)
    h_low.Draw("E SAME")
    h_up.Draw("E SAME")

    f.SetRange(x_min2, x_max2)
    f.SetLineColor(ROOT.kBlue + 1)
    f.SetLineWidth(3)
    f.SetLineStyle(1)
    f.Draw("SAME")

    line_lo = ROOT.TLine(window_lo, 0, window_lo, ymax)
    line_hi = ROOT.TLine(window_hi, 0, window_hi, ymax)
    for ln in (line_lo, line_hi):
        ln.SetLineStyle(3)
        ln.SetLineColor(ROOT.kGray + 2)
        ln.Draw("SAME")

    leg = ROOT.TLegend(0.56, 0.74, 0.93, 0.90)
    leg.AddEntry(h_low, f"{BKG_MODE.upper()} sidebands", "lep")
    leg.AddEntry(f, "exponential fit", "l")
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.045)
    leg.Draw()

    c.RedrawAxis()
    c.SaveAs(str(out_path))
    print(f"saved {out_path}")


def main():
    setLHCbPlotStyle()
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    CACHE_DIR.mkdir(parents=True, exist_ok=True)

    print(f"=== Sideband fit: {BKG_MODE.upper()} ===")

    window_lo, window_hi = get_windows()
    models = bu.load_models(str(MODEL_STEM), N_SPLITS)
    print(f"loaded {len(models)} fold models from {MODEL_STEM}_fold*.json")

    mass, score = score_full_sideband(models, window_hi, window_lo)
    mass = mass[score > BDT_CUT]

    h, f, R, lam, chi2_ndf = fit_sideband(mass, window_lo, window_hi)
    plot_sideband(h, f, window_lo, window_hi, PLOT_DIR / "sideband_fit.pdf")

    if BKG_MODE == "ss":
        n_sb = len(mass)
        print(f"SS sideband extrapolated bkg in window = {n_sb * R:.1f}")
        mass2, score2 = score_middle(models, window_hi, window_lo)
        mass2 = mass2[score2 > 0.85]
        print(f"SS direct bkg in window = {len(mass2):.1f}")
        score = score2

    np.savez(
        CACHE_FILE,
        R=R,
        lam=lam,
        chi2_ndf=chi2_ndf,
        score_sb=score,
        bdt_cut=BDT_CUT,
        window_lo=window_lo,
        window_hi=window_hi)
    print(f"saved cache -> {CACHE_FILE}")


if __name__ == "__main__":
    main()
