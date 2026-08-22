from pathlib import Path
import numpy as np
import ROOT
import sys

sys.path.append("../utils")
from lhcbPlotStyle import setLHCbPlotStyle

from bdt_run3 import BKG_MODE, OUT_DIR, skim_dir
from sb_fit import store_count
from make_fit_tuples import GEN_3PI, GEN_3PIPI0, BR_3PI, BR_3PIPI0, _load_gen

ROOT.gROOT.SetBatch(True)

A_PUNZI = 5.0
A_HALF = A_PUNZI / 2.0
N_GRID = 101  # threshold grid
T_START = 0.990  # threshold start point
N_MIN = 10  # min background events after cut to trust optimum

BR = {0: BR_3PI, 1: BR_3PIPI0}
BRE = {0: 0.0005, 1: 0.0005}
BR_TOT = BR_3PI + BR_3PIPI0
GEN = {0: GEN_3PI, 1: GEN_3PIPI0}
COV_BR = 0.0

MC_COMBINED = OUT_DIR / "cache" / "scored_mc_combined.root"
PLOT_DIR = OUT_DIR / "plots" / BKG_MODE
CACHE_FILE = OUT_DIR / "cache" / f"sideband_fit_{BKG_MODE}.npz"


def load_cache():
    if not CACHE_FILE.exists():
        raise FileNotFoundError(f"{CACHE_FILE} not found")
    d = np.load(CACHE_FILE)
    R = float(d["R"])
    score_sb = np.asarray(d["score_sb"], dtype=float)
    bdt_cut = float(d["bdt_cut"])
    print(
        f"loaded cache: R={R:.4f}, BDT CUT={bdt_cut:.3f}, {len(score_sb)} background scores,  "
        f"window=[{float(d['window_lo']):.1f}, {float(d['window_hi']):.1f}] MeV"
    )
    return score_sb, R, bdt_cut, float(d["window_lo"]), float(d["window_hi"])


def _eff_one(sorted_scores, n_gen, thresholds):
    n_tot = len(sorted_scores)
    n_pass = n_tot - np.searchsorted(sorted_scores, thresholds, side="right")
    eff = n_pass / n_gen
    eff_bdt = n_pass / n_tot
    err = np.sqrt(np.maximum(eff * (1.0 - eff), 0.0) / n_gen)
    err_bdt = np.sqrt(np.maximum(eff * (1.0 - eff), 0.0) / n_tot)
    return eff, eff_bdt, err, err_bdt, n_pass


def signal_efficiency(thresholds, lo, hi):

    if not MC_COMBINED.exists():
        raise FileNotFoundError(
            f"{MC_COMBINED} not found -- run make_fit_tuples.py")

    rdf = ROOT.RDataFrame(
        "tree", str(MC_COMBINED)).Filter(f"Lb_M > {lo} && Lb_M <{hi}")
    arr = rdf.AsNumpy(["bdt", "mode"])
    bdt = np.asarray(arr["bdt"], dtype=float)
    mode = np.asarray(arr["mode"], dtype=int)

    e, de, e_bdt, de_bdt = {}, {}, {}, {}
    for m in (0, 1):
        sel = (mode == m) & np.isfinite(bdt)
        s = np.sort(bdt[sel])
        n_gen_m = _load_gen(GEN[m])
        e[m], e_bdt[m], de[m], de_bdt[m], n_pass_m = _eff_one(
            s, n_gen_m, thresholds)

        if m == 0:
            n_bdt_0 = n_pass_m
        elif m == 1:
            n_bdt_1 = n_pass_m

        print(f"  mode {m}: {len(s)} in-window / n_gen={n_gen_m:.0f}, "
              f"BR={BR[m]:.4f}+/-{BRE[m]:.4f}")

    eff = (BR[0] * e[0] + BR[1] * e[1]) / BR_TOT
    eff_bdt = (BR[0] * e_bdt[0] + BR[1] * e_bdt[1]) / BR_TOT

    stat = np.sqrt((BR[0] / BR_TOT)**2 * de[0]**2 +
                   (BR[1] / BR_TOT)**2 * de[1]**2)
    stat_bdt = np.sqrt((BR[0] / BR_TOT)**2 * de_bdt[0]**2 +
                       (BR[1] / BR_TOT)**2 * de_bdt[1]**2)

    delta = e[0] - e[1]
    d0 = BR[1] * delta / BR_TOT**2
    d1 = -BR[0] * delta / BR_TOT**2
    syst = np.maximum(
        ((d0 * BRE[0])**2 + (d1 * BRE[1])**2 + 2.0 * d0 * d1 * COV_BR), 0.0)
    syst = np.sqrt(syst)

    return eff, eff_bdt, stat, stat_bdt, syst, n_bdt_0, n_bdt_1


def punzi_scan(eff, eff_err, n_bkg, R_eff, thresholds):
    B = R_eff * n_bkg.astype(float)
    dB = R_eff * np.sqrt(n_bkg.astype(float))
    D = A_HALF + np.sqrt(B)
    fom = eff / D

    with np.errstate(divide="ignore", invalid="ignore"):
        term_eps = eff_err / D
        term_B = np.where(B > 0, eff * dB / (2.0 * D**2 * np.sqrt(B)), 0.0)
    fom_err = np.sqrt(term_eps**2 + term_B**2)

    valid = n_bkg >= N_MIN
    if not valid.any():
        raise RuntimeError(f"no threshold has >= {N_MIN} background events")
        valid = np.ones_like(n_bkg, dtype=bool)
    i = int(np.argmax(np.where(valid, fom, -np.inf)))
    return dict(
        i=i,
        t=thresholds[i],
        eff=eff[i],
        B=B[i],
        dB=dB[i],
        n=int(n_bkg[i]),
        fom=fom[i],
        R_eff=R_eff,
        fom_err=fom_err[i],
        fom_err_arr=fom_err,
        fom_arr=fom,
        B_arr=B,
        valid=valid)


def survivors(scores, thresholds):
    s = np.sort(np.asarray(scores, dtype=float))
    return len(s) - np.searchsorted(s, thresholds, side="right")


def plot_fom(thresholds,
             fom,
             fom_err,
             valid,
             t_opt,
             out_path,
             label,
             opt=False):
    c = ROOT.TCanvas("c1", "c1", 600, 500)
    c.SetTopMargin(0.06)
    c.SetRightMargin(0.05)
    c.SetBottomMargin(0.14)
    c.SetLeftMargin(0.16)
    c.SetTicks(1, 1)
    c.cd()

    y = np.where(np.isfinite(fom), fom, 0.0)
    e = np.where(np.isfinite(fom_err), fom_err, 0.0)
    y_up = y + e
    y_dn = np.maximum(y - e, 0.0)

    g = ROOT.TGraph(
        len(thresholds), np.ascontiguousarray(thresholds, dtype=float),
        np.ascontiguousarray(y, dtype=float))
    g.SetLineWidth(3)
    g.SetLineColor(ROOT.kAzure + 1)
    g.SetMarkerStyle(20)
    g.SetMarkerColor(ROOT.kAzure + 1)
    g.Draw("AL")
    g.GetXaxis().SetTitle("BDT cut t")
    g.GetYaxis().SetTitle(label)
    for ax in (g.GetXaxis(), g.GetYaxis()):
        ax.SetTitleSize(0.055)
        ax.SetLabelSize(0.045)
    g.GetXaxis().SetTitleOffset(1.10)
    g.GetXaxis().SetNdivisions(505)
    g.GetYaxis().SetTitleOffset(1.30)
    g.GetYaxis().SetMaxDigits(3)
    ymax = 1.25 * float(np.max(y_up[valid]))
    g.GetHistogram().SetMinimum(0.0)
    g.GetHistogram().SetMaximum(ymax)

    if opt:
        g_up = ROOT.TGraph(
            len(thresholds), np.ascontiguousarray(thresholds, dtype=float),
            np.ascontiguousarray(y_up, dtype=float))
        g_dn = ROOT.TGraph(
            len(thresholds), np.ascontiguousarray(thresholds, dtype=float),
            np.ascontiguousarray(y_dn, dtype=float))
        for gg in (g_up, g_dn):
            gg.SetLineWidth(2)
            gg.SetLineStyle(2)
            gg.SetLineColor(ROOT.kAzure + 1)
            gg.Draw("L SAME")

    line = ROOT.TLine(t_opt, 0.0, t_opt, ymax)
    line.SetLineStyle(2)
    line.SetLineColor(ROOT.kRed + 1)
    line.SetLineWidth(2)
    line.Draw("SAME")

    if opt:
        leg = ROOT.TLegend(0.20, 0.68, 0.55, 0.83)
        leg.AddEntry(g_up, "#pm 1#sigma", "l")
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextFont(42)
        leg.SetTextSize(0.042)
        leg.Draw()

    lat = ROOT.TLatex()
    lat.SetNDC()
    lat.SetTextFont(42)
    lat.SetTextSize(0.045)
    lat.DrawLatex(0.20, 0.86, f"optimum  BDT > {t_opt:.4f}")

    c.Modified()
    c.Update()
    c.RedrawAxis()
    c.SaveAs(str(out_path))
    print(f"saved {out_path}")


TABLE_ROWS = [
    ("Generated", "n_gen"),
    ("HLT2 trigger", "n_reco"),
    ("Signal window", "n_win"),
    ("Fiducial cuts", "n_vtx"),
    ("delZ cuts", "n_end_vz"),
    ("HLT1 cuts", "n_hlt1"),
    ("PID cuts", "n_pid"),
    ("Mass vetoes", "n_veto"),
    ("BDT selection", "n_bdt"),
]


def efficiency_table():
    """2-column table: selection name, successive efficiency N_i / N_{i-1}.
    Final row is the total efficiency sig_bdt / total."""
    d0 = np.load(GEN[0])
    d1 = np.load(GEN[1])

    rows = [(lbl, key) for lbl, key in TABLE_ROWS
            if key in d0.files and key in d1.files]

    def binom_err(k, n):
        if n <= 0:
            return float("nan")
        e = k / n
        return (e * (1.0 - e) / n)**0.5

    N0 = [float(d0[key]) for _, key in rows]
    N1 = [float(d1[key]) for _, key in rows]
    header = f"{'Selection':<20}{'eps(3pi nu)':>22}{'eps(3pi pi0 nu)':>26}"
    sep = "-" * len(header)
    lines = [header, sep]

    for i in range(1, len(rows)):
        lbl = rows[i][0]
        e0 = N0[i] / N0[i - 1] if N0[i - 1] > 0 else float("nan")
        e1 = N1[i] / N1[i - 1] if N1[i - 1] > 0 else float("nan")
        de0 = binom_err(N0[i], N0[i - 1])
        de1 = binom_err(N1[i], N1[i - 1])
        lines.append(
            f"{lbl:<20}{e0:>11.5f} +/- {de0:<7.5f}{e1:>13.5f} +/- {de1:<7.5f}")
    lines.append(sep)
    tot0 = N0[-1] / N0[0] if N0[0] > 0 else float("nan")
    tot1 = N1[-1] / N1[0] if N1[0] > 0 else float("nan")
    dtot0 = binom_err(N0[-1], N0[0])
    dtot1 = binom_err(N1[-1], N1[0])
    lines.append(
        f"{'TOTAL per-mode':<20}{tot0:>11.3e} +/- {dtot0:<7.3e}{tot1:>13.3e} +/- {dtot1:<7.3e}"
    )

    table = "\n".join(lines)
    print("\n" + table)


def print_gen_counts():  #f
    """Dump everything currently stored in gen_counts.npz (raw), so you can see
    exactly which counts exist before the table is built."""
    if not GEN_COUNTS.exists():
        print(f"{GEN_COUNTS} does not exist yet")
        return
    d = np.load(GEN_COUNTS)
    print(f"\ncontents of {GEN_COUNTS}:")
    for k in d.files:
        v = d[k]
        try:
            print(f"   {k:<12} = {float(v):,.0f}")
        except (TypeError, ValueError):
            print(f"   {k:<12} = {v}")


def main():
    setLHCbPlotStyle()
    PLOT_DIR.mkdir(parents=True, exist_ok=True)

    score_sb, R, bdt_cut, window_lo, window_hi = load_cache()
    thresholds = np.linspace(T_START, 1.0, N_GRID)

    if BKG_MODE == "os":
        R_eff = R
        desc = f"OS sidebands x R (extrapolated), R={R:.4f}"
    elif BKG_MODE == "ss":
        R_eff = 1.0
        desc = "SS in-window direct count (no extrapolation)"

    eff, eff_bdt, stat, stat_bdt, syst, n_bdt_0, n_bdt_1 = signal_efficiency(
        thresholds, window_lo, window_hi)
    n_bkg = survivors(score_sb, thresholds)
    res = punzi_scan(eff, stat, n_bkg, R_eff, thresholds)
    res_bdt = punzi_scan(eff_bdt, stat_bdt, n_bkg, R_eff, thresholds)

    print("\n" + "=" * 89)
    print(f"PUNZI RESULT  (a={A_PUNZI:.0f}, {BKG_MODE} bkg)")
    print("=" * 89)
    print(f"  optimal cut            BDT > {res['t']:.4f}")
    print(
        f"  signal efficiency      eps = {res['eff']:.3e} +/- {stat[res['i']]:.2e}"
    )
    print(
        f"  signal efficiency(bdt) eps = {res_bdt['eff']:.4f} +/- {stat_bdt[res['i']]:.4f}"
    )
    if BKG_MODE == "os":
        print(
            f"  expected bkg           B   = {res['B']:.2f} +/- {res['dB']:.2f} "
            f"({res['n']} sideband events x R={R_eff:.4f})")
    else:
        print(
            f"  expected bkg           B   = {res['B']:.2f} +/- {res['dB']:.2f} "
            f"({res['n']} events in signal window)")
    print(
        f"  FOM                        = {res['fom']:.3e} +/- {res['fom_err']:.2e}"
    )
    print("=" * 89)

    print("Neighbourhood:")
    lo = max(0, res["i"] - 10)
    hi = min(len(thresholds), res["i"] + 11)
    for j in range(lo, hi):
        note = "" if res["valid"][j] else "(N/A)"
        print(
            f"  t={thresholds[j]:.4f}  eps={eff[j]:.3e}  eps(bdt)={eff_bdt[j]:.4f}  B={res['B_arr'][j]:6.2f}  FOM={res['fom_arr'][j]:.3e} +/- {res['fom_err_arr'][j]:.2e} {note}"
        )
    print("=" * 89)

    store_count(n_bdt_0[res["i"]], "n_bdt", GEN[0])
    store_count(n_bdt_1[res["i"]], "n_bdt", GEN[1])

    plot_fom(
        thresholds,
        res["fom_arr"],
        res["fom_err_arr"],
        res["valid"],
        res["t"],
        PLOT_DIR / "punzi_fom.pdf",
        "Punzi FOM #varepsilon / (a/2 + #sqrt{B})",
        opt=True)
    plot_fom(
        thresholds,
        eff,
        stat,
        res["valid"],
        res["t"],
        PLOT_DIR / "punzi_eff.pdf",
        "Signal efficiency #varepsilon",
        opt=False)

    efficiency_table()

    print(
        f"Final combined efficiency = {res['eff']:.3e} +/- {stat[res['i']]:.2e} (stat.) +/- {syst[res['i']]:.2e} (syst.)"
    )


if __name__ == "__main__":
    main()
