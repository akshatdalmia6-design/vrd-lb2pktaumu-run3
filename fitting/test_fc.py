from tabulate import tabulate

import ROOT
import numpy as np

# Import the model
from LFV_model import LFVModel, compileBaseCode
from Plotting import findIntervalFromFC

import json
import sys

from multiprocessing import Pool
import functools

CL_LEVELS = {
    "95CL": {
        "alpha": 1.0 - 0.95,
        "label": "95% CL",
        "color": ROOT.kRed + 1,
        "style": ROOT.kDashed
    },
    "90CL": {
        "alpha": 1.0 - 0.90,
        "label": "90% CL",
        "color": ROOT.kGreen + 2,
        "style": ROOT.kDotted
    },
}

PRIMARY_CL = "95CL"


def _interpolate_crossing(x1, y1, x2, y2, thr):

    if y2 == y1:
        return None, None
    t = (thr - y1) / (y2 - y1)
    return x1 + t * (x2 - x1), t


def find_fc_limits(hist, alpha):

    n = hist.GetNbinsX()
    x = np.array([hist.GetXaxis().GetBinCenter(i) for i in range(1, n + 1)])
    y = np.array([hist.GetBinContent(i) for i in range(1, n + 1)])
    e = np.array([hist.GetBinError(i) for i in range(1, n + 1)])

    ipeak = int(np.argmax(y))
    out = {
        "best": float(x[ipeak]),
        "lower": None,
        "upper": None,
        "lower_err": None,
        "upper_err": None,
        "n_bins_above": int((y > alpha).sum())
    }

    # ---- upper edge: walk right from the peak ----
    for i in range(ipeak, n - 1):
        if y[i] >= alpha and y[i + 1] < alpha:
            xc, t = _interpolate_crossing(x[i], y[i], x[i + 1], y[i + 1],
                                          alpha)
            if xc is None:
                break
            slope = (y[i + 1] - y[i]) / (x[i + 1] - x[i])
            # error on the interpolated p-value at the crossing
            dp = np.hypot((1.0 - t) * e[i], t * e[i + 1])
            out["upper"] = float(xc)
            out["upper_err"] = float(abs(dp / slope)) if slope != 0 else None
            break

    # ---- lower edge: walk left from the peak ----
    for i in range(ipeak, 0, -1):
        if y[i] >= alpha and y[i - 1] < alpha:
            xc, t = _interpolate_crossing(x[i], y[i], x[i - 1], y[i - 1],
                                          alpha)
            if xc is None:
                break
            slope = (y[i - 1] - y[i]) / (x[i - 1] - x[i])
            dp = np.hypot((1.0 - t) * e[i], t * e[i - 1])
            out["lower"] = float(xc)
            out["lower_err"] = float(abs(dp / slope)) if slope != 0 else None
            break

    # If the curve never falls below alpha on the left, the interval runs to
    # the physical boundary -> a pure upper limit.
    if out["lower"] is None and y[0] >= alpha:
        out["lower"] = float(hist.GetXaxis().GetXmin())
        out["lower_err"] = 0.0

    return out


# Global wrapper for 1D FH scan point
def scan_point_fh_worker(args_tuple):
    # Unpack including the new global_min_nll
    data, pdf, fh_var, val_fh, n_toys, output_name, min_nll_fh = args_tuple

    # Create local instance
    # Before f"{output_name}_{val_fh}_temp.root"
    fc = ROOT.FeldmanCousins(data, pdf, fh_var, "")
    ROOT.RooRandom.randomGenerator().SetSeed(42)

    # INJECT THE MINIMUM HERE
    fc.setMinNLLFH(min_nll_fh)

    # Alternatively, use the new C++ method
    fc.setRandomSeed(42)

    # Now pointInFH will use the correct m_min_nll_FH for delta_nll calculation
    pval, perr = fc.pointInFH(val_fh, n_toys)

    return val_fh, pval, perr


# Global wrapper for 1D scan in AFB for 2D model
def scan_point_2d_afb_worker(args_tuple):
    # Unpack including the new global_min_nll
    data, pdf, fh_var, afb_var, val_afb, n_toys, output_name, min_nll_afb, start_params = args_tuple

    # Set values on the shared PDF object before constructor
    set_variables(pdf, start_params)
    ROOT.RooRandom.randomGenerator().SetSeed(42)

    # Create local instance
    # Before f"{output_name}_{val_afb}_temp.root"
    fc = ROOT.FeldmanCousins(data, pdf, fh_var, afb_var, "")

    # INJECT THE MINIMUM HERE
    fc.setMinNLLAFB(min_nll_afb)

    # Alternatively, use the new C++ method
    fc.setRandomSeed(42)

    # Now pointInAFB will use the correct m_min_nll_AFB for delta_nll calculation
    pval, perr = fc.pointInAFB(val_afb, n_toys)

    return val_afb, pval, perr


# Global wrapper for 1D scan in FH for 2D model
def scan_point_2d_fh_worker(args_tuple):
    # Unpack including the new global_min_nll
    data, pdf, fh_var, afb_var, val_fh, n_toys, output_name, min_nll_fh, start_params = args_tuple

    # Set values on the shared PDF object before constructor
    set_variables(pdf, start_params)
    ROOT.RooRandom.randomGenerator().SetSeed(42)

    # Create local instance
    # Before f"{output_name}_{val_fh}_temp.root"
    fc = ROOT.FeldmanCousins(data, pdf, fh_var, afb_var, "")

    # INJECT THE MINIMUM HERE
    fc.setMinNLLFH(min_nll_fh)

    # Alternatively, use the new C++ method
    fc.setRandomSeed(42)

    # Now pointInFH will use the correct m_min_nll_FH for delta_nll calculation
    pval, perr = fc.pointInFH(val_fh, n_toys)

    return val_fh, pval, perr


# Global wrapper for 2D scan point
def scan_point_2d_worker(args_tuple):
    # Unpack values
    data, pdf, fh_var, afb_var, val_fh, val_afb, n_toys, output_name, min_nll, start_params = args_tuple

    # Set values on the shared PDF object before constructor
    set_variables(pdf, start_params)
    ROOT.RooRandom.randomGenerator().SetSeed(42)

    # Create local instance
    fc = ROOT.FeldmanCousins(data, pdf, fh_var, afb_var, "")

    # INJECT THE MINIMUM HERE
    fc.setMinNLL(min_nll)

    # Alternatively, use the new C++ method
    fc.setRandomSeed(42)

    # Perform scans
    pval, perr = fc.pointInPlane(val_fh, val_afb, n_toys)

    return val_fh, val_afb, pval, perr


def get_cache_variables(pdf):
    # 1. Obtain current parameters from the model/pdf
    all_vars = pdf.getVariables()
    starting_params = {}

    # RooArgSet is now iterable directly in Python
    for var in all_vars:
        if not var.isConstant():
            starting_params[var.GetName()] = var.getVal()

    print("Floating parameters at start of FC scan:")
    for k, v in starting_params.items():
        print(f"   {k:<20} = {v:.6g}")

    return starting_params


def set_variables(pdf, start_params):
    # (In some ROOT versions, this helps the constructor's fit)
    vars = pdf.getVariables()
    for name, value in start_params.items():
        v = vars.find(name)
        if v: v.setVal(value)


def plot_cl_scan(histNsig, bins, limits, plot_name):

    canvas = ROOT.TCanvas('canvas_nsig', 'canvas_nsig', 600, 400)

    histNsig.SetMaximum(1.05)
    histNsig.GetYaxis().SetTitle("CL")
    histNsig.Draw("E")

    keep = []

    # threshold lines, correctly labelled
    for key, settings in CL_LEVELS.items():
        ln = ROOT.TLine(bins[1], settings["alpha"], bins[2], settings["alpha"])
        ln.SetLineColor(settings["color"])
        ln.SetLineStyle(settings["style"])
        ln.SetLineWidth(2)
        ln.Draw("SAME")
        keep.append(ln)

    histNsig.Draw("SAME HISTC")

    # vertical marker at the extracted upper limit
    lim = limits[PRIMARY_CL]
    if lim["upper"] is not None:
        vl = ROOT.TLine(lim["upper"], 0.0, lim["upper"], 1.05)
        vl.SetLineColor(ROOT.kRed + 1)
        vl.SetLineStyle(ROOT.kSolid)
        vl.SetLineWidth(2)
        vl.Draw("SAME")
        keep.append(vl)

        mk = ROOT.TMarker(lim["upper"], CL_LEVELS[PRIMARY_CL]["alpha"], 29)
        mk.SetMarkerColor(ROOT.kRed + 1)
        mk.SetMarkerSize(2.0)
        mk.Draw("SAME")
        keep.append(mk)

    leg = ROOT.TLegend(0.52, 0.62, 0.92, 0.90)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.038)
    leg.AddEntry(histNsig, "FC scan", "lp")
    for key, settings in CL_LEVELS.items():
        leg.AddEntry(keep[list(CL_LEVELS).index(key)], settings["label"], "l")
    if lim["upper"] is not None:
        err = lim["upper_err"]
        txt = f"#it{{N}}_{{sig}} < {lim['upper']:.2f}"
        leg.AddEntry(0, txt, "")
    leg.Draw()
    keep.append(leg)

    canvas.Update()
    canvas.SaveAs(f"{plot_name}_nsig_1DCL.pdf")
    print(f"Created: {plot_name}_nsig_1DCL.pdf")


def plot_confidence_belt(fc_main, bins, cl, obs_nsig, limits, plot_name):

    belt = fc_main.getBeltGraph()
    if not belt:
        print("plot_confidence_belt: no belt graph available, skipping")
        return

    canvas = ROOT.TCanvas("canvas_band", "", 600, 500)
    canvas.SetTicks(1, 1)
    canvas.SetLeftMargin(0.14)
    canvas.SetBottomMargin(0.14)

    frame = canvas.DrawFrame(bins[1], bins[1], bins[2], bins[2])
    frame.GetXaxis().SetTitle("#hat{#it{N}}_{sig}  (measured)")
    frame.GetYaxis().SetTitle("#it{N}_{sig}  (true)")
    frame.GetXaxis().SetTitleSize(0.050)
    frame.GetYaxis().SetTitleSize(0.050)
    frame.GetYaxis().SetTitleOffset(1.25)

    belt.SetFillColorAlpha(ROOT.kAzure - 9, 0.60)
    belt.SetLineColor(ROOT.kAzure + 2)
    belt.SetLineWidth(2)
    belt.Draw("F SAME")
    belt.Draw("L SAME")

    keep = []

    # diagonal: where measured == true
    diag = ROOT.TLine(bins[1], bins[1], bins[2], bins[2])
    diag.SetLineStyle(ROOT.kDotted)
    diag.SetLineColor(ROOT.kGray + 2)
    diag.Draw("SAME")
    keep.append(diag)

    # observed measured value -> read the interval off vertically
    vl = None
    do_vl = False
    if do_vl:
        vl = ROOT.TLine(obs_nsig, bins[1], obs_nsig, bins[2])
        vl.SetLineColor(ROOT.kRed + 1)
        vl.SetLineWidth(2)
        vl.Draw("SAME")
        keep.append(vl)

    lim = limits[PRIMARY_CL]
    hl = None
    if lim["upper"] is not None:
        x_stop = obs_nsig if do_vl else bins[2]
        hl = ROOT.TLine(bins[1], lim["upper"], x_stop, lim["upper"])
        hl.SetLineColor(ROOT.kRed + 1)
        hl.SetLineStyle(ROOT.kDashed)
        hl.SetLineWidth(2)
        hl.Draw("SAME")
        keep.append(hl)

    leg = ROOT.TLegend(0.50, 0.16, 0.92, 0.38)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.036)
    leg.AddEntry(belt, f"{int(cl * 100)}% acceptance belt", "f")
    if vl is not None:
        leg.AddEntry(vl, f"observed #hat{{#it{{N}}}}_{{sig}} = {obs_nsig:.2f}",
                     "l")
    if hl is not None:
        leg.AddEntry(hl, f"upper limit = {lim['upper']:.2f}", "l")
    leg.Draw()
    keep.append(leg)

    canvas.Update()
    canvas.SaveAs(f"{plot_name}_nsig_band_{int(cl * 100)}CL.pdf")
    print(f"Created: {plot_name}_nsig_band_{int(cl * 100)}CL.pdf")


#
# Feldman-Cousins approach
#
def doFeldmanCousinsScan(
        data,
        model,
        bins=None,
        nFits=1000,
        plot_name="test",
        output_name=None,
        set_param=None,
        processes=10,
):
    '''
    Performs Feldman-Cousins scan for a given model and data.
    '''

    # Use jitting to load extra functionality
    compileBaseCode()

    pdf = model.getCombinedPdf()

    Nsig = model.getNsig()

    print("Performing Feldman Cousins scan")

    if set_param is not None:
        # Set again original parameters
        model.setParam(set_param)

    # ROOT.FeldmanCousins.debug = True
    ROOT.gErrorIgnoreLevel = ROOT.kSysError
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)

    # 1. Capture the "Best Fit" parameters from the model/pdf
    start_params = get_cache_variables(pdf)

    if len(start_params) == 0:
        raise RuntimeError(
            "No floating parameters: every parameter is fixed by --param. "
            "At least the combinatorial yield must float.")

    ROOT.RooRandom.randomGenerator().SetSeed(42)

    # 2. Pre-scan for global minimum
    fc_main = ROOT.FeldmanCousins(data, pdf, Nsig, f"{output_name}_main.root")

    min_nll = fc_main.getMinNLL()
    print(f"Unconditional NLL minimum = {min_nll:.4f}")

    # Sets the scan binning
    if bins is None:
        bins = [20, 0, 20]

    fc_main.setBins(*bins)
    #

    print("Performing 1D scan in Nsig.")

    histNsig = fc_main.makeNsigHistogram(nFits)

    limits = {
        key: find_fc_limits(histNsig, settings["alpha"])
        for key, settings in CL_LEVELS.items()
    }

    val_nsig = fc_main.getBestNsig()

    plot_cl_scan(histNsig, bins, limits, plot_name)

    print(
        "\nPerforming 1D scan in Nsig with confidence band of Nsig true vs. Nobs."
    )
    cl_band = 1.0 - CL_LEVELS[PRIMARY_CL]["alpha"]
    h_band = fc_main.makeConfidenceBandPlot(nFits, cl_band)
    plot_confidence_belt(fc_main, bins, cl_band, val_nsig, limits, plot_name)

    # Print the intervals and the values with uncertaities
    low_nsig, upp_nsig = findIntervalFromFC(histNsig)

    print("\n" + "=" * 66)
    print("FELDMAN-COUSINS INTERVAL ON N_sig")
    print("=" * 66)
    print(f"  constrained best fit   N_sig_hat = {val_nsig:.3f}")
    print(f"  scan range             [{bins[1]}, {bins[2]}] in {bins[0]} bins")
    print(f"  toys per scan point    {nFits}")
    for key in ("95CL", "90CL"):
        lim = limits[key]
        lbl = CL_LEVELS[key]["label"]
        if lim["upper"] is None:
            print(f"  {lbl:<8} NO CROSSING -- the CL curve stays above "
                  f"{CL_LEVELS[key]['alpha']:.2f} out to N_sig = {bins[2]}. "
                  f"Widen the scan range.")
            continue
        err = lim["upper_err"]
        errs = f" +/- {err:.3f} (toy stat.)" if err is not None else ""
        lo = lim["lower"]
        if lo is not None and lo > bins[1] + 1e-9:
            print(f"  {lbl:<8} interval  [{lo:.3f}, {lim['upper']:.3f}]"
                  f"   (two-sided){errs}")
        else:
            print(f"  {lbl:<8} N_sig < {lim['upper']:.3f}{errs}")
    print(f"  (bin-edge cross-check: [{low_nsig:.3f}, {upp_nsig:.3f}])")
    print("=" * 66 + "\n")

    _lo = limits[PRIMARY_CL][
        "lower"] if limits[PRIMARY_CL]["lower"] is not None else 0.0
    _up = limits[PRIMARY_CL][
        "upper"] if limits[PRIMARY_CL]["upper"] is not None else float('nan')

    print('\n')
    print(
        tabulate([['Nsig', '%.4f -- %.4f' % (_lo, _up)]],
                 headers=['Observable', 'Interval'],
                 tablefmt='presto'))

    print('\n')
    print(
        tabulate([[
            'Nsig',
            '{%.4f}^{+%.4f}_{-%.4f}' %
            (val_nsig, _up - val_nsig, val_nsig - _lo)
        ]],
                 headers=['Observable', 'Value'],
                 tablefmt='presto'))
    print('\n')

    if output_name:

        lim = limits[PRIMARY_CL]
        output = {
            "cl": CL_LEVELS[PRIMARY_CL]["label"],
            "alpha": CL_LEVELS[PRIMARY_CL]["alpha"],
            "low_nsig": lim["lower"],
            "upp_nsig": lim["upper"],
            "upp_nsig_err_toystat": lim["upper_err"],
            "val_nsig": val_nsig,
            "n_toys_per_point": nFits,
            "scan_bins": bins,
            # both levels, for convenience
            "limits": {k: limits[k]
                       for k in limits},
        }

        output_file = f"{output_name}.json"

        with open(output_file, 'w') as json_file:
            json.dump(output, json_file, indent=4)
        print(f"Created: {output_file}")

        output_root_file = ROOT.TFile(f"{output_name}_hists.root", "recreate")
        output_root_file.cd()
        histNsig.Write()
        h_band.Write()
        belt = fc_main.getBeltGraph()
        if belt:
            belt.Write()
        output_root_file.Close()
        print(f"Created: {output_name}_hists.root")

    return limits
