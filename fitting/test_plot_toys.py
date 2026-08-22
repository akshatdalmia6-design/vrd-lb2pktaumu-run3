import ROOT
import json
import numpy as np


def getMeanAndAverage(_list):

    np_list = np.array(_list)

    if np_list.size == 0:
        return 0.0, 0.0

    return np.mean(np_list), np.std(np_list)


def _range_around(mean, std, nsigma=5.0, floor=1.0):

    half = max(nsigma * std, floor)
    return mean - half, mean + half


def fitHistogram(histogram, positive_range=False):

    gaus_func = ROOT.TF1(f"gaus_func_{histogram.GetName()}", "gaus",
                         histogram.GetXaxis().GetXmin(),
                         histogram.GetXaxis().GetXmax())
    gaus_func.SetParameter(1, histogram.GetMean())
    gaus_func.SetParameter(2, histogram.GetRMS())
    gaus_func.SetLineColor(ROOT.kBlue)

    if positive_range:
        gaus_func.FixParameter(1, 0.01)
        histogram.Fit(gaus_func, "S", "", 0, histogram.GetXaxis().GetXmax())
    else:
        histogram.Fit(gaus_func, "S")

    ROOT.gPad.Update()
    stats = histogram.FindObject("stats")
    if stats:
        stats.SetX1NDC(0.60)  # Bottom-left X
        stats.SetY1NDC(0.50)  # Bottom-left Y
        stats.SetX2NDC(0.92)  # Top-right X
        stats.SetY2NDC(0.92)  # Top-right Y
        stats.SetFillStyle(0)
    ROOT.gPad.Modified()
    ROOT.gPad.Update()


def plotToyResults(input_file, plot_name, reparam=False):
    '''
    Plots toy results from a given input file.
    '''

    results = dict()

    with open(input_file, 'r') as f:
        # Load json file
        results = json.load(f)

    if len(results['toys']) == 0:
        print("plotToyResults: no toys in the input file")
        return

    nsig_mean, nsig_std = getMeanAndAverage(
        [toy['Nsig']['value'] for toy in results['toys'].values()])
    ncom_mean, ncom_std = getMeanAndAverage(
        [toy['Ncom']['value'] for toy in results['toys'].values()])
    exp_mean, exp_std = getMeanAndAverage(
        [toy['Expcom']['value'] for toy in results['toys'].values()])

    hist_status = ROOT.TH1D('hist_status', 'hist_status', 6, -0.5, 5.5)
    hist_status.SetXTitle('status')

    hist_quality = ROOT.TH1D('hist_quality', 'hist_quality', 6, -0.5, 5.5)
    hist_quality.SetXTitle('cov. quality')

    nsig_lo, nsig_hi = _range_around(nsig_mean, nsig_std, floor=5.0)
    ncom_lo, ncom_hi = _range_around(ncom_mean, ncom_std, floor=5.0)
    exp_lo, exp_hi = _range_around(
        exp_mean, exp_std, floor=abs(exp_mean) * 0.5 + 1e-6)

    hist_nsig = ROOT.TH1D('hist_nsig', 'hist_nsig', 100, nsig_lo, nsig_hi)
    hist_nsig.SetXTitle('#it{n}_{sig}')
    hist_nsig.SetYTitle('Pseudoexperiments')

    hist_ncom = ROOT.TH1D('hist_ncom', 'hist_ncom', 100, ncom_lo, ncom_hi)
    hist_ncom.SetXTitle('#it{n}_{com}')
    hist_ncom.SetYTitle('Pseudoexperiments')

    hist_exp = ROOT.TH1D('hist_exp', 'hist_exp', 100, exp_lo, exp_hi)
    hist_exp.SetXTitle('#it{e}_{0}')
    hist_exp.SetYTitle('Pseudoexperiments')

    hist_nsig_pull = ROOT.TH1D('hist_nsig_pull', 'hist_nsig_pull', 100, -10.0,
                               10.0)
    hist_nsig_pull.SetXTitle('#it{n}_{sig} pull')
    hist_nsig_pull.SetYTitle('Pseudoexperiments')

    hist_ncom_pull = ROOT.TH1D('hist_ncom_pull', 'hist_ncom_pull', 100, -10.0,
                               10.0)
    hist_ncom_pull.SetXTitle('#it{n}_{com} pull')
    hist_ncom_pull.SetYTitle('Pseudoexperiments')

    hist_exp_pull = ROOT.TH1D('hist_exp_pull', 'hist_exp_pull', 100, -10.0,
                              10.0)
    hist_exp_pull.SetXTitle('#it{e}_{0} pull')
    hist_exp_pull.SetYTitle('Pseudoexperiments')

    val_nsig = results['inputs']['Nsig']
    val_ncom = results['inputs']['Ncom']
    val_exp = results['inputs']['Expcom']

    all_toys = 0
    converged_toys = 0
    good_fit_toys = 0
    useful_seeds = list()
    bad_seeds = list()
    typical_seeds = list()

    # Loop over toys
    for toy in results['toys'].values():
        hist_status.Fill(toy['status'])
        hist_quality.Fill(toy['quality'])

        hist_nsig.Fill(toy['Nsig']['value'])
        hist_ncom.Fill(toy['Ncom']['value'])
        hist_exp.Fill(toy['Expcom']['value'])

        all_toys += 1

        good_toy = 0
        good_pull = 0

        if toy['status'] == 0 and toy['quality'] == 3:
            good_fit_toys += 1
        else:
            bad_seeds.append(toy['seed'])

        pulls = []

        if toy['Nsig']['error'] > 0:
            pull_nsig = (
                toy['Nsig']['value'] - val_nsig) / toy['Nsig']['error']
            hist_nsig_pull.Fill(pull_nsig)
            good_toy += 1
            pulls.append(pull_nsig)
            if abs(toy['Nsig']['value'] - val_nsig) < 0.001:
                good_pull += 1

        if toy['Ncom']['error'] > 0:
            pull_ncom = (
                toy['Ncom']['value'] - val_ncom) / toy['Ncom']['error']
            hist_ncom_pull.Fill(pull_ncom)
            good_toy += 1
            pulls.append(pull_ncom)
            if abs(toy['Ncom']['value'] - val_ncom) < 0.1:
                good_pull += 1

        if toy['Expcom']['error'] > 0:
            pull_exp = (
                toy['Expcom']['value'] - val_exp) / toy['Expcom']['error']
            hist_exp_pull.Fill(pull_exp)
            good_toy += 1
            pulls.append(pull_exp)
            if abs(toy['Expcom']['value'] - val_exp) < 0.001:
                good_pull += 1

        if good_toy == 3:
            converged_toys += 1
        if good_pull == 3:
            useful_seeds.append(toy['seed'])

        if (toy['status'] == 0 and toy['quality'] == 3 and len(pulls) == 3
                and all(abs(x) < 0.5 for x in pulls)):
            typical_seeds.append(toy['seed'])

    canvas_quality = ROOT.TCanvas('canvas_quality', 'canvas_quality', 1200,
                                  400)
    canvas_quality.Divide(2)

    canvas_quality.cd(1)
    hist_status.Draw()

    canvas_quality.cd(2)
    hist_quality.Draw()

    canvas_quality.SaveAs(f"{plot_name}_quality.pdf")

    canvas_yields = ROOT.TCanvas('canvas_yields', 'canvas_yields', 1800, 400)
    canvas_yields.Divide(3)

    canvas_yields.cd(1)
    hist_nsig.Draw()

    canvas_yields.cd(2)
    hist_ncom.Draw()

    canvas_yields.cd(3)
    hist_exp.Draw()

    canvas_yields.SaveAs(f"{plot_name}_yields.pdf")

    canvas_nsig = ROOT.TCanvas('canvas_nsig', 'canvas_nsig')
    hist_nsig_pull.Draw()
    fitHistogram(hist_nsig_pull)
    canvas_nsig.Update()

    canvas_nsig.SaveAs(f"{plot_name}_nsig.pdf")

    canvas_ncom = ROOT.TCanvas('canvas_ncom', 'canvas_ncom')
    hist_ncom_pull.Draw()
    fitHistogram(hist_ncom_pull)
    canvas_ncom.Update()

    canvas_ncom.SaveAs(f"{plot_name}_ncom.pdf")

    canvas_exp = ROOT.TCanvas('canvas_exp', 'canvas_exp')
    hist_exp_pull.Draw()
    fitHistogram(hist_exp_pull)
    canvas_exp.Update()

    canvas_exp.SaveAs(f"{plot_name}_exp.pdf")

    converged_fraction = converged_toys / all_toys

    print(f"\nFraction of converged toys (error > 0): {converged_fraction}")

    print(f"Fraction with status==0 and covQual==3: "
          f"{good_fit_toys / all_toys:.3f}")

    if val_nsig == 0:
        print("\nNOTE: nsig=0 is on the sig:yield >= 0 boundary, so the n_sig "
              "pull is truncated by construction.")

    print(f"\nUseful seeds = {useful_seeds}")
    print(f"\nTypical seeds (|pull| < 0.5, status==0, covQual==3, "
          f"{len(typical_seeds)} of {all_toys}) = {typical_seeds[:30]}" +
          (" ..." if len(typical_seeds) > 30 else ""))
    print(
        f"\nSeeds of NON-converged toys ({len(bad_seeds)}) = {bad_seeds[:30]}"
        + (" ..." if len(bad_seeds) > 30 else ""))


if __name__ == '__main__':

    # Set the plotting style
    import lhcbStyle
    lhcbStyle.applyStyle()
    import os

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--reparam', action='store_true', help='reparameterise pdf')
    parser.add_argument('--input', type=str, required=True, help='input file')
    args = parser.parse_args()

    file_suffix = os.path.basename(args.input).split(".j")[0]

    plotToyResults(
        input_file=args.input,
        plot_name=f"plots/{file_suffix}",
        reparam=args.reparam)
