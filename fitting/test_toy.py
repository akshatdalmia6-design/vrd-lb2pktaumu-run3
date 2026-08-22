from tabulate import tabulate

# Use jitting to load extra functionality
import ROOT
import json
import sys
import numpy as np

# Import the model
from LFV_model import LFVModel, compileBaseCode
from Plotting import drawComponents, draw1Dprojection
from test_fc import doFeldmanCousinsScan


def bootstrap_dataset(original_data, seed=1234):
    """
    Bootstraps a RooDataSet with a specific seed for reproducibility.
    """

    # 1. Set random number generator
    my_rng = ROOT.TRandom3(seed)

    dat_entries = original_data.numEntries()
    toy_entries = my_rng.Poisson(dat_entries)

    # 2. Set container
    bs_data = original_data.emptyClone(original_data.GetName() +
                                       f"_bootstrapped_{seed}")

    # 2. Resample
    for eventid in range(toy_entries):
        idx = my_rng.Integer(dat_entries)
        row = original_data.get(idx)
        bs_data.add(row)

    return bs_data


def generateToy(
        model,
        nsig,
        ncom,
        set_param,
        fix_param,
        seed,
        reparam=False,
        alt_model=None,
        alt_set_param=None,
        alt_fix_param=None,
        bootstrap=False,
        data=None,
        alt_cut=None,
):
    '''
    Generates a single toy and propagates results in dictionary.
    '''

    if bootstrap and data is None:
        print("Provide data for bootstrapping.")
        sys.exit(1)

    ROOT.RooRandom.randomGenerator().SetSeed(seed)

    model.setParam(set_param)
    model.fixParam(fix_param)

    # Generate a dataset
    if bootstrap:
        dataset = bootstrap_dataset(data, seed)
    else:
        dataset = model.generateDataset(nsig=nsig, ncom=ncom, seed=seed)

    print("Dataset for toy ")
    dataset.Print()

    pdf = model.getCombinedPdf()
    result = pdf.fitTo(dataset, ROOT.RooFit.Save(), ROOT.RooFit.PrintLevel(-1))

    if result is None:
        raise RuntimeError(f"toy seed {seed}: fit returned no result ")

    toy_results = {
        'Nsig': {
            'value': model.getValue('sig:yield'),
            'error': model.getError('sig:yield')
        },
        'Ncom': {
            'value': model.getValue('com:yield'),
            'error': model.getError('com:yield')
        },
        'Expcom': {
            'value': model.getValue('com:mass:e0'),
            'error': model.getError('com:mass:e0')
        },
        'seed': seed,
        'status': result.status(),
        'quality': result.covQual()
    }

    return toy_results, dataset


def generateToys(
        model,
        nsig,
        ncom,
        set_param,
        fix_param,
        output_name,
        reparam=False,
        ntoys=1,
        alt_model=None,
        alt_set_param=None,
        alt_fix_param=None,
        bootstrap=False,
        data=None,
        alt_cut=None,
):
    '''
    Generate several toys
    '''

    model.setParam(set_param)
    model.fixParam(fix_param)

    if 'com:mass:e0' not in set_param and 'com:mass:e0' not in fix_param:
        raise KeyError(
            "com:mass:e0 is in neither the set nor the fix parameter list. ")

    results = {}
    results['toys'] = {}
    results['inputs'] = {
        'Nsig': nsig,
        'Ncom': ncom,
        'Expcom': model.getValue('com:mass:e0'),
    }

    n_failed = 0

    for i in range(ntoys):
        print(f"Generating toy {i}", end='\r')
        try:
            results['toys']['toy%i' % i], toy_data = generateToy(
                model,
                nsig,
                ncom,
                set_param,
                fix_param,
                seed=i + 1,
                reparam=reparam,
                alt_model=alt_model,
                alt_set_param=alt_set_param,
                alt_fix_param=alt_fix_param,
                bootstrap=bootstrap,
                data=data,
                alt_cut=alt_cut,
            )
        except RuntimeError as exc:
            print(f"\ntoy {i} failed: {exc}")
            n_failed += 1

    print(f"\nGenerated {len(results['toys'])} toys ({n_failed} failed)")

    # Save the result
    output_file = f"{output_name}_toys.json"
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=4)

    return output_file


# Examples for pre-approval talk
# python test_toy.py --param parameters-low-qsq.json --nsig 93 --nmid 63 --ncom 196 --FH 0.0 --AFB 0.0 --plot
# python test_toy.py --param parameters-high-qsq.json --nsig 97 --nmid 25 --ncom 333 --FH 0.0 --AFB 0.0 --plot

if __name__ == '__main__':

    import os

    # Set the plotting style
    import lhcbStyle
    lhcbStyle.applyStyle()
    lhcbStyle.printLHCb()
    lhcbStyle.lhcbLatex.SetNDC()
    lhcbStyle.lhcbLatex.SetTextSize(lhcbStyle.lhcbTSize)

    # Use jitting to load extra functionality
    compileBaseCode()

    # Job options parsing
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--seed', type=int, default=1234, help='random number seed')
    parser.add_argument(
        '--nsig', type=float, default=116., help='signal yield')
    parser.add_argument(
        '--ncom', type=float, default=100., help='combinatorial yield')
    parser.add_argument(
        '--param', type=str, required=True, help='parameter file')
    parser.add_argument(
        '--ntoys', dest='ntoys', type=int, default=1, help='number of toys')
    parser.add_argument('--output', type=str, help='output file')
    parser.add_argument(
        '--plot', action='store_true', help='Plots fit projections.')
    parser.add_argument(
        '--silent', action='store_true', help='Silences RooFit.')
    parser.add_argument(
        '--mode',
        type=str,
        default="signal",
        choices=["signal", "norm"],
        help='Choose decay mode.')

    parser.add_argument(
        '--doFeldmanCousinsScan',
        action="store_true",
        help='Does Feldman-Cousins scan.')

    parser.add_argument(
        '--fcbins',
        type=int,
        nargs=3,
        default=[40, 0, 20],
        metavar=('NBINS', 'NSIGMIN', 'NSIGMAX'),
        help='FC scan binning: nbins, min, max')
    parser.add_argument(
        '--fctoys', type=int, default=1000, help='toys per FC scan point')

    args = parser.parse_args()

    ROOT.RooRandom.randomGenerator().SetSeed(args.seed)

    from utilities import silenceRooFit, createFolder

    if args.silent:
        silenceRooFit()

    plotDir = createFolder("./plots")

    binsData = 20

    file_suffix = os.path.basename(args.param).split(".")[0]
    output = f"test_toy_{file_suffix}_nsig{args.nsig}_ncom{args.ncom}"

    output = args.output if args.output else output

    # Set the model parameters
    from DefaultParameters import (getFixParameters, getSetParameters)

    set_param = getSetParameters(args.param)
    fix_param = getFixParameters(args.param)

    # Instantiate the model
    model = LFVModel(mode=args.mode)

    # Print model parameters
    model.setParam(set_param)
    model.fixParam(fix_param)
    print("\nParameters for toy generation:")
    model.printParam()

    _floating = [k for k in model.var if not model.var[k].isConstant()]
    print(f"\nFLOATING in the fit ({len(_floating)}): {_floating}")
    print(f"FIXED by --param   ({len(model.var) - len(_floating)} of "
          f"{len(model.var)} parameters)\n")
    if 'com:yield' not in _floating:
        print("WARNING: com:yield is FIXED.")

    test_toy_results, toy_data = generateToy(
        model=model,
        nsig=args.nsig,
        ncom=args.ncom,
        set_param=set_param,
        fix_param=fix_param,
        seed=args.seed)

    print(f"\nGenerated toy data sample:")
    toy_data.Print()

    # Print model parameters
    print("\nParameters after fit to single toy:")
    model.printParam()
    model.getCombinedPdf().Print()

    if args.plot:

        draw1Dprojection(
            model.mass,
            model.getCombinedPdf(),
            toy_data,
            nameOfPlot=f"{plotDir}/{output}",
            components=[[
                'sig:mass:C', ROOT.kGreen + 3, ROOT.kDashed, "Signal"
            ], ['com:mass', ROOT.kBlue - 7, 4, "Combinatorial"]],
            lineSettings=[ROOT.kBlue, ROOT.kSolid],
            dotSettings=[ROOT.RooAbsData.SumW2],
            nBins=binsData,
            plotPulls=True,
            logAxis=False,
            yMin=None,
            yMax=None,
            plotLegend="right")

    if args.doFeldmanCousinsScan:

        print("\nNOTE: running FC on a single toy dataset (seed "
              f"{args.seed}).")

        doFeldmanCousinsScan(
            data=toy_data,
            model=model,
            bins=list(args.fcbins),
            nFits=args.fctoys,
            plot_name=f"{plotDir}/{output}",
            output_name=f"./results/FC_scan_{output}",
            set_param=set_param)

    if args.ntoys > 1:
        toys_output = generateToys(
            model=model,
            nsig=args.nsig,
            ncom=args.ncom,
            set_param=set_param,
            fix_param=fix_param,
            output_name=f"./results/test_toys_{output}",
            ntoys=args.ntoys)

        from test_plot_toys import plotToyResults

        plotToyResults(
            input_file=toys_output,
            plot_name=f"{plotDir}/test_toys_{output}",
            # reparam=args.reparam
        )

#    if args.doProfileScan:
#
#        from test_profile import doProfileScan
#
#        doProfileScan(
#            data=dataset,
#            pdf=pdf,
#            FH=model.getFH(),
#            AFB=AFB,
#            bins=FH_AFB_bins,
#            plot_name=f"{plotDir}/{output}",
#            output_name=f"../results/profile_scan_{output}",
#            drawGuideLines=False,
#            yMax=5,
#            angle_var=args.angle_var)
