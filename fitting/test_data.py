from tabulate import tabulate

import ROOT
import json
import sys
import os
import numpy as np
from pathlib import Path

# Import the model
from LFV_model import LFVModel
from SampleAdmin import FitSamples
from Plotting import drawComponents, drawComponentLines, draw1Dprojection, draw2Dprojection, calculateNormFractions
from DefaultParameters import (
    getFixParameters, getSetParameters, getNuisanceParameters,
    getVetoHistogram, getEfficiencyParameters, getDefaultEfficiencyParameters,
    getDefaultVetoHistogram, saveParameters, print_correlation_matrix,
    print_covariance_matrix)
from utilities import silenceRooFit, createFolder

# from test_profile import doProfileScan

# Set the plotting style
import lhcbStyle
lhcbStyle.applyStyle()
lhcbStyle.printLHCb()
lhcbStyle.lhcbLatex.SetNDC()
lhcbStyle.lhcbLatex.SetTextSize(lhcbStyle.lhcbTSize)

# Job options parsing
import argparse

parser = argparse.ArgumentParser()
parser.add_argument(
    '--seed', type=int, default=1234, help='random number seed')
parser.add_argument(
    '--year',
    dest='year',
    type=str,
    default="2024",
    choices=["2016", "2017", "2018", "2024", "all"],
    help="Fit type. Choose year.")
parser.add_argument(
    '--polarity',
    dest='polarity',
    type=str,
    default="Both",
    choices=["Up", "Down", "Both"],
    help="Fit type. Choose polarity.")
parser.add_argument(
    '--mode',
    type=str,
    default="norm",
    choices=["signal", "norm"],
    help='Choose decay mode.')
parser.add_argument(
    '--PIDsel',
    dest='PIDsel',
    type=str,
    default="default",
    choices=["default", "tight"],
    help="Hadron PID selection. Choose between default and tight.")
parser.add_argument(
    '--bdtcut',
    dest='bdtcut',
    type=float,
    default=0.87,
    help='BDT threshold for the SIGNAL channel; candidates with '
    'bdt > bdtcut are kept (take this from the Punzi optimum). '
    'Ignored for the norm channel.')
parser.add_argument(
    '-fNBs',
    '--fitNBootstrappedSamples',
    dest='fitNBootstrappedSamples',
    type=int,
    default=0,
    help='Number of bootstrapped toy samples.')
parser.add_argument(
    '-pj',
    '--processes',
    dest='processes',
    type=int,
    default=10,
    help='Number of bootstrapped toy samples.')
parser.add_argument(
    '-tam',
    '--fitToysAltModel',
    dest="fitToysAltModel",
    default="default",
    type=str,
    choices=["default", "tightPIDsel"],
    help='Alternative model for toys.')
parser.add_argument(
    '--forcePosBkg',
    dest="com_pos",
    action="store_true",
    help='Forces positive definite combinatorial angular PDF.')
parser.add_argument('--fitMC', action="store_true", help='Fits mc.')
parser.add_argument(
    '--fitDoubleMidMC', action="store_true", help='Fits doubly mid mc.')
parser.add_argument('--fitOnlyMC', action="store_true", help='Fits mc.')
parser.add_argument(
    '--fitSideband', action="store_true", help='Does 1D fits of sideband.')
parser.add_argument(
    '--fitMass', action="store_true", help='Does 1D fit of mass.')
parser.add_argument(
    '--doProfileScan',
    action="store_true",
    help='Does profile Likelihood scan.')
parser.add_argument(
    '--reparam', action='store_true', help='Reparameterise pdf.')
parser.add_argument('--param', type=str, help='Parameter file.')
parser.add_argument('--silent', action='store_true', help='Silences RooFit.')
args = parser.parse_args()

ROOT.RooRandom.randomGenerator().SetSeed(args.seed)

CACHE_FILE = Path(
    " "
)


def _read_window():
    if not CACHE_FILE.exists():
        raise FileNotFoundError(
            f"{CACHE_FILE} not found -- run sb_fit.py first")
    d = np.load(CACHE_FILE)
    return float(d["window_lo"]), float(d["window_hi"])


if args.silent:
    silenceRooFit()

plotDir = createFolder("./plots")
paramDir = createFolder("./parameters")

canvas = ROOT.TCanvas('canvas', 'canvas', 600, 400)
binsMC = [82, 80]
binsData = [82, 80] if args.mode != "pimumu" else [20, 20]
binsSideData = [15, 20]

plotLineSettings = {
    'sig': [ROOT.kGreen + 3, ROOT.kDashed],
    'com': [ROOT.kBlue - 7, 4],  # ROOT.kDashDotted],
    'total': [ROOT.kBlue, ROOT.kSolid],
}

plotLineSettings['com'].append("Combinatorial")
plotLineSettings['total'].append("Total fit")
if args.mode == "norm":
    plotLineSettings['sig'].append(
        #        "#it{B}^{#plus} #rightarrow #it{J/#kern[-0.2]{#lower[0.1]{#psi}}}(#it{#mu}^{#plus}#it{#mu}^{#minus}) #it{K}^{#plus} "
        "#it{#Lambda}^{0}_{#it{b}} #rightarrow #it{J/#kern[-0.2]{#lower[0.1]{#psi}}} #it{K}^{#minus} #it{p}"
    )
else:
    plotLineSettings['sig'].append(
        "#it{#Lambda}^{0}_{#it{b}} #rightarrow #it{p}#it{K}^{#minus}#it{#tau}#it{#mu}"
    )
#  plotLineSettings['mid'].append(
#        "#it{B}^{#plus} #rightarrow #it{J/#kern[-0.2]{#lower[0.1]{#psi}}}(#it{#mu}^{#plus}#it{#mu}^{#minus}) #it{#pi}^{#plus}"
#     "#it{B}^{#plus} #rightarrow #it{J/#kern[-0.2]{#lower[0.1]{#psi}}} #it{#pi}^{#plus}"
#  )

massCompsToPlot = [['sig:mass', *plotLineSettings['sig']],
                   ['com:mass', *plotLineSettings['com']]]

totalCompsToPlot = [['sig:total', *plotLineSettings['sig']],
                    ['com:total', *plotLineSettings['com']]]

plotLineSettings['sig'][1] = 7
# plotLineSettings['com'][1] = 10
# plotLineSettings['mid'][1] = 2

file_suffix = f"test_fit_{args.mode}_{args.year}_{args.polarity}"  #_{args.q2bin}"

if args.PIDsel == "tight":

    file_suffix += f"_{args.PIDsel}PIDsel"

if args.com_pos:

    file_suffix += f"_posdefbkg"

default_param_file = f"parameters/{file_suffix}.json"

if args.mode == "signal":
    win_lo, win_hi = _read_window()
    MASS_LO, MASS_HI = 4800., 7000.
    SIG_LO, SIG_HI = win_lo, win_hi
else:
    MASS_LO, MASS_HI = 5400., 6000.
    SIG_LO, SIG_HI = MASS_LO, MASS_HI

# Instantiate the model
model = LFVModel(
    reparam=args.reparam,
    mode=args.mode,
    com_pos=args.com_pos,
    mass_range=(MASS_LO, MASS_HI))

# Set regions for studies in sections
#model.mass.setRange("sigregion", 5260, 5300)
model.mass.setRange("fullregion", MASS_LO, MASS_HI)
if args.mode == "signal":
    model.mass.setRange("sigwindow", SIG_LO, SIG_HI)
    model.mass.setRange("sb_low", MASS_LO, SIG_LO)
    model.mass.setRange("sb_high", SIG_HI, MASS_HI)
    model.mass.setRange("sideband", MASS_LO, MASS_HI)
else:
    model.mass.setRange("sideband", 5400, model.mass.getMax())

print("Setting up fit sample admin and fit samples")
samples = FitSamples(
    mode=args.mode,
    year=args.year,
    Polarity=args.polarity,
    fitVars=model.getFitVars(),
    setData=not args.fitOnlyMC,
    setMC=args.fitMC,
    PIDsel=args.PIDsel,
    bdt_cut=(args.bdtcut if args.mode == "signal" else None))

if args.fitMC:
    print("Doing fits to simulation")

    # Fix nuisance params
    model.fixParam({'mass:scaling': 1.0, 'mass:offset': 0.0})

    sigMC = samples.getMCSample(sample="sig")

    # Fit mass pdfs
    sig_mass_pdf = model.getSignalMassPdf()
    mc_fit_opts = [ROOT.RooFit.Extended(False), ROOT.RooFit.Save()]
    if sigMC.isWeighted():
        mc_fit_opts.append(ROOT.RooFit.SumW2Error(True))
    fit_sig_mass_pdf = sig_mass_pdf.fitTo(sigMC, *mc_fit_opts)
    sigMC.Print()
    model.mass.Print()

    draw1Dprojection(
        model.mass,
        sig_mass_pdf,
        sigMC,
        nameOfPlot=f"{plotDir}/{file_suffix}_sig_sim_mass",
        components=[],
        lineSettings=plotLineSettings['sig'],
        dotSettings=[ROOT.RooAbsData.SumW2],
        nBins=binsMC[0],
        plotPulls=True,
        logAxis=False,
        yMin=None,
        yMax=None,
        plotLegend="right",
        dataLabel="Simulation",
        modelLabel="Signal model",
        lhcbText="LHCb simulation",
        lumiText=None)

    fit_sig_mass_pdf.Print()

    # Fix mass parameters
    model.fixParams(keyword="mass")

    # Fix mass shift and scaling from JpsiK for other modes
    if args.mode in ["signal"]:
        norm_param_file = (
            f"parameters/test_fit_norm_{args.year}_{args.polarity}.json")
        if not os.path.exists(norm_param_file):
            sys.exit(f"[signal] missing {norm_param_file}; run "
                     f"`python test_data.py --mode norm --fitMC` first.")
        nuisance_params = getNuisanceParameters(filename=norm_param_file)
        print(nuisance_params)
        model.fixParam(nuisance_params, fixed=True)
    else:
        model.fixParam({
            'mass:scaling': 1.0,
            'mass:offset': 0.0,
        },
                       fixed=False)
    # Release nuisance
    model.fixParams(keyword='com:mass:e0', fixed=False)

if args.fitOnlyMC:
    sys.exit(0)

data = samples.getDataSample()

if not args.fitMC:

    inputParamFile = args.param if args.param else default_param_file

    print(f"Setting the model parameters using file: {inputParamFile}")

    set_param = getSetParameters(inputParamFile)
    fix_param = getFixParameters(inputParamFile)
    model.setParam(set_param)
    model.fixParam(fix_param)

# Fix first order coefficient for the combinatorial background for the jpsik mode.
if args.mode == "norm":

    model.fixParam({
        'com:angle:p1': 0.0,
    }, fixed=True)

if args.fitSideband:

    # model.fixParam({
    #         "com:angle:p1": 0.0,
    #     })
    if args.mode == "signal":
        sb_cut = f"(Lb_M < {SIG_LO} || Lb_M > {SIG_HI})"
        sb_range = "sb_low,sb_high"
    else:
        sb_cut = "Lb_M > 5400"
        sb_range = "sideband"

    sideband_data = data.reduce(ROOT.RooFit.Cut(sb_cut))
    sideband_data.SetName(data.GetName() + "_sideband")

    # Fit combinatorial mass pdf
    com_mass_pdf = model.getComMassPdf()
    fit_com_mass_pdf = com_mass_pdf.fitTo(sideband_data,
                                          ROOT.RooFit.Extended(False),
                                          ROOT.RooFit.Range(sb_range),
                                          ROOT.RooFit.Save())

    draw1Dprojection(
        model.mass,
        com_mass_pdf,
        sideband_data,
        nameOfPlot=f"{plotDir}/{file_suffix}_sideband_mass",
        components=[],
        lineSettings=plotLineSettings['com'],
        dotSettings=[ROOT.RooAbsData.Poisson],
        nBins=binsSideData[0],
        plotPulls=True,
        logAxis=False,
        yMin=None,
        yMax=None,
        varRange=sb_range,
        plotLegend="right",
        dataLabel="Data (sidebands)",
        modelLabel="Combinatorial")

    fit_com_mass_pdf.Print()

    comb_exp = sideband_data.numEntries() / calculateNormFractions(
        model.mass, model.pdf[f"com:mass"], sb_range)
    print("Events in sideband = ", sideband_data.numEntries())
    print("Expected comb_exp in full range = ", comb_exp)

    model_data = model.printParam()

    sig_veto = (None, None)
    mid_veto = (None, None)

    saveParameters(
        data=model_data,
        sig_veto=sig_veto,
        mid_veto=mid_veto,
        eff_param_path="",
        output_file=f"parameters/{file_suffix}.json")

    sys.exit()

if args.fitMass:
    ##### 1D mass fit to data #######
    print("Performing 1D mass fit to data.")

    total_mass_pdf = model.getCombinedMassPdf()

    total_mass_pdf_fit = total_mass_pdf.fitTo(data, ROOT.RooFit.Extended(True),
                                              ROOT.RooFit.Save())

    yMax = 150 if args.mode == "pimumu" else None

    draw1Dprojection(
        model.mass,
        total_mass_pdf,
        data,
        nameOfPlot=f"{plotDir}/{file_suffix}_data_mass",
        components=massCompsToPlot,
        dotSettings=[ROOT.RooAbsData.Poisson],
        nBins=binsData[0],
        plotPulls=True,
        logAxis=False,
        yMin=None,
        yMax=yMax)

    data.Print()
    total_mass_pdf.Print()
    total_mass_pdf_fit.Print()

    # comb_exp = data.sumEntries("Lb_M > 5700") / calculateNormFractions(
    #    model.mass, model.pdf[f"com:mass"], "sideband")
    # print("Events in sideband = ", data.sumEntries("Lb_M > 5600"))
    # print("Expected comb_exp in full range = ", comb_exp)
    '''

    # For JpsiK release nuisance again
    if args.mode == "sig":
        model.fixParams(keyword='mass:scaling', fixed=False)
        model.fixParams(keyword='mass:offset', fixed=False)
    '''

# Set the veto histograms for rare mode, not needed for control modes
sig_veto = (None, None)
mid_veto = (None, None)

model.printParam()

# sys.exit()

# Fit total pdf to data
total_pdf = model.getCombinedPdf()

total_pdf.Print()

total_pdf_fit = total_pdf.fitTo(data, ROOT.RooFit.Extended(True),
                                ROOT.RooFit.Save())
massComps = {}
# for i_comp in range(len(massCompsToPlot)):
#    comp_name = massCompsToPlot[i_comp][0].split(":")[0]

#    i_comp_norm = model.var[f"{comp_name}:yield"].getVal(
#    ) * calculateNormFractions(model.mass, model.pdf[f"{comp_name}:mass"])

#    print(comp_name, i_comp_norm, model.var[f"{comp_name}:yield"].getVal())

#    model.var[f"{comp_name}:yield"].setVal(i_comp_norm)

#    massCompsToPlot[i_comp].append(i_comp_norm)
#    massComps[comp_name] = [
#            model.pdf[f"{comp_name}:mass"], *massCompsToPlot[i_comp][1:]
#        ]

draw1Dprojection(
    model.mass,
    total_pdf,
    data,
    nameOfPlot=f"{plotDir}/{file_suffix}_data_total",
    components=massCompsToPlot,
    lineSettings=plotLineSettings['total'],
    dotSettings=[ROOT.RooAbsData.SumW2],
    nBins=binsMC[0],
    plotPulls=True,
    logAxis=False,
    yMin=None,
    yMax=None,
    plotLegend="right")
'''
yMax = 160
yMax2 = 50
if args.mode == "jpsik":
    yMax = 450000
    yMax2 = 40000
elif args.mode == "jpsipi":
    yMax = 18000
    yMax2 = 3000
'''

# Print fit results
data.Print()
total_pdf_fit.Print()

# Print model parameters
model_data = model.printParam()

# Covariances and correlations
print_covariance_matrix(total_pdf_fit)
print_correlation_matrix(total_pdf_fit)

if args.fitMC:
    saveParameters(
        data=model_data,
        sig_veto=sig_veto,
        mid_veto=mid_veto,
        eff_param_path="",
        output_file=f"parameters/{file_suffix}.json")

sys.exit(0)
