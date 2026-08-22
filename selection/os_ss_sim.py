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
'''
variables = {
    "VTXISO_1T_3p0_NParts": ["tau_VTXISO_tau_OneTrack_DChi2Vtx_3p0_NParts", 20, 0, 80, "N_{parts} (1-track, #chi^{2}_{vtx}<3.0)", ""],
    "VTXISO_1T_3p0_DCHI2": ["tau_VTXISO_tau_OneTrack_DChi2Vtx_3p0_Smallest_DELTACHI2", 80, -5, 100, "Smallest #Delta#chi^{2} (1-track, #chi^{2}_{vtx}<3.0)", ""],
    "VTXISO_1T_3p0_CHI2": ["tau_VTXISO_tau_OneTrack_DChi2Vtx_3p0_Smallest_CHI2", 80, -0.5, 3, "Smallest #chi^{2} (1-track, #chi^{2}_{vtx}<3.0)", ""],
   # "VTXISO_1T_3p0_MASS": ["tau_VTXISO_tau_OneTrack_DChi2Vtx_3p0_Smallest_DELTACHI2_MASS", 80, 3000, 11000, "Mass from smallest #Delta#chi^{2} (1-track, #chi^{2}_{vtx}<3.0)", "[MeV/#it{c}^{2}]"],

    "VTXISO_2T_3p0_NParts": ["tau_VTXISO_tau_TwoTracks_DChi2Vtx_3p0_NParts", 20, 0, 80, "N_{parts} (2-track, #chi^{2}_{vtx}<3.0)", ""],
    "VTXISO_2T_3p0_DCHI2": ["tau_VTXISO_tau_TwoTracks_DChi2Vtx_3p0_Smallest_DELTACHI2", 80, -5, 100, "Smallest #Delta#chi^{2} (2-track, #chi^{2}_{vtx}<3.0)", ""],
    "VTXISO_2T_3p0_CHI2": ["tau_VTXISO_tau_TwoTracks_DChi2Vtx_3p0_Smallest_CHI2", 80, -0.5, 3, "Smallest #chi^{2} (2-track, #chi^{2}_{vtx}<3.0)", ""],
   # "VTXISO_2T_3p0_MASS": ["tau_VTXISO_tau_TwoTracks_DChi2Vtx_3p0_Smallest_DELTACHI2_MASS", 80, 4000, 12000, "Mass from smallest #Delta#chi^{2} (2-track, #chi^{2}_{vtx}<3.0)", "[MeV/#it{c}^{2}]"],

    "VTXISO_1T_9p0_NParts": ["tau_VTXISO_tau_OneTrack_DChi2Vtx_9p0_NParts", 20, 0, 80, "N_{parts} (1-track, #chi^{2}_{vtx}<9.0)", ""],
    "VTXISO_1T_9p0_DCHI2": ["tau_VTXISO_tau_OneTrack_DChi2Vtx_9p0_Smallest_DELTACHI2", 80, -5, 100, "Smallest #Delta#chi^{2} (1-track, #chi^{2}_{vtx}<9.0)", ""],
    "VTXISO_1T_9p0_CHI2": ["tau_VTXISO_tau_OneTrack_DChi2Vtx_9p0_Smallest_CHI2", 80, -1, 9, "Smallest #chi^{2} (1-track, #chi^{2}_{vtx}<9.0)", ""],
   # "VTXISO_1T_9p0_MASS": ["tau_VTXISO_tau_OneTrack_DChi2Vtx_9p0_Smallest_DELTACHI2_MASS", 80, 4000, 15000, "Mass from smallest #Delta#chi^{2} (1-track, #chi^{2}_{vtx}<9.0)", "[MeV/#it{c}^{2}]"],

    "VTXISO_2T_9p0_NParts": ["tau_VTXISO_tau_TwoTracks_DChi2Vtx_9p0_NParts", 20, 0, 80, "N_{parts} (2-track, #chi^{2}_{vtx}<9.0)", ""],
    "VTXISO_2T_9p0_DCHI2": ["tau_VTXISO_tau_TwoTracks_DChi2Vtx_9p0_Smallest_DELTACHI2", 80, -5, 100, "Smallest #Delta#chi^{2} (2-track, #chi^{2}_{vtx}<9.0)", ""],
    "VTXISO_2T_9p0_CHI2": ["tau_VTXISO_tau_TwoTracks_DChi2Vtx_9p0_Smallest_CHI2", 80, -1, 9, "Smallest #chi^{2} (2-track, #chi^{2}_{vtx}<9.0)", ""],
   # "VTXISO_2T_9p0_MASS": ["tau_VTXISO_tau_TwoTracks_DChi2Vtx_9p0_Smallest_DELTACHI2_MASS", 80, 4000, 11000, "Mass from smallest #Delta#chi^{2} (2-track, #chi^{2}_{vtx}<9.0)", "[MeV/#it{c}^{2}]"],

    "VTXISO_1T_25p0_NParts": ["tau_VTXISO_tau_OneTrack_DChi2Vtx_25p0_NParts", 20, 0, 80, "N_{parts} (1-track, #chi^{2}_{vtx}<25.0)", ""],
    "VTXISO_1T_25p0_DCHI2": ["tau_VTXISO_tau_OneTrack_DChi2Vtx_25p0_Smallest_DELTACHI2", 80, -5, 100, "Smallest #Delta#chi^{2} (1-track, #chi^{2}_{vtx}<25.0)", ""],
    "VTXISO_1T_25p0_CHI2": ["tau_VTXISO_tau_OneTrack_DChi2Vtx_25p0_Smallest_CHI2", 80, -1, 25, "Smallest #chi^{2} (1-track, #chi^{2}_{vtx}<25.0)", ""],
   # "VTXISO_1T_25p0_MASS": ["tau_VTXISO_tau_OneTrack_DChi2Vtx_25p0_Smallest_DELTACHI2_MASS", 80, 4000, 12000, "Mass from smallest #Delta#chi^{2} (1-track, #chi^{2}_{vtx}<25.0)", "[MeV/#it{c}^{2}]"],

    "VTXISO_2T_25p0_NParts": ["tau_VTXISO_tau_TwoTracks_DChi2Vtx_25p0_NParts", 20, 0, 80, "N_{parts} (2-track, #chi^{2}_{vtx}<25.0)", ""],
    "VTXISO_2T_25p0_DCHI2": ["tau_VTXISO_tau_TwoTracks_DChi2Vtx_25p0_Smallest_DELTACHI2", 80, -5, 100, "Smallest #Delta#chi^{2} (2-track, #chi^{2}_{vtx}<25.0)", ""],
    "VTXISO_2T_25p0_CHI2": ["tau_VTXISO_tau_TwoTracks_DChi2Vtx_25p0_Smallest_CHI2", 80, -1, 25, "Smallest #chi^{2} (2-track, #chi^{2}_{vtx}<25.0)", ""],
   # "VTXISO_2T_25p0_MASS": ["tau_VTXISO_tau_TwoTracks_DChi2Vtx_25p0_Smallest_DELTACHI2_MASS", 80, 4000, 11000, "Mass from smallest #Delta#chi^{2} (2-track, #chi^{2}_{vtx}<25.0)", "[MeV/#it{c}^{2}]"],

    "VTXISO_1T_100p0_NParts": ["tau_VTXISO_tau_OneTrack_DChi2Vtx_100p0_NParts", 20, 0, 80, "N_{parts} (1-track, #chi^{2}_{vtx}<100.0)", ""],
    "VTXISO_1T_100p0_DCHI2": ["tau_VTXISO_tau_OneTrack_DChi2Vtx_100p0_Smallest_DELTACHI2", 80, -5, 100, "Smallest #Delta#chi^{2} (1-track, #chi^{2}_{vtx}<100.0)", ""],
    "VTXISO_1T_100p0_CHI2": ["tau_VTXISO_tau_OneTrack_DChi2Vtx_100p0_Smallest_CHI2", 80, -5, 25, "Smallest #chi^{2} (1-track, #chi^{2}_{vtx}<100.0)", ""],
   # "VTXISO_1T_100p0_MASS": ["tau_VTXISO_tau_OneTrack_DChi2Vtx_100p0_Smallest_DELTACHI2_MASS", 80, 3500, 13000, "Mass from smallest #Delta#chi^{2} (1-track, #chi^{2}_{vtx}<100.0)", "[MeV/#it{c}^{2}]"],

    "VTXISO_2T_100p0_NParts": ["tau_VTXISO_tau_TwoTracks_DChi2Vtx_100p0_NParts", 20, 0, 80, "N_{parts} (2-track, #chi^{2}_{vtx}<100.0)", ""],
    "VTXISO_2T_100p0_DCHI2": ["tau_VTXISO_tau_TwoTracks_DChi2Vtx_100p0_Smallest_DELTACHI2", 80, -5, 100, "Smallest #Delta#chi^{2} (2-track, #chi^{2}_{vtx}<100.0)", ""],
    "VTXISO_2T_100p0_CHI2": ["tau_VTXISO_tau_TwoTracks_DChi2Vtx_100p0_Smallest_CHI2", 80, -5, 25, "Smallest #chi^{2} (2-track, #chi^{2}_{vtx}<100.0)", ""],
   # "VTXISO_2T_100p0_MASS": ["tau_VTXISO_tau_TwoTracks_DChi2Vtx_100p0_Smallest_DELTACHI2_MASS", 80, 4000, 16000, "Mass from smallest #Delta#chi^{2} (2-track, #chi^{2}_{vtx}<100.0)", "[MeV/#it{c}^{2}]"],

    "VTXISO_1T_NO_NParts": ["tau_VTXISO_tau_OneTrack_DChi2Vtx_NO_NParts", 20, 0, 80, "N_{parts} (1-track, no cut)", ""],
    "VTXISO_1T_NO_DCHI2": ["tau_VTXISO_tau_OneTrack_DChi2Vtx_NO_Smallest_DELTACHI2", 80, -20, 200, "Smallest #Delta#chi^{2} (1-track, no cut)", ""],
    "VTXISO_1T_NO_CHI2": ["tau_VTXISO_tau_OneTrack_DChi2Vtx_NO_Smallest_CHI2", 80, -20, 150, "Smallest #chi^{2} (1-track, no cut)", ""],
   # "VTXISO_1T_NO_MASS": ["tau_VTXISO_tau_OneTrack_DChi2Vtx_NO_Smallest_DELTACHI2_MASS", 80, 3500, 13000, "Mass from smallest #Delta#chi^{2} (1-track, no cut)", "[MeV/#it{c}^{2}]"],

    "VTXISO_2T_NO_NParts": ["tau_VTXISO_tau_TwoTracks_DChi2Vtx_NO_NParts", 20, 0, 80, "N_{parts} (2-track, no cut)", ""],
    "VTXISO_2T_NO_DCHI2": ["tau_VTXISO_tau_TwoTracks_DChi2Vtx_NO_Smallest_DELTACHI2", 80, -20, 150, "Smallest #Delta#chi^{2} (2-track, no cut)", ""],
    "VTXISO_2T_NO_CHI2": ["tau_VTXISO_tau_TwoTracks_DChi2Vtx_NO_Smallest_CHI2", 80, -20, 150, "Smallest #chi^{2} (2-track, no cut)", ""],
   # "VTXISO_2T_NO_MASS": ["tau_VTXISO_tau_TwoTracks_DChi2Vtx_NO_Smallest_DELTACHI2_MASS", 80, 4000, 16000, "Mass from smallest #Delta#chi^{2} (2-track, no cut)", "[MeV/#it{c}^{2}]"],
}

'''

variables = {

    # ---------- DR2 < 0.05 ----------
    #4    "CONE_C_0p05_CMULT":  ["Lb_HEAD_CC_B_ChargedIso_DR2_0p05_CMULT", 40, 0, 40, "Charged cone mult. (#DeltaR^{2}<0.05)", ""],
    #4    "CONE_C_0p05_PTASY":  ["Lb_HEAD_CC_B_ChargedIso_DR2_0p05_PTASY", 80, -1.05, 1.05, "Charged cone #it{p}_{T} asym. (#DeltaR^{2}<0.05)", ""],
    #4    "CONE_C_0p05_PASY":   ["Lb_HEAD_CC_B_ChargedIso_DR2_0p05_PASY", 80, -1.05, 1.05, "Charged cone #it{p} asym. (#DeltaR^{2}<0.05)", ""],
    #4    "CONE_C_0p05_MAXP":   ["Lb_HEAD_CC_B_ChargedIso_DR2_0p05_Max_P", 80, 0, 60000, "Charged cone max #it{p} (#DeltaR^{2}<0.05)", "[MeV/#it{c}]"],
    #4    "CONE_C_0p05_MAXPT":  ["Lb_HEAD_CC_B_ChargedIso_DR2_0p05_Max_PT", 80, 0, 10000, "Charged cone max #it{p}_{T} (#DeltaR^{2}<0.05)", "[MeV/#it{c}]"],
    #4    "CONE_N_0p05_CMULT":  ["Lb_HEAD_NC_B_NeutralIso_DR2_0p05_CMULT", 40, 0, 40, "Neutral cone mult. (#DeltaR^{2}<0.05)", ""],
    #4    "CONE_N_0p05_PTASY":  ["Lb_HEAD_NC_B_NeutralIso_DR2_0p05_PTASY", 80, -1.05, 1.05, "Neutral cone #it{p}_{T} asym. (#DeltaR^{2}<0.05)", ""],
    #4    "CONE_N_0p05_PASY":   ["Lb_HEAD_NC_B_NeutralIso_DR2_0p05_PASY", 80, -1.05, 1.05, "Neutral cone #it{p} asym. (#DeltaR^{2}<0.05)", ""],
    #4    "CONE_N_0p05_MAXP":   ["Lb_HEAD_NC_B_NeutralIso_DR2_0p05_Max_P", 80, 0, 60000, "Neutral cone max #it{p} (#DeltaR^{2}<0.05)", "[MeV/#it{c}]"],
    #4    "CONE_N_0p05_MAXPT":  ["Lb_HEAD_NC_B_NeutralIso_DR2_0p05_Max_PT", 80, 0, 10000, "Neutral cone max #it{p}_{T} (#DeltaR^{2}<0.05)", "[MeV/#it{c}]"],

    # ---------- DR2 < 0.1 ----------
    #4    "CONE_C_0p1_CMULT":   ["Lb_HEAD_CC_B_ChargedIso_DR2_0p1_CMULT", 40, 0, 40, "Charged cone mult. (#DeltaR^{2}<0.1)", ""],
    #4    "CONE_C_0p1_PTASY":   ["Lb_HEAD_CC_B_ChargedIso_DR2_0p1_PTASY", 80, -1.05, 1.05, "Charged cone #it{p}_{T} asym. (#DeltaR^{2}<0.1)", ""],
    #4    "CONE_C_0p1_PASY":    ["Lb_HEAD_CC_B_ChargedIso_DR2_0p1_PASY", 80, -1.05, 1.05, "Charged cone #it{p} asym. (#DeltaR^{2}<0.1)", ""],
    #4    "CONE_C_0p1_MAXP":    ["Lb_HEAD_CC_B_ChargedIso_DR2_0p1_Max_P", 80, 0, 60000, "Charged cone max #it{p} (#DeltaR^{2}<0.1)", "[MeV/#it{c}]"],
    #4    "CONE_C_0p1_MAXPT":   ["Lb_HEAD_CC_B_ChargedIso_DR2_0p1_Max_PT", 80, 0, 10000, "Charged cone max #it{p}_{T} (#DeltaR^{2}<0.1)", "[MeV/#it{c}]"],
    #4    "CONE_N_0p1_CMULT":   ["Lb_HEAD_NC_B_NeutralIso_DR2_0p1_CMULT", 40, 0, 40, "Neutral cone mult. (#DeltaR^{2}<0.1)", ""],
    #4    "CONE_N_0p1_PTASY":   ["Lb_HEAD_NC_B_NeutralIso_DR2_0p1_PTASY", 80, -1.05, 1.05, "Neutral cone #it{p}_{T} asym. (#DeltaR^{2}<0.1)", ""],
    #4    "CONE_N_0p1_PASY":    ["Lb_HEAD_NC_B_NeutralIso_DR2_0p1_PASY", 80, -1.05, 1.05, "Neutral cone #it{p} asym. (#DeltaR^{2}<0.1)", ""],
    #4    "CONE_N_0p1_MAXP":    ["Lb_HEAD_NC_B_NeutralIso_DR2_0p1_Max_P", 80, 0, 60000, "Neutral cone max #it{p} (#DeltaR^{2}<0.1)", "[MeV/#it{c}]"],
    #4    "CONE_N_0p1_MAXPT":   ["Lb_HEAD_NC_B_NeutralIso_DR2_0p1_Max_PT", 80, 0, 10000, "Neutral cone max #it{p}_{T} (#DeltaR^{2}<0.1)", "[MeV/#it{c}]"],

    # ---------- DR2 < 0.15 ----------
    #4    "CONE_C_0p15_CMULT":  ["Lb_HEAD_CC_B_ChargedIso_DR2_0p15_CMULT", 40, 0, 40, "Charged cone mult. (#DeltaR^{2}<0.15)", ""],
    #4    "CONE_C_0p15_PTASY":  ["Lb_HEAD_CC_B_ChargedIso_DR2_0p15_PTASY", 80, -1.05, 1.05, "Charged cone #it{p}_{T} asym. (#DeltaR^{2}<0.15)", ""],
    #4    "CONE_C_0p15_PASY":   ["Lb_HEAD_CC_B_ChargedIso_DR2_0p15_PASY", 80, -1.05, 1.05, "Charged cone #it{p} asym. (#DeltaR^{2}<0.15)", ""],
    #4    "CONE_C_0p15_MAXP":   ["Lb_HEAD_CC_B_ChargedIso_DR2_0p15_Max_P", 80, 0, 60000, "Charged cone max #it{p} (#DeltaR^{2}<0.15)", "[MeV/#it{c}]"],
    #4    "CONE_C_0p15_MAXPT":  ["Lb_HEAD_CC_B_ChargedIso_DR2_0p15_Max_PT", 80, 0, 10000, "Charged cone max #it{p}_{T} (#DeltaR^{2}<0.15)", "[MeV/#it{c}]"],
    #4    "CONE_N_0p15_CMULT":  ["Lb_HEAD_NC_B_NeutralIso_DR2_0p15_CMULT", 40, 0, 40, "Neutral cone mult. (#DeltaR^{2}<0.15)", ""],
    #4    "CONE_N_0p15_PTASY":  ["Lb_HEAD_NC_B_NeutralIso_DR2_0p15_PTASY", 80, -1.05, 1.05, "Neutral cone #it{p}_{T} asym. (#DeltaR^{2}<0.15)", ""],
    #4    "CONE_N_0p15_PASY":   ["Lb_HEAD_NC_B_NeutralIso_DR2_0p15_PASY", 80, -1.05, 1.05, "Neutral cone #it{p} asym. (#DeltaR^{2}<0.15)", ""],
    #4    "CONE_N_0p15_MAXP":   ["Lb_HEAD_NC_B_NeutralIso_DR2_0p15_Max_P", 80, 0, 60000, "Neutral cone max #it{p} (#DeltaR^{2}<0.15)", "[MeV/#it{c}]"],
    #4    "CONE_N_0p15_MAXPT":  ["Lb_HEAD_NC_B_NeutralIso_DR2_0p15_Max_PT", 80, 0, 10000, "Neutral cone max #it{p}_{T} (#DeltaR^{2}<0.15)", "[MeV/#it{c}]"],

    # ---------- DR2 < 0.2 ----------
    #4    "CONE_C_0p2_CMULT":   ["Lb_HEAD_CC_B_ChargedIso_DR2_0p2_CMULT", 40, 0, 40, "Charged cone mult. (#DeltaR^{2}<0.2)", ""],
    #4    "CONE_C_0p2_PTASY":   ["Lb_HEAD_CC_B_ChargedIso_DR2_0p2_PTASY", 80, -1.05, 1.05, "Charged cone #it{p}_{T} asym. (#DeltaR^{2}<0.2)", ""],
    #4    "CONE_C_0p2_PASY":    ["Lb_HEAD_CC_B_ChargedIso_DR2_0p2_PASY", 80, -1.05, 1.05, "Charged cone #it{p} asym. (#DeltaR^{2}<0.2)", ""],
    #4    "CONE_C_0p2_MAXP":    ["Lb_HEAD_CC_B_ChargedIso_DR2_0p2_Max_P", 80, 0, 60000, "Charged cone max #it{p} (#DeltaR^{2}<0.2)", "[MeV/#it{c}]"],
    #4    "CONE_C_0p2_MAXPT":   ["Lb_HEAD_CC_B_ChargedIso_DR2_0p2_Max_PT", 80, 0, 10000, "Charged cone max #it{p}_{T} (#DeltaR^{2}<0.2)", "[MeV/#it{c}]"],
    #4    "CONE_N_0p2_CMULT":   ["Lb_HEAD_NC_B_NeutralIso_DR2_0p2_CMULT", 40, 0, 40, "Neutral cone mult. (#DeltaR^{2}<0.2)", ""],
    #4    "CONE_N_0p2_PTASY":   ["Lb_HEAD_NC_B_NeutralIso_DR2_0p2_PTASY", 80, -1.05, 1.05, "Neutral cone #it{p}_{T} asym. (#DeltaR^{2}<0.2)", ""],
    #4    "CONE_N_0p2_PASY":    ["Lb_HEAD_NC_B_NeutralIso_DR2_0p2_PASY", 80, -1.05, 1.05, "Neutral cone #it{p} asym. (#DeltaR^{2}<0.2)", ""],
    #4    "CONE_N_0p2_MAXP":    ["Lb_HEAD_NC_B_NeutralIso_DR2_0p2_Max_P", 80, 0, 60000, "Neutral cone max #it{p} (#DeltaR^{2}<0.2)", "[MeV/#it{c}]"],
    #4    "CONE_N_0p2_MAXPT":   ["Lb_HEAD_NC_B_NeutralIso_DR2_0p2_Max_PT", 80, 0, 10000, "Neutral cone max #it{p}_{T} (#DeltaR^{2}<0.2)", "[MeV/#it{c}]"],

    # ---------- DR2 < 0.5 ----------
    "CONE_C_0p5_CMULT": [
        "Lb_HEAD_CC_B_ChargedIso_DR2_0p5_CMULT", 40, 0, 30,
        "Charged cone mult. (#DeltaR^{2}<0.5)", ""
    ],
    "CONE_C_0p5_PTASY": [
        "Lb_HEAD_CC_B_ChargedIso_DR2_0p5_PTASY", 80, -0.6, 1.05,
        "Charged cone #it{p}_{T} asym. (#DeltaR^{2}<0.5)", ""
    ],
    "CONE_C_0p5_PASY": [
        "Lb_HEAD_CC_B_ChargedIso_DR2_0p5_PASY", 80, -0.4, 1.05,
        "Charged cone #it{p} asym. (#DeltaR^{2}<0.5)", ""
    ],
    "CONE_C_0p5_MAXP": [
        "Lb_HEAD_CC_B_ChargedIso_DR2_0p5_Max_P", 80, 0, 150000,
        "Charged cone max #it{p} (#DeltaR^{2}<0.5)", "[MeV/#it{c}]"
    ],
    "CONE_C_0p5_MAXPT": [
        "Lb_HEAD_CC_B_ChargedIso_DR2_0p5_Max_PT", 80, 0, 10000,
        "Charged cone max #it{p}_{T} (#DeltaR^{2}<0.5)", "[MeV/#it{c}]"
    ],
    "CONE_N_0p5_CMULT": [
        "Lb_HEAD_NC_B_NeutralIso_DR2_0p5_CMULT", 40, 0, 40,
        "Neutral cone mult. (#DeltaR^{2}<0.5)", ""
    ],
    "CONE_N_0p5_PTASY": [
        "Lb_HEAD_NC_B_NeutralIso_DR2_0p5_PTASY", 80, -0.6, 1.05,
        "Neutral cone #it{p}_{T} asym. (#DeltaR^{2}<0.5)", ""
    ],
    "CONE_N_0p5_PASY": [
        "Lb_HEAD_NC_B_NeutralIso_DR2_0p5_PASY", 80, -1.05, 1.05,
        "Neutral cone #it{p} asym. (#DeltaR^{2}<0.5)", ""
    ],
    "CONE_N_0p5_MAXP": [
        "Lb_HEAD_NC_B_NeutralIso_DR2_0p5_Max_P", 80, 0, 60000,
        "Neutral cone max #it{p} (#DeltaR^{2}<0.5)", "[MeV/#it{c}]"
    ],
    "CONE_N_0p5_MAXPT": [
        "Lb_HEAD_NC_B_NeutralIso_DR2_0p5_Max_PT", 80, 0, 10000,
        "Neutral cone max #it{p}_{T} (#DeltaR^{2}<0.5)", "[MeV/#it{c}]"
    ],

    # --PID--
    #4    "proton_pID_p": ["proton_PID_P",  80, -10, 250, "#Delta log #it{L(p)} (p)", " "],
    #4    "proton_pID_diff": ["proton_PID_P-proton_PID_K",  80, -50, 250, "#Delta log #it{L(p)} - #Delta log #it{L(K)} (p)", " "],
    #4    "kaon_pID_K": ["kaon_PID_K", 80, -20, 200, "#Delta log #it{L(K)} (K^{-})", " "],
    #4    "mu_pID_mu": ["mu_PID_MU", 80, -20, 40, "#Delta log #it{L(#mu)} (#mu)", " "],
    #4    "pi1_pID_K": ["pi1_PID_K", 80, -250, 200, "#Delta log #it{L(K)} (#pi 1)", " "],
    #4    "pi1_pID_p": ["pi1_PID_P", 80, -250, 200, "#Delta log #it{L(p)} (#pi 1)", " "],
    #4    "pi2_pID_K": ["pi2_PID_K", 80, -250, 200, "#Delta log #it{L(K)} (#pi 2)", " "],
    #4    "pi2_pID_p": ["pi2_PID_P", 80, -250, 200, "#Delta log #it{L(p)} (#pi 2)", " "],
    #4    "pi3_pID_K": ["pi3_PID_K", 80, -250, 200, "#Delta log #it{L(K)} (#pi 3)", " "],
    #4    "pi3_pID_p": ["pi3_PID_P", 80, -250, 200, "#Delta log #it{L(p)} (#pi 3)", " "]

    # --ProbNN--
    "proton_probNN_P": ["proton_PROBNN_P", 80, -0.1, 1.1, "ProbNN_P (p)", " "],
    "proton_probNN_K": ["proton_PROBNN_K", 80, -0.1, 1.1, "ProbNN_K (p)", " "],
    "proton_probNN_PI":
    ["proton_PROBNN_PI", 80, -0.1, 1.1, "ProbNN_PI (p)", " "],
    "kaon_probNN_K": ["kaon_PROBNN_K", 80, -0.1, 1.1, "ProbNN_K (K)", " "],
    "kaon_probNN_PI": ["kaon_PROBNN_PI", 80, -0.1, 1.1, "ProbNN_PI (K)", " "],
    "mu_probNN_MU": ["mu_PROBNN_MU", 80, -0.1, 1.1, "ProbNN_MU (#mu)", " "],
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

            plotfile_name = f"./plots/{sample}_bdt_conevars_sim_vs_{bkg_name}_1.pdf"
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

        radii = ["0p05", "0p1", "0p15", "0p2", "0p5"]

        groups = {
            **{
                f" Lb Cone radius DR2 < {r.replace('p', '.')}": (lambda v, r=r: f"_{r}_" in v)
                for r in radii
            }, " PID variables": lambda v: "pID" in v,
            " PROBNN variables": lambda v: "probNN" in v
        }

        for group_name, belongs in groups.items():
            print("\n" + "=" * 50)
            print(group_name)
            print("=" * 50)
            print(f"{'Variable':<22} {'sep vs SS':>10} {'sep vs OS sb':>12}")
            print("-" * 50)
            rows = [(v, s) for v, s in separations["ssdata"].items()
                    if belongs(v)]
            for var, sep_ss in sorted(rows, key=lambda x: x[1], reverse=True):
                sep_os = separations["osdata_sideband"][var]
                print(f"{var:<22} {sep_ss:10.4f} {sep_os:12.4f}")
            print("=" * 50)

    c1.Close()
