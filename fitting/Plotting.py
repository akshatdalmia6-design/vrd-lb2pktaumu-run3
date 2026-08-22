import ROOT
import json
import ctypes
import math
import sys

# Set the plotting style
import lhcbStyle


def getDotsWithoutZero(dots):
    """
    For plots with residuals.
    """

    xValDot = ctypes.c_double(-1.E30)
    yValDot = ctypes.c_double(-1.E30)

    iBin = 0

    dataHist = ROOT.RooHist()

    while iBin < dots.GetN():
        dots.GetPoint(iBin, xValDot, yValDot)

        if yValDot.value > 0:

            dataHist.addBinWithXYError(
                xValDot.value,
                yValDot.value,
                # 0,
                # 0,
                dots.GetErrorXlow(iBin),
                dots.GetErrorXhigh(iBin),
                dots.GetErrorYlow(iBin),
                dots.GetErrorYhigh(iBin))

        iBin += 1

    dataHist.SetFillColor(ROOT.kWhite)
    dataHist.SetMarkerStyle(dots.GetMarkerStyle())
    dataHist.SetMarkerSize(dots.GetMarkerSize())
    dataHist.SetLineWidth(dots.GetLineWidth())

    return dataHist


def drawComponents(plot, pdf, components):
    '''
    Plotting tool for stacking components
    '''
    pdf.plotOn(plot, ROOT.RooFit.PrintEvalErrors(0),
               ROOT.RooFit.EvalErrorValue(0.))
    for i in range(len(components)):
        _component = ','.join(
            [components[j][0] for j in range(i, len(components))])
        print(_component)
        pdf.plotOn(plot, ROOT.RooFit.Components(_component),
                   ROOT.RooFit.FillColor(components[i][1]),
                   ROOT.RooFit.DrawOption('F'), ROOT.RooFit.PrintEvalErrors(0),
                   ROOT.RooFit.EvalErrorValue(0.))
    pdf.plotOn(plot, ROOT.RooFit.PrintEvalErrors(0),
               ROOT.RooFit.EvalErrorValue(0.))
    return


def drawComponentLines(plot, pdf, components, lineSettings=None, linewidth=3):
    '''
    Plotting tool for line components
    '''

    if lineSettings is None:
        lineSettings = [ROOT.kBlue, ROOT.kSolid]

    pdf.plotOn(plot, ROOT.RooFit.LineColor(lineSettings[0]),
               ROOT.RooFit.LineStyle(lineSettings[1]),
               ROOT.RooFit.LineWidth(linewidth),
               ROOT.RooFit.PrintEvalErrors(0), ROOT.RooFit.EvalErrorValue(0.),
               ROOT.RooFit.Name("model"))
    print("components=", components)
    for component, color, style, label in components:
        pdf.plotOn(plot, ROOT.RooFit.Components(component),
                   ROOT.RooFit.LineColor(color), ROOT.RooFit.LineStyle(style),
                   ROOT.RooFit.LineWidth(linewidth),
                   ROOT.RooFit.PrintEvalErrors(0),
                   ROOT.RooFit.EvalErrorValue(0.), ROOT.RooFit.Name(component))
    # pdf.plotOn(plot, ROOT.RooFit.PrintEvalErrors(0),
    #            ROOT.RooFit.EvalErrorValue(0.))
    return


def drawComponentLinesNorm(plot,
                           pdf,
                           components,
                           lineSettings=None,
                           linewidth=3):
    '''
    Plotting tool for line components
    '''

    NumEvent = ROOT.RooAbsReal.NumEvent

    if lineSettings is None:
        lineSettings = [ROOT.kBlue, "Total fit", ROOT.kSolid]

    pdf.plotOn(plot, ROOT.RooFit.LineColor(lineSettings[0]),
               ROOT.RooFit.LineStyle(lineSettings[1]),
               ROOT.RooFit.LineWidth(linewidth),
               ROOT.RooFit.PrintEvalErrors(0), ROOT.RooFit.EvalErrorValue(0.),
               ROOT.RooFit.Name("model"),
               ROOT.RooFit.Normalization(lineSettings[3], NumEvent))
    for component, color, style, legend, normalization in components.values():
        component.plotOn(plot, ROOT.RooFit.LineColor(color),
                         ROOT.RooFit.LineStyle(style),
                         ROOT.RooFit.LineWidth(linewidth),
                         ROOT.RooFit.PrintEvalErrors(0),
                         ROOT.RooFit.EvalErrorValue(0.),
                         ROOT.RooFit.Normalization(normalization, NumEvent),
                         ROOT.RooFit.Name(component.GetName()))
    # pdf.plotOn(plot, ROOT.RooFit.PrintEvalErrors(0),
    #            ROOT.RooFit.EvalErrorValue(0.))
    return


def setup1DRooPlot(plotTotal, rooFitVar, offset=1.1):

    dots_nozero = getDotsWithoutZero(plotTotal.getHist("dataSet0"))
    dots_nozero.SetName("dataSet")
    plotTotal.addPlotable(dots_nozero, "P")
    plotTotal.remove("dataSet0")

    binwidth = plotTotal.GetXaxis().GetBinWidth(1)
    binwidth = round(binwidth, -int(math.floor(math.log10(abs(binwidth)) - 2)))
    plotTotal.GetYaxis().SetTitle("Candidates per {0} {1}".format(
        binwidth, rooFitVar.getUnit()))
    plotTotal.GetYaxis().SetTitleOffset(offset)

    xlabel = str(rooFitVar.getTitle())
    if str(rooFitVar.getUnit()):
        xlabel += f" [{rooFitVar.getUnit()}]"

    plotTotal.GetXaxis().SetTitle(xlabel)
    plotTotal.GetXaxis().SetLabelOffset(0.018)

    return


def setupPullHist(plot_res, xlabel):
    '''
    Sets up pull histo. 
    '''
    plot_res.SetTitle("")
    plot_res.GetXaxis().SetTitle(xlabel)
    plot_res.GetXaxis().SetTitleSize(0.2)
    plot_res.GetXaxis().SetTitleOffset(0.9)
    plot_res.GetXaxis().SetTickSize(0.07)
    plot_res.GetYaxis().SetTickSize(0.024)
    plot_res.GetYaxis().SetTitle("#splitline{Normalised}{ Residuals}")
    plot_res.GetYaxis().SetTitleSize(0.16)
    plot_res.GetYaxis().SetTitleOffset(0.45)
    plot_res.GetXaxis().SetLabelSize(0.145)
    plot_res.GetYaxis().SetLabelSize(0.120)


regions = {"sigregion": "B_M > 5220 && B_M < 5340"}


def calculateNormFractions(mass, massmodel, region="sigregion"):
    '''
    Calculates normalisation for each component.
    '''

    model_int_allreg_1D = massmodel.createIntegral(
        ROOT.RooArgSet(mass), ROOT.RooFit.Range("fullregion"))
    model_int_reg_1D = massmodel.createIntegral(
        ROOT.RooArgSet(mass), ROOT.RooFit.Range(region))

    frac_model_B_M = model_int_reg_1D.getVal() / model_int_allreg_1D.getVal()

    print(f"frac_model_{region}_B_M = {frac_model_B_M}")

    return frac_model_B_M


def draw1Dprojection(rooFitVar,
                     rooFitModel,
                     dataSet,
                     nameOfPlot,
                     components=None,
                     lineSettings=None,
                     dotSettings=None,
                     nBins=82,
                     plotPulls=True,
                     logAxis=False,
                     yMin=None,
                     yMax=None,
                     varRange: str = None,
                     var_region=None,
                     plotLegend=None,
                     dataLabel="Data",
                     modelLabel="Total fit",
                     lhcbText="LHCb",
                     lumiText="7.5 fb^{#minus1}"):
    '''
    Plots a 1D fit projection.
    '''

    if lineSettings is None:
        lineSettings = [ROOT.kBlue, ROOT.kSolid]
    if dotSettings is None:
        dotSettings = [ROOT.RooAbsData.Poisson]
    if plotPulls:
        nameOfPlot += "WithPulls"

    options = [ROOT.RooFit.Name("dataSet0")]
    if varRange:
        print(f"Plot range: {varRange}")
        options += [ROOT.RooFit.NormRange(varRange)]

    plotTotal = dataSet.plotOn(
        rooFitVar.frame(nBins), ROOT.RooFit.DataError(dotSettings[0]), *options
        # ROOT.RooFit.XErrorSize(0) Doesn't work for angular distribution residuals somehow.
    )

    if var_region is None:
        drawComponentLines(
            plotTotal,
            rooFitModel,
            components=components,
            lineSettings=lineSettings)
    else:
        drawComponentLinesNorm(
            plotTotal,
            rooFitModel,
            components=components,
            lineSettings=lineSettings,
            linewidth=4)

    setup1DRooPlot(plotTotal, rooFitVar)

    canvas_size = [700, 640] if plotPulls else [700, 550]

    plotCanvas = ROOT.TCanvas(
        f"{rooFitVar.GetName}_{rooFitModel.GetName()}_{dataSet.GetName()}", "",
        *canvas_size)
    plotCanvas.SetTopMargin(0.07)
    plotCanvas.SetLeftMargin(0.25)
    plotCanvas.SetRightMargin(0.05)
    plotCanvas.SetBottomMargin(0)
    plotCanvas.Divide(1)
    plotCanvas.cd(1)

    # Wish to sep

    # For pull guiding lines
    xMin = plotTotal.GetXaxis().GetXmin()
    xMax = plotTotal.GetXaxis().GetXmax()
    guidelines = {
        0: [ROOT.kGray + 3, ROOT.kDashed],
        3: [ROOT.kRed, ROOT.kSolid],
        -3: [ROOT.kRed, ROOT.kSolid]
    }

    pad_fit_ylow = 0.277 if plotPulls else 0
    pad_fit = ROOT.TPad(
        f"p1_{rooFitVar.GetName}_{rooFitModel.GetName()}_{dataSet.GetName()}",
        "p1", 0, pad_fit_ylow, 1, 1, 0)

    if logAxis:
        plotTotal.SetMinimum(0.001)
        pad_fit.SetLogy()

    if yMin: plotTotal.SetMinimum(yMin)
    if yMax: plotTotal.SetMaximum(yMax)

    pad_fit.Draw()

    pad_fit.SetTopMargin(0.06)
    pad_fit.SetLeftMargin(0.15)
    pad_fit.SetBottomMargin(0.15)

    if plotPulls:
        pad_fit.SetBottomMargin(0.025)
        pad_res = ROOT.TPad(
            f"p2_{rooFitVar.GetName}_{rooFitModel.GetName()}_{dataSet.GetName()}",
            "p2", 0, 0, 1, 0.276, 0)
        pad_res.Draw()
        xlabel = plotTotal.GetXaxis().GetTitle()
        plotTotal.GetXaxis().SetTitle("")
        plotTotal.GetXaxis().SetLabelSize(0)

    pad_fit.cd()

    plotTotal.Draw("P")

    plotTotal.Print()

    keep = []

    if plotLegend is not None:

        n_entries = 2 + len(components or [])
        y_hi = 0.83
        y_lo = y_hi - 0.075 * n_entries  # legend height follows content
        x_lo, x_hi = (0.20, 0.57) if plotLegend == "left" else (0.57, 0.94)

        leg = ROOT.TLegend(x_lo, y_lo, x_hi, y_hi)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextFont(lhcbStyle.lhcbFont)
        leg.SetTextSize(lhcbStyle.lhcbTSize * 0.93)

        leg.AddEntry(plotTotal.findObject("dataSet"), dataLabel, 'lpe')
        leg.AddEntry(plotTotal.findObject("model"), modelLabel, 'l')
        for component, color, style, label in (components or []):
            leg.AddEntry(plotTotal.findObject(component), label, "l")
        leg.Draw("SAME")
        keep.append(leg)

        tag = lhcbStyle.lhcbLatex.Clone()
        tag.SetNDC(True)
        tag.SetTextAlign(13)  # left / top
        tag.DrawLatex(
            0.21, 0.92, f"#splitline{{{lhcbText}}}{{{lumiText}}}"
            if lumiText else lhcbText)
        keep.append(tag)
        '''

        if plotLegend == "left":
            leg = ROOT.TLegend(0.20, 0.48, 0.57, 0.83)
        else:
            leg = ROOT.TLegend(0.57, 0.48, 0.94, 0.83)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)

        leg.SetTextFont(lhcbStyle.lhcbFont)
        leg.SetTextSize(lhcbStyle.lhcbTSize * 0.93)

        leg.AddEntry(plotTotal.findObject("dataSet"), "Data", 'lpe')
        leg.AddEntry(plotTotal.findObject("model"), "Total fit", 'l')

        # for component, color, style, label, norm in components.values():
        for component, color, style, label in components:
            # leg.AddEntry(plotTotal.findObject(component.GetName()), label, "l")
            leg.AddEntry(plotTotal.findObject(component), label, "l")

        leg.Draw("SAME")

        lhcbStyle.lhcbNames[0].SetFillStyle(0)
        lhcbStyle.lhcbNames[0].SetX1(0.20)
        lhcbStyle.lhcbNames[0].SetY1(0.81)
        lhcbStyle.lhcbNames[0].SetX2(0.33)
        lhcbStyle.lhcbNames[0].SetY2(0.94)

        lhcbStyle.lhcbLatex.SetTextAlign(12)      
        lhcbStyle.lhcbLatex.SetTextSize(lhcbStyle.lhcbTSize * 0.85)

        lhcbStyle.lhcbLatex.SetText(0.21, 0.755, "7.5 fb^{#minus1}")
        lhcbStyle.lhcbNames[0].Draw("SAME")
        lhcbStyle.lhcbLatex.Draw("SAME")
        '''

    if plotLegend and "pimumu" in nameOfPlot:

        B_Mlabel = lhcbStyle.lhcbLatex.Clone()
        B_Mlabel.SetTextSize(0.056)
        B_Mlabeltext = "5260 < #it{m}(#it{#pi}^{#plus}#it{#mu}^{#plus}#it{#mu}^{#minus}) < 5300 MeV/#it{c}^{2}"
        B_Mlabel.SetText(0.285, 0.975, B_Mlabeltext)
        B_Mlabel.Draw("SAME")

        q2label = lhcbStyle.lhcbLatex.Clone()
        q2label.SetTextSize(0.055)
        q2labeltext = "  1.1 < #it{q}^{2} < 6.0 GeV^{2}/#it{c}^{4}" if "low" in nameOfPlot else "15.0 < #it{q}^{2} < 22.0 GeV^{2}/#it{c}^{4}"
        q2label.SetText(0.50, 0.872, q2labeltext)
        q2label.Draw("SAME")

    if plotPulls:

        pad_res.SetTopMargin(0.01)
        pad_res.SetLeftMargin(0.15)
        pad_res.SetBottomMargin(0.4)

        pad_res.cd()

        residuals = plotTotal.pullHist("dataSet", "model")

        plot_res = rooFitVar.frame()
        plot_res.addPlotable(residuals, "P")

        setupPullHist(plot_res, xlabel)

        plot_res.SetAxisRange(-5, 5, "Y")

        plot_res.Draw()

        for val, settings in guidelines.items():
            # if val >= residuals.GetYaxis().GetXmax(
            # ) or val <= residuals.GetYaxis().GetXmin():
            #     continue
            guidelines[val].append(ROOT.TLine(xMin, val, xMax, val))
            guidelines[val][-1].SetLineColor(settings[0])
            guidelines[val][-1].SetLineStyle(settings[1])
            guidelines[val][-1].Draw("SAME")

        residuals.Draw("P SAME")

    plotCanvas.SaveAs(f"{nameOfPlot}.pdf")

    if not plotPulls:

        plotCanvas.SaveAs(f"{nameOfPlot}.png")
        plotCanvas.SaveAs(f"{nameOfPlot}.eps")
        plotCanvas.SaveAs(f"{nameOfPlot}.C")

    plotCanvas = None
    pad_fit = None
    pad_res = None


def draw2Dprojection(
        rooFitVar1,
        rooFitVar2,
        rooFitModel,
        dataSet,
        nameOfPlot,
        components=None,
        lineSettings=None,
        dotSettings=None,
        nBins=None,
        plotPulls=True,
        logAxis=False,
        yMin=None,
        yMax=None,
        yMax2=None,
        plotLegend=False,
):
    '''
    Plots two 1D fit projections for a 2D model.
    '''

    if lineSettings is None:
        lineSettings = [ROOT.kBlue, ROOT.kSolid]
    if dotSettings is None:
        dotSettings = [ROOT.RooAbsData.Poisson]
    if plotPulls:
        nameOfPlot += "WithPulls"
    if nBins is None:
        nBins = [82, 40]

    linewidth = 1

    # First observable
    plotTotal1 = dataSet.plotOn(
        rooFitVar1.frame(nBins[0]),
        ROOT.RooFit.DataError(dotSettings[0]),
        ROOT.RooFit.LineWidth(linewidth),
        ROOT.RooFit.Name("dataSet0"),
        # ROOT.RooFit.XErrorSize(0) Doesn't work for angular distribution residuals somehow.
    )

    drawComponentLines(
        plotTotal1,
        rooFitModel,
        components=components,
        lineSettings=lineSettings,
        linewidth=linewidth + 1)

    setup1DRooPlot(plotTotal1, rooFitVar1, 1.2)

    # Second observable
    plotTotal2 = dataSet.plotOn(
        rooFitVar2.frame(nBins[1]),
        ROOT.RooFit.DataError(dotSettings[0]),
        ROOT.RooFit.LineWidth(linewidth),
        ROOT.RooFit.Name("dataSet0"),
        # ROOT.RooFit.XErrorSize(0) Doesn't work for angular distribution residuals somehow.
    )

    drawComponentLines(
        plotTotal2,
        rooFitModel,
        components=components,
        lineSettings=lineSettings,
        linewidth=linewidth + 1)

    setup1DRooPlot(plotTotal2, rooFitVar2, 1.2)

    # Now plotting
    canvas_size = [2 * 700, 640] if plotPulls else [2 * 700, 550]

    plotCanvas = ROOT.TCanvas(
        f"{rooFitVar1.GetName}_{rooFitVar2.GetName}_{rooFitModel.GetName()}_{dataSet.GetName()}",
        "", *canvas_size)
    plotCanvas.SetTopMargin(0.07)
    plotCanvas.SetLeftMargin(0.25)
    plotCanvas.SetRightMargin(0.05)
    plotCanvas.SetBottomMargin(0)
    plotCanvas.Divide(2)

    previous_HistLineWidth = ROOT.gStyle.GetHistLineWidth()
    previous_FrameLineWidth = ROOT.gStyle.GetFrameLineWidth()
    previous_LineWidth = ROOT.gStyle.GetLineWidth()
    previous_TickLength = ROOT.gStyle.GetTickLength("X")

    ROOT.gStyle.SetHistLineWidth(2)  #  Force histogram lines
    ROOT.gStyle.SetFrameLineWidth(1)  #  Thinner axis boxes
    ROOT.gStyle.SetLineWidth(1)  #  Thinner default lines
    ROOT.gStyle.SetTickLength(
        0.02, "x")  #  Make ticks shorter so they don't look huge

    # For pull guiding lines
    xMin1 = plotTotal1.GetXaxis().GetXmin()
    xMax1 = plotTotal1.GetXaxis().GetXmax()
    xMin2 = plotTotal2.GetXaxis().GetXmin()
    xMax2 = plotTotal2.GetXaxis().GetXmax()
    guidelines = {
        0: [ROOT.kGray + 3, ROOT.kDashed],
        3: [ROOT.kRed, ROOT.kSolid],
        -3: [ROOT.kRed, ROOT.kSolid]
    }

    pad_fit_ylow = 0.277 if plotPulls else 0

    plotCanvas.cd(1)

    pad_fit1 = ROOT.TPad(
        f"p1_{rooFitVar1.GetName}_{rooFitModel.GetName()}_{dataSet.GetName()}",
        "p1", 0, pad_fit_ylow, 1, 1, 0)

    if logAxis:
        plotTotal1.SetMinimum(0.001)
        pad_fit1.SetLogy()

    if yMin: plotTotal1.SetMinimum(yMin)
    if yMax: plotTotal1.SetMaximum(yMax)

    pad_fit1.Draw()

    pad_fit1.SetTopMargin(0.06)
    pad_fit1.SetLeftMargin(0.16)
    pad_fit1.SetBottomMargin(0.15)

    if plotPulls:
        pad_fit1.SetBottomMargin(0.025)
        pad_res1 = ROOT.TPad(
            f"p2_{rooFitVar1.GetName}_{rooFitModel.GetName()}_{dataSet.GetName()}",
            "p2", 0, 0, 1, 0.276, 0)
        pad_res1.Draw()
        xlabel = plotTotal1.GetXaxis().GetTitle()
        plotTotal1.GetXaxis().SetTitle("")
        plotTotal1.GetXaxis().SetLabelSize(0)

    pad_fit1.cd()

    plotTotal1.Draw("P")

    plotTotal1.Print()

    if plotLegend:

        legx1 = 0.45 if "jpsi" in nameOfPlot else 0.57

        leg = ROOT.TLegend(legx1, 0.48, 0.94, 0.83)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)

        leg.SetTextFont(lhcbStyle.lhcbFont)
        leg.SetTextSize(lhcbStyle.lhcbTSize * 0.93)

        leg.AddEntry(plotTotal1.findObject("dataSet"), "Data", 'lpe')
        leg.AddEntry(plotTotal1.findObject("model"), "Total fit", 'l')

        for component, color, style, label in components:
            leg.AddEntry(plotTotal1.findObject(component), label, "l")

        leg.Draw("SAME")

        lhcbStyle.lhcbNames[0].SetFillStyle(0)
        lhcbStyle.lhcbNames[0].SetX1(0.20)
        lhcbStyle.lhcbNames[0].SetY1(0.81)
        lhcbStyle.lhcbNames[0].SetX2(0.33)
        lhcbStyle.lhcbNames[0].SetY2(0.94)

        lhcbStyle.lhcbLatex.SetText(0.35, 0.88, "7.5 fb^{#minus1}")
        lhcbStyle.lhcbNames[0].Draw("SAME")
        lhcbStyle.lhcbLatex.Draw("SAME")

    if plotLegend and "pimumu" in nameOfPlot:

        q2label = lhcbStyle.lhcbLatex.Clone()
        q2label.SetTextSize(0.055)
        q2labeltext = "  1.1 < #it{q}^{2} < 6.0 GeV^{2}/#it{c}^{4}" if "low" in nameOfPlot else "15.0 < #it{q}^{2} < 22.0 GeV^{2}/#it{c}^{4}"
        q2label.SetText(0.50, 0.872, q2labeltext)
        q2label.Draw("SAME")

    if plotPulls:

        pad_res1.SetTopMargin(0.01)
        pad_res1.SetLeftMargin(0.16)
        pad_res1.SetBottomMargin(0.4)

        pad_res1.cd()

        residuals1 = plotTotal1.pullHist("dataSet", "model")

        plot_res1 = rooFitVar1.frame()
        plot_res1.addPlotable(residuals1, "P")

        setupPullHist(plot_res1, xlabel)

        plot_res1.SetAxisRange(-5, 5, "Y")

        plot_res1.Draw()

        for val, settings in guidelines.items():
            guidelines[val].append(ROOT.TLine(xMin1, val, xMax1, val))
            guidelines[val][-1].SetLineColor(settings[0])
            guidelines[val][-1].SetLineStyle(settings[1])
            guidelines[val][-1].Draw("SAME")

        residuals1.Draw("P SAME")

    plotCanvas.cd(2)

    pad_fit2 = ROOT.TPad(
        f"p1_{rooFitVar2.GetName}_{rooFitModel.GetName()}_{dataSet.GetName()}",
        "p1", 0, pad_fit_ylow, 1, 1, 0)

    if logAxis:
        plotTotal2.SetMinimum(0.001)
        pad_fit2.SetLogy()

    if yMin: plotTotal2.SetMinimum(yMin)
    if yMax2: plotTotal2.SetMaximum(yMax2)

    pad_fit2.Draw()

    pad_fit2.SetTopMargin(0.06)
    pad_fit2.SetLeftMargin(0.16)
    pad_fit2.SetBottomMargin(0.15)

    if plotPulls:
        pad_fit2.SetBottomMargin(0.025)
        pad_res2 = ROOT.TPad(
            f"p2_{rooFitVar2.GetName}_{rooFitModel.GetName()}_{dataSet.GetName()}",
            "p2", 0, 0, 1, 0.276, 0)
        pad_res2.Draw()
        xlabel = plotTotal2.GetXaxis().GetTitle()
        plotTotal2.GetXaxis().SetTitle("")
        plotTotal2.GetXaxis().SetLabelSize(0)

    pad_fit2.cd()

    plotTotal2.Draw("P")

    if plotLegend:
        lhcbStyle.lhcbNames[0].Draw("SAME")
        lhcbStyle.lhcbLatex.Draw("SAME")

    if plotLegend and "pimumu" in nameOfPlot:
        q2label.Draw("SAME")

    if plotPulls:

        pad_res2.SetTopMargin(0.01)
        pad_res2.SetLeftMargin(0.16)
        pad_res2.SetBottomMargin(0.4)

        pad_res2.cd()

        residuals2 = plotTotal2.pullHist("dataSet", "model")

        plot_res2 = rooFitVar2.frame()
        plot_res2.addPlotable(residuals2, "P")

        setupPullHist(plot_res2, xlabel)

        plot_res2.SetAxisRange(-5, 5, "Y")

        plot_res2.Draw()

        for val, settings in guidelines.items():
            guidelines[val].append(ROOT.TLine(xMin2, val, xMax2, val))
            guidelines[val][-1].SetLineColor(settings[0])
            guidelines[val][-1].SetLineStyle(settings[1])
            guidelines[val][-1].Draw("SAME")

        residuals2.Draw("P SAME")

    plotCanvas.Update()

    plotCanvas.SaveAs(f"{nameOfPlot}.pdf")

    if not plotPulls:

        plotCanvas.SaveAs(f"{nameOfPlot}.png")
        plotCanvas.SaveAs(f"{nameOfPlot}.eps")
        plotCanvas.SaveAs(f"{nameOfPlot}.C")

    plotCanvas = None
    pad_fit1 = None
    pad_res1 = None
    pad_fit2 = None
    pad_res2 = None
    ROOT.gStyle.SetHistLineWidth(previous_HistLineWidth)
    ROOT.gStyle.SetFrameLineWidth(previous_FrameLineWidth)
    ROOT.gStyle.SetLineWidth(previous_LineWidth)
    ROOT.gStyle.SetTickLength(
        previous_TickLength,
        "x")  #  Make ticks shorter so they don't look huge


def findIntervalFromDLL(histogram, value=0.5):
    '''
    Find interval from DLL distribution
    '''
    nmax = histogram.GetNbinsX()
    xmin = 0
    xmax = nmax

    for i in range(0, nmax):
        if xmin == 0 and histogram.GetBinContent(i + 1) < value:
            xmin = i + 1
    for i in range(0, nmax):
        if xmax == nmax and histogram.GetBinContent(nmax - i) < value:
            xmax = nmax - i

    print("xmin:", xmin, " xmax:", xmax)

    return (histogram.GetXaxis().GetBinLowEdge(xmin),
            histogram.GetXaxis().GetBinUpEdge(xmax))


def findIntervalFromFC(histogram, CL=0.683):
    '''
    Find interval from FC prob distribution
    '''
    nmax = histogram.GetNbinsX()
    xmin = 1
    xmax = nmax

    for i in range(0, nmax):
        if xmin == 1 and histogram.GetBinContent(i + 1) > 1 - CL:
            xmin = i + 1
    for i in range(0, nmax):
        if xmax == nmax and histogram.GetBinContent(nmax - i) > 1 - CL:
            xmax = nmax - i

    print("xmin:", xmin, " xmax:", xmax)

    return (histogram.GetXaxis().GetBinLowEdge(xmin),
            histogram.GetXaxis().GetBinUpEdge(xmax))
