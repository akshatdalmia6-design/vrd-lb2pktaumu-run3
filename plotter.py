import ROOT
import sys
sys.path.append("../utils")
from lhcbPlotStyle import setLHCbPlotStyle
setLHCbPlotStyle()
input_dir = "/afs/cern.ch/work/f/fabudine/public/vrd-lb2pktaumu/tuples"


def make_plots(sig_path, bkg_path, vars_plot, n_bins=50):
    filter_code = "tau_M>800 && tau_M<1500 && Lb_DIRA_OWNPV>0.9998 && tau_ENDVERTEX_Z-Lb_ENDVERTEX_Z>-1 && tau_ENDVERTEX_Z-Lc_ENDVERTEX_Z > -1 && Lc_M > 2300 && proton_IPCHI2_OWNPV > 25 && kaon_IPCHI2_OWNPV > 25 && mu_IPCHI2_OWNPV > 25 && pi1_IPCHI2_OWNPV > 25 && pi2_IPCHI2_OWNPV > 25 && pi3_IPCHI2_OWNPV > 25 && (mu_L0MuonDecision_TOS || Lb_L0HadronDecision_TOS || Lb_L0HadronDecision_TIS) && (Lb_Hlt1TrackMVADecision_TOS || Lb_Hlt1TwoTrackMVADecision_TOS || Lb_Hlt1TrackMuonMVADecision_TOS) && (Lb_Hlt2Topo2BodyDecision_TOS || Lb_Hlt2Topo3BodyDecision_TOS || Lb_Hlt2Topo4BodyDecision_TOS || Lb_Hlt2TopoMu2BodyDecision_TOS || Lb_Hlt2TopoMu3BodyDecision_TOS || Lb_Hlt2TopoMu4BodyDecision_TOS) && Lc_FDCHI2_OWNPV > 25"

    df_sig = ROOT.RDataFrame(
        "LbTuple1/DecayTree",
        sig_path).Filter(filter_code + " && Lb_BKGCAT == 50")
    df_bkg = ROOT.RDataFrame("LbTuple1/DecayTree",
                             bkg_path).Filter(filter_code)

    histos = {}
    for var, title in vars_plot:
        h_model = (f"{var}", f";{title};", n_bins, 0, 0)
        histos[f"{var}_s"] = df_sig.Histo1D(h_model, var)
        histos[f"{var}_b"] = df_bkg.Histo1D(h_model, var)

    for var, title in vars_plot:
        c = ROOT.TCanvas(f"{var}", "", 800, 700)
        h_s = histos[f"{var}_s"].GetValue()
        h_b = histos[f"{var}_b"].GetValue()

        h_s.Scale(1.0 / h_s.Integral())
        h_b.Scale(1.0 / h_b.Integral())

        h_s.SetLineColor(ROOT.kBlue)
        h_b.SetLineColor(ROOT.kRed)

        h_s.Draw("HIST")
        h_b.Draw("HIST SAME")

        c.SaveAs(f"plot_{var}.pdf")

    print(f"Final Signal Entries:     {df_sig.Count().GetValue()}")
    print(f"Final Background Entries: {df_bkg.Count().GetValue()}")


if __name__ == '__main__':
    sig_path = f"{input_dir}/mcLb_pKtaumu_3pi_2018_MagUp_mass1.root"
    bkg_path = f"{input_dir}/Lb_pKtaumuSS_3pi_2018_MagUp_filtered_mass1.root"
    vars_plot = [("tau_M", "tau mass_2"), ("Lb_M", "Lb mass"),
                 ("Lc_M", "Lc mass [MeV]"),
                 ("proton_IPCHI2_OWNPV", "Proton IP #chi^{2}"),
                 ("kaon_IPCHI2_OWNPV", "Kaon IP #chi^{2}"),
                 ("mu_IPCHI2_OWNPV", "Muon IP #chi^{2}"),
                 ("Lb_DIRA_OWNPV", "Lb DIRA"),
                 ("Lc_FDCHI2_OWNPV", "Lc FD #chi^{2}")]

    make_plots(sig_path, bkg_path, vars_plot)
