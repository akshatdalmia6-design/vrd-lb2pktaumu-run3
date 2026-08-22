from ROOT import RDataFrame, gInterpreter, EnableImplicitMT
import argparse
import os

# Enable multi-threading for faster processing (optional, but highly recommended)
EnableImplicitMT()

# JIT compile the helper functions in C++ to handle the Lorentz vector math at native speed
gInterpreter.Declare("""
#include <Math/Vector4D.h>
#include <cmath>

namespace MassCalcs {
    // Standard particle mass constants (in MeV)
    const double M_p   = 938.272;
    const double M_K   = 493.677;
    const double M_pi  = 139.57;
    const double M_e   = 0.511;
    const double M_tau = 1776.86;

    using LorentzVector = ROOT::Math::PtEtaPhiMVector;

    // Helper to build a LorentzVector from ROOT branches (assuming PX, PY, PZ exist)
    inline ROOT::Math::PxPyPzMVector makeVec(double px, double py, double pz, double m) {
        return ROOT::Math::PxPyPzMVector(px, py, pz, m);
    }

    // Example calculation for invariant mass of N particles
    // We use PxPyPzMVector as it matches your original PX, PY, PZ, M input structure
    double invMass2(double px1, double py1, double pz1, double m1,
                    double px2, double py2, double pz2, double m2) {
        auto v1 = makeVec(px1, py1, pz1, m1);
        auto v2 = makeVec(px2, py2, pz2, m2);
        return (v1 + v2).M();
    }

    double invMass3(double px1, double py1, double pz1, double m1,
                    double px2, double py2, double pz2, double m2,
                    double px3, double py3, double pz3, double m3) {
        auto v1 = makeVec(px1, py1, pz1, m1);
        auto v2 = makeVec(px2, py2, pz2, m2);
        auto v3 = makeVec(px3, py3, pz3, m3);
        return (v1 + v2 + v3).M();
    }

    double invMass4(double px1, double py1, double pz1, double m1,
                    double px2, double py2, double pz2, double m2,
                    double px3, double py3, double pz3, double m3,
                    double px4, double py4, double pz4, double m4) {
        auto v1 = makeVec(px1, py1, pz1, m1);
        auto v2 = makeVec(px2, py2, pz2, m2);
        auto v3 = makeVec(px3, py3, pz3, m3);
        auto v4 = makeVec(px4, py4, pz4, m4);
        return (v1 + v2 + v3 + v4).M();
    }

    double invMass5(double px1, double py1, double pz1, double m1,
                    double px2, double py2, double pz2, double m2,
                    double px3, double py3, double pz3, double m3,
                    double px4, double py4, double pz4, double m4,
                    double px5, double py5, double pz5, double m5) {
        auto v1 = makeVec(px1, py1, pz1, m1);
        auto v2 = makeVec(px2, py2, pz2, m2);
        auto v3 = makeVec(px3, py3, pz3, m3);
        auto v4 = makeVec(px4, py4, pz4, m4);
        auto v5 = makeVec(px5, py5, pz5, m5);
        return (v1 + v2 + v3 + v4 + v5).M();
    }

    double invMass6(double px1, double py1, double pz1, double m1,
                    double px2, double py2, double pz2, double m2,
                    double px3, double py3, double pz3, double m3,
                    double px4, double py4, double pz4, double m4,
                    double px5, double py5, double pz5, double m5,
                    double px6, double py6, double pz6, double m6) {
        auto v1 = makeVec(px1, py1, pz1, m1);
        auto v2 = makeVec(px2, py2, pz2, m2);
        auto v3 = makeVec(px3, py3, pz3, m3);
        auto v4 = makeVec(px4, py4, pz4, m4);
        auto v5 = makeVec(px5, py5, pz5, m5);
        auto v6 = makeVec(px6, py6, pz6, m6);
        return (v1 + v2 + v3 + v4 + v5 + v6).M();
    }

    // For squared invariant masses (scaled by 1e-6)
    double invMassSqScaled2(double px1, double py1, double pz1, double m1,
                            double px2, double py2, double pz2, double m2) {
        auto v1 = makeVec(px1, py1, pz1, m1);
        auto v2 = makeVec(px2, py2, pz2, m2);
        return (v1 + v2).M2() * 1e-6;
    }
}
""")


def mass_vetoes(lepton="mu"):
    PDG_D0_M = 1864.84
    PDG_Dp_M = 1869.66
    PDG_Kstz_M = 895.55
    PDG_phi_M = 1019.46
    PDG_Lb_M = 5619.60

    # --- 1. OPPOSITE-SIGN VETOES (not is_PplusEplus) ---
    d0_k3pi = f"abs(k3pi_M-{PDG_D0_M})>30"
    d0_kpi_os = f"abs(kpi1_M-{PDG_D0_M})>30 && abs(kpi3_M-{PDG_D0_M})>30"
    dplus_kpipi = f"abs(kpi1pi3_M-{PDG_Dp_M})>30"
    lcstar_2625 = "pk3pi_M-pkpi1_M>350 && pk3pi_M-pkpi3_M>350"
    lc_pk3pi = "pk3pi_M>2310"
    lc_pkpi_os = "pkpi1_M>2320 && pkpi3_M>2320"
    kst_kpi_os = f"abs(kpi1_M-{PDG_Kstz_M})>60 && abs(kpi3_M-{PDG_Kstz_M})>60"

    Vetoes_OS = f"({d0_k3pi} && {d0_kpi_os} && {dplus_kpipi} && {lcstar_2625} && {lc_pk3pi} && {lc_pkpi_os} && {kst_kpi_os})"

    # --- 2. SAME-SIGN VETOES (is_PplusEplus) ---
    d0_kpi_ss = f"abs(kpi2_M-{PDG_D0_M})>30"
    dstar_d0_pi_kmunu = "kepi2_M-ke_M>155"
    dstar_d0_pi_kpi_misid = "kepi2_kpipi_M-ke_kpi_M>155"
    lc_pke = "pke_M>2320"
    lc_pkpi_ss = "pkpi2_M>2320"
    kst_kpi_ss = f"abs(kpi2_M-{PDG_Kstz_M})>60"
    d0_k3pi_misd = f"abs(p3pi_k3pi_M-{PDG_D0_M})>30"

    Vetoes_SS = f"({d0_kpi_ss} && {dstar_d0_pi_kmunu} && {dstar_d0_pi_kpi_misid} && {lc_pke} && {lc_pkpi_ss} && {kst_kpi_ss} && {d0_k3pi_misd})"

    # --- 3. COMMON VETOES ---
    phi = f"abs(pk_kk_M-{PDG_phi_M})>12"
    lb_pk3pie_misd = f"abs(pk3pie_pk3pipi_M-{PDG_Lb_M})>60"
    lb_pk3pie = f"abs(pk3pie_M-{PDG_Lb_M})>60"
    common_vetoes = f"({phi} && {lb_pk3pie_misd} && {lb_pk3pie})"

    # --- 4. DYNAMIC SELECTION ---
    # We dynamically select the right veto group depending on the charge multiplication
    # (proton_CHARGE * lepton_CHARGE) > 0 means same-sign, else opposite-sign
    charge_condition = f"(proton_CHARGE * {lepton}_CHARGE) > 0"
    dynamic_vetoes = f"({charge_condition} ? {Vetoes_SS} : {Vetoes_OS})"

    return f"({dynamic_vetoes} && {common_vetoes})"


def run_veto_pipeline(df=None, lepton="mu", pi1="pi1", pi2="pi2", pi3="pi3"):
    # Setup inputs

    # Particle prefix mappings to match loop properties
    p = "proton"
    k = "kaon"
    tau = "tau"

    # Shortcuts for passing branch coordinates to C++ helper
    def xyz(name):
        return f"{name}_PX, {name}_PY, {name}_PZ"

    # Mass aliases
    m_p, m_k, m_pi, m_e, m_tau = "MassCalcs::M_p", "MassCalcs::M_K", "MassCalcs::M_pi", "MassCalcs::M_e", "MassCalcs::M_tau"

    # Define standard masses
    df = df.Define("pkpi1_M", f"MassCalcs::invMass3({xyz(p)}, {m_p}, {xyz(k)}, {m_k}, {xyz(pi1)}, {m_pi})") \
           .Define("pkpi2_M", f"MassCalcs::invMass3({xyz(p)}, {m_p}, {xyz(k)}, {m_k}, {xyz(pi2)}, {m_pi})") \
           .Define("pkpi3_M", f"MassCalcs::invMass3({xyz(p)}, {m_p}, {xyz(k)}, {m_k}, {xyz(pi3)}, {m_pi})") \
           .Define("pke_M",   f"MassCalcs::invMass3({xyz(p)}, {m_p}, {xyz(k)}, {m_k}, {xyz(lepton)}, {m_e})") \
           .Define("kpi1_M",  f"MassCalcs::invMass2({xyz(k)}, {m_k}, {xyz(pi1)}, {m_pi})") \
           .Define("kpi2_M",  f"MassCalcs::invMass2({xyz(k)}, {m_k}, {xyz(pi2)}, {m_pi})") \
           .Define("kpi3_M",  f"MassCalcs::invMass2({xyz(k)}, {m_k}, {xyz(pi3)}, {m_pi})") \
           .Define("k3pi_M",  f"MassCalcs::invMass4({xyz(k)}, {m_k}, {xyz(pi1)}, {m_pi}, {xyz(pi2)}, {m_pi}, {xyz(pi3)}, {m_pi})") \
           .Define("pk3pi_M", f"MassCalcs::invMass5({xyz(p)}, {m_p}, {xyz(k)}, {m_k}, {xyz(pi1)}, {m_pi}, {xyz(pi2)}, {m_pi}, {xyz(pi3)}, {m_pi})") \
           .Define("kpi1pi2_M", f"MassCalcs::invMass3({xyz(k)}, {m_k}, {xyz(pi1)}, {m_pi}, {xyz(pi2)}, {m_pi})") \
           .Define("kpi1pi3_M", f"MassCalcs::invMass3({xyz(k)}, {m_k}, {xyz(pi1)}, {m_pi}, {xyz(pi3)}, {m_pi})") \
           .Define("kpi2pi3_M", f"MassCalcs::invMass3({xyz(k)}, {m_k}, {xyz(pi2)}, {m_pi}, {xyz(pi3)}, {m_pi})") \
           .Define("pk3pie_M",  f"MassCalcs::invMass6({xyz(p)}, {m_p}, {xyz(k)}, {m_k}, {xyz(pi1)}, {m_pi}, {xyz(pi2)}, {m_pi}, {xyz(pi3)}, {m_pi}, {xyz(lepton)}, {m_e})") \
           .Define("pktaue_M",  f"MassCalcs::invMass4({xyz(p)}, {m_p}, {xyz(k)}, {m_k}, {xyz(tau)}, {m_tau}, {xyz(lepton)}, {m_e})") \
           .Define("pk_M",      f"MassCalcs::invMass2({xyz(p)}, {m_p}, {xyz(k)}, {m_k})") \
           .Define("pi1pi2_M",  f"MassCalcs::invMass2({xyz(pi1)}, {m_pi}, {xyz(pi2)}, {m_pi})") \
           .Define("pi1pi3_M",  f"MassCalcs::invMass2({xyz(pi1)}, {m_pi}, {xyz(pi3)}, {m_pi})") \
           .Define("pi2pi3_M",  f"MassCalcs::invMass2({xyz(pi2)}, {m_pi}, {xyz(pi3)}, {m_pi})") \
           .Define("sq_pi1pi2_M", f"MassCalcs::invMassSqScaled2({xyz(pi1)}, {m_pi}, {xyz(pi2)}, {m_pi})") \
           .Define("sq_pi1pi3_M", f"MassCalcs::invMassSqScaled2({xyz(pi1)}, {m_pi}, {xyz(pi3)}, {m_pi})") \
           .Define("sq_pi2pi3_M", f"MassCalcs::invMassSqScaled2({xyz(pi2)}, {m_pi}, {xyz(pi3)}, {m_pi})") \
           .Define("kepi2_M", f"MassCalcs::invMass3({xyz(k)}, {m_k}, {xyz(lepton)}, {m_e}, {xyz(pi2)}, {m_pi})") \
           .Define("ke_M",    f"MassCalcs::invMass2({xyz(k)}, {m_k}, {xyz(lepton)}, {m_e})")

    # Define Mis-ID masses
    df = df.Define("pk_kk_M",          f"MassCalcs::invMass2({xyz(p)}, {m_k}, {xyz(k)}, {m_k})") \
           .Define("pke_pkpi_M",       f"MassCalcs::invMass3({xyz(p)}, {m_p}, {xyz(k)}, {m_k}, {xyz(lepton)}, {m_pi})") \
           .Define("p3pi_k3pi_M",      f"MassCalcs::invMass4({xyz(p)}, {m_k}, {xyz(pi1)}, {m_pi}, {xyz(pi2)}, {m_pi}, {xyz(pi3)}, {m_pi})") \
           .Define("pk3pie_pk3pipi_M", f"MassCalcs::invMass6({xyz(p)}, {m_p}, {xyz(k)}, {m_k}, {xyz(pi1)}, {m_pi}, {xyz(pi2)}, {m_pi}, {xyz(pi3)}, {m_pi}, {xyz(lepton)}, {m_pi})") \
           .Define("ke_kpi_M",         f"MassCalcs::invMass2({xyz(k)}, {m_k}, {xyz(lepton)}, {m_pi})") \
           .Define("kepi2_kpipi_M",    f"MassCalcs::invMass3({xyz(k)}, {m_k}, {xyz(lepton)}, {m_pi}, {xyz(pi2)}, {m_pi})")

    vetoes_expr = mass_vetoes(lepton=lepton)

    # Filter by vetoes and generate passVetoes flag branch
    df = df.Define("passVetoes", vetoes_expr)

    return df


if __name__ == '__main__':

    import sys

    parser = argparse.ArgumentParser(
        description=
        'Update the ntuple with masses calculated with alternative mass hypotheses and apply vetoes.'
    )
    parser.add_argument(
        '--inFile', type=str, required=True, help='input file name')
    parser.add_argument(
        '--tree',
        type=str,
        default='Hlt2RD_LbToPKTauMu_TauTo3Pi_OS_DTF3/DecayTree',
        help='output file name')
    parser.add_argument(
        '--outFile', type=str, required=True, help='output file name')
    parser.add_argument(
        '--saveFiltered',
        action='store_true',
        default=False,
        help=
        'save intermediate ntuple (unsupported/unnecessary in RDataFrame scheme)'
    )

    args = parser.parse_args()

    df = RDataFrame(args.tree, args.inFile)

    # Determine whether we are processing "PplusEplus" mode
    is_PplusEplus = "PplusEplus" in args.inFile

    df = run_veto_pipeline(df)

    df = df.Filter("Lb_BKGCAT == 50")

    print("\nEntries after tm_cut = ", df.Count().GetValue())

    df = df.Filter("passVetoes")

    print("\nEntries after veto = ", df.Count().GetValue())

    if not args.saveFiltered:
        sys.exit(0)

    # Filter the dataset and snapshot (write to disk)
    df.Snapshot("DecayTree", args.outFile)

    print(f"\nFile processed and saved successfully to: {outFile}")
