import ROOT

path_run2 = "/afs/cern.ch/work/f/fabudine/public/vrd-lb2pktaumu/tuples/mcLb_pKtaumu_3pi_2018_MagUp.root"
path_run3 = "/eos/lhcb/grid/prod/lhcb/anaprod/lhcb/MC/2024/TUPLE.ROOT/00384527/0000/00384527_00000001_1.tuple.root"

chain2 = ROOT.TChain("LbTuple1/DecayTree")
chain2.Add(path_run2)
chain3 = ROOT.TChain("MCDT_tuple1/MCDecayTree")
chain3.Add(path_run3)

tm_cut = "Lb_BKGCAT == 50"
cuts_run2 = [
    "proton_PT > 250", "kaon_PT > 250", "kaon_P > 2000", "pi1_PT > 250 ",
    "pi1_P > 2000", "pi2_PT > 250", "pi2_P > 2000", "pi3_PT > 250 ",
    "pi3_P > 2000", "tau_PT > 1000", "mu_isMuon == 1",
    "proton_IPCHI2_OWNPV > 9", "kaon_IPCHI2_OWNPV > 9", "mu_IPCHI2_OWNPV > 9",
    "pi1_MINIPCHI2 > 16", "pi2_MINIPCHI2 > 16", "pi3_MINIPCHI2 > 16",
    "tau_DIRA_OWNPV > 0.99"
]
cuts_run3 = [
    "(pi1_IPCHI2_OWNPV > 25 || pi2_IPCHI2_OWNPV > 25 || pi3_IPCHI2_OWNPV > 25)",
    "(pi1_PT > 800 || pi2_PT > 800 || pi3_PT > 800)", "proton_PIDp > 2",
    "(proton_PIDp - proton_PIDK) > 0", "kaon_PIDK > 8", "mu_PIDmu > 2",
    "pi1_PIDK < 5", "pi1_PIDe < 5", "pi1_PIDp < 10", "pi1_PIDmu < 10",
    "pi2_PIDK < 5", "pi2_PIDe < 5", "pi2_PIDp < 10", "pi2_PIDmu < 10",
    "pi3_PIDK < 5", "pi3_PIDe < 5", "pi3_PIDp < 10", "pi3_PIDmu < 10",
    "proton_P > 8500", "L1520_M > 1300 && L1520_M < 5620", "L1520_PT > 1000",
    "L1520_ENDVERTEX_CHI2/1.0 < 9", "L1520_FDCHI2_OWNPV > 16", "mu_PT > 500",
    "mu_P > 3000", "tau_M > 900 && tau_M < 1550", "tau_FDCHI2_OWNPV > 100",
    "tau_ENDVERTEX_CHI2/3.0 < 6", "Lb_M > 2000 && Lb_M < 10000", "Lb_P > 3000",
    "Lb_PT > 1000", "Lb_FDCHI2_OWNPV > 120", "Lb_ENDVERTEX_CHI2/5.0 < 80",
    "Lb_IPCHI2_OWNPV < 40", "Lb_DIRA_OWNPV > 0.995"
]

total_run2 = chain2.GetEntries(tm_cut)
print(f"Total entries in Run2: {total_run2}\n")

for cut in cuts_run2:
    pas = chain2.GetEntries(f"{cut} && {tm_cut}")
    print(f"Cut: {cut}")
    print(f"  Pass: {pas}\n")
'''
all_cuts2 = " && ".join(f"({c})" for c in cuts_run2)
pas2 = chain2.GetEntries(all_cuts2)
print(f"After all Run2 cuts:")
print(f"  Pass: {pas2}\n")


total_run3 = chain3.GetEntries()
print(f"Total entries in Run3: {total_run3}\n")
'''

results = []

for cut in cuts_run3:
    pas = chain2.GetEntries(f"{cut} && {tm_cut}")
    print(f"Cut: {cut} and {tm_cut}")
    print(f"  Pass: {pas}\n")
    results.append((cut, pas))

all_cuts3 = " && ".join(f"({c})" for c in cuts_run3) + f" && ({tm_cut})"
pas3 = chain2.GetEntries(all_cuts3)
print(f"After all Run3 cuts:")
print(f"  Pass: {pas3}\n")

cut_width = max(len(c) for c, _ in results)
cut_width = max(cut_width, len("All cuts combined"), len("Efficiency"))

print()
print("─" * (cut_width + 18))
print(f"{f'Cuts with {tm_cut}':<{cut_width}}   {'% cut off':>12}")
print("─" * (cut_width + 18))

rows = [(cut, 100.0 * (total_run2 - pas) / total_run2 if total_run2 else 0.0)
        for cut, pas in results]
rows.sort(key=lambda x: x[1], reverse=True)

for cut, cutoff_pct in rows:
    print(f"{cut:<{cut_width}}   {cutoff_pct:>11.2f} %")
print("─" * (cut_width + 18))
efficiency = pas3 / 502314.0
print(f"{'Efficiency':<{cut_width}}   {efficiency:>13.4f}")
print("─" * (cut_width + 18))
efficiency_drop = 100.0 * (0.007 - efficiency) / 0.007
print(f"{'Efficiency Drop Off':<{cut_width}}   {efficiency_drop:>11.2f} %")
print("─" * (cut_width + 18))
