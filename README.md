# vrd-lb2pktaumu-run3
================================================================================
 Lb -> pK tau mu  (Run 3, 2024)   --   analysis code guide
================================================================================
 
Search for the LFV decay Lb0 -> p K- tau mu, with tau -> 3pi(pi0) nu, using
2024 data. Normalisation channel: Lb0 -> J/psi(mumu) K- p.
 
This file explains what each script does and, more importantly, in WHICH ORDER
to run them. Several scripts read caches written by others, so the order is not
optional. Read section 2 before running anything.
 
 
--------------------------------------------------------------------------------
0. ENVIRONMENT
--------------------------------------------------------------------------------
Needs:
  - ROOT >= 6.30 with PyROOT + RooFit  (lb-conda default, or an LCG view)
  - python: numpy, pandas, scipy, scikit-learn, xgboost, matplotlib, tabulate
  - ../utils/lhcbPlotStyle.py  and  lhcbStyle.py  on the path
    (every plotting script does sys.path.append("../utils"))
 
Build the C++ Feldman-Cousins class before any limit setting:
 
    mkdir -p build && cd build
    cmake ..            # picks up $ROOTSYS
    make 
    cd ..
 
This produces build/libFitModel.so, loaded by LFV_model.compileBaseCode().
Sources: FeldmanCousins.cpp / FeldmanCousins.h (+ LinkDef.h for the dictionary).
CMakeLists.txt expects them under src/ and include/.
 
 
--------------------------------------------------------------------------------
1. WHAT EACH FILE IS
--------------------------------------------------------------------------------
Bookkeeping / helpers
  run3_samples_v1r5791.py  Dict of all EOS tuple paths (MC + data, per week).
                           Edit here when new productions arrive.
  simplePlots.py           get_all_input_files() (path globbing) + quick
                           N-1 / variable plots. Also a scratchpad.
  utilities.py             createFolder(), silenceRooFit().
  compare_23.py            Run 2 vs Run 3 variable comparisons. Source of the
                           plotInCanvas() plotting style used elsewhere.
  cuts_check.py            Cut-by-cut yields, Run 2 vs Run 3 selections.
  os_ss_sim.py             OS data / SS data / signal MC overlays, separation
  os_ss_sim2.py            power ranking of BDT input candidates.
  n_skim.py                Light skimmer used for the normalisation channel.
 
Selection chain
  skim_files_comp.py       Skims raw tuples -> local ROOT files. Truth matching,
                           DTF convergence, Lb quality cuts, MC weights, and the
                           generated-event counts. RUN IN TWO PASSES (see below).
  bdt_run3.py              BDT training driver: features, sample building,
                           k-fold training, scored tuple. Main config lives at
                           the top (skim_dir, BKG_MODE, x_min/x_max, NSIGMA).
  bdt_util.py              XGBoost helpers: train_kfold, ROC, importance,
                           correlation matrices, save/load models.
  sb_fit.py                Fits the mass sidebands with an exponential, defines
                           the +-3 RMS90 signal window, writes the cache npz
                           that everything downstream reads.
  punzi_scan.py            Scans the BDT threshold, Punzi FOM eps/(a/2+sqrt(B)),
                           signal efficiency (BR-weighted over the two tau
                           modes) and the cut-flow efficiency table.
  mass_vetoes.py           D0/Jpsi/etc. mass vetoes, computed with JIT'd C++.
  make_fit_tuples.py       Builds the final flat tuples for the fit: signal MC
                           (both tau modes, weighted, mode-tagged) and OS data
                           (sidebands only unless unblinding).
 
Fitting / limit
  LFV_model.py             The RooFit model (signal + combinatorial mass PDFs),
                           mass ranges, and the x-axis label. mode="norm" or
                           "signal".
  SampleAdmin.py           Loads the fit tuples into RooDataSets (year,
                           polarity, PID selection, BDT cut).
  DefaultParameters.py     Default/fixed/set parameter dictionaries, JSON
                           read/write, covariance and correlation printing.
  Plotting.py              draw1Dprojection() etc: fit projections with pull
                           panels, legends, LHCb label.
  test_data.py             MAIN fit driver: fits MC, sidebands, and data.
  test_toy.py              Toy generation + fits; entry point for the FC scan.
  test_fc.py               Feldman-Cousins scan, CL curve, confidence belt,
                           upper limits.
  test_plot_toys.py        Pull/residual summaries from a toy JSON (ROOT).
  plot_toy_distr.py        Same but matplotlib.  usage: python plot_toy_distr.py toys.json
 
 
--------------------------------------------------------------------------------
2. RUN ORDER  (the important bit)
--------------------------------------------------------------------------------
There is a circular-looking dependency: the skim's second pass needs the signal
window, but the window comes from a fit that needs the skims. It is resolved by
splitting skim_files_comp.py into two passes, marked in its __main__ block as
"run first 1 out of 2" and "run second 2 out of 2". Comment out the block you
are not running.
 
  STEP 1   skim_files_comp.py      (pass 1: the four skim() calls)
           -> skims/taumu_sim_OS.root        signal MC, tau->3pi, truth-matched
              skims/taumu_pi0_sim_OS.root    signal MC, tau->3pi pi0
              skims/taumu_dat_SS.root        SS data (background proxy A)
              skims/taumu_dat_OS.root        OS data (background proxy B)
           Prints n_gen and the MC weights. Slow, needs EOS access.
 
  STEP 2   python bdt_run3.py
           -> models  <MODEL_STEM>_fold*.json
              scored tuple <SCORED_TUPLE>
              plots/: roc.pdf, importance.pdf, corr_*.pdf, bdt_corr.pdf
           BKG_MODE at the top selects the background proxy (ss / os / combined).
 
  STEP 3   python sb_fit.py
           -> cache/sideband_fit_<BKG_MODE>.npz   {R, lam, window_lo, window_hi,
                                                   bdt_cut, score_sb}
              plots/sideband_fit.pdf
           The window is peak +- NSIGMA*RMS90, from get_windows().
 
  STEP 4   skim_files_comp.py      (pass 2: the counting block)
           Reads the window from the cache via punzi_scan.load_cache().
           -> skims/taumu_gen_counts.npz, skims/taumu_pi_gen_counts.npz
              {n_gen, n_reco, n_win}
 
  STEP 5   python punzi_scan.py
           -> optimal BDT threshold, signal efficiency +- stat +- syst,
              cut-flow table, plots/punzi_fom.pdf, plots/punzi_eff.pdf
           Also writes n_bdt back into the two gen_counts npz files.
           TAKE THE OPTIMAL THRESHOLD FROM HERE -> it becomes --bdtcut later.
 
  STEP 6   python make_fit_tuples.py
           -> OUT_MC   combined signal MC tuple (bdt, Lb_M, mc_weight, mode)
              OUT_OS   OS data tuple
           OS_REGION = "sideband" keeps the box blind. Set to "full" only when
           you are ready to unblind; it prints a warning when you do.
 
  STEP 7   Normalisation fit:
             python test_data.py --mode norm --fitMC --fitMass --silent
           -> plots/test_fit_norm_2024_Both_*.pdf
              parameters/test_fit_norm_2024_Both.json
           Fit MC first (shape), then data. The JSON is the parameter file
           consumed by the signal fit and by the toys.
 
  STEP 8   Signal fit:
             python test_data.py --mode signal --bdtcut <from step 5> \
                                 --fitMC --fitSideband --fitMass \
                                 --param parameters/test_fit_norm_2024_Both.json
           -> plots/test_fit_signal_2024_Both_{sig_sim_mass,sideband_mass,
                                               data_mass,data_total}WithPulls.pdf
              parameters/test_fit_signal_2024_Both.json
 
  STEP 9   Toys / validation:
             python test_toy.py --param parameters/test_fit_signal_2024_Both.json \
                                --nsig 0 --ncom <expected> --ntoys 1000 \
                                --output toys_nsig0 --silent
             python test_plot_toys.py --input results/toys_nsig0.json
             python plot_toy_distr.py results/toys_nsig0.json
           Check pull means/widths are ~0/1 before trusting the limit.
 
  STEP 10  Feldman-Cousins limit:
             python test_toy.py --param parameters/test_fit_signal_2024_Both.json \
                                --nsig <observed> --ncom <expected> \
                                --doFeldmanCousinsScan --fcbins 40 0 20 \
                                --fctoys 1000 --plot
           -> plots/*_nsig_1DCL.pdf          CL vs N_sig, with the 90/95% lines
              plots/*_nsig_band_95CL.pdf     confidence belt
              results/FC_scan_*.json/root
           The printed "N_sig < x" at 95% CL is the number to quote.
 
 
--------------------------------------------------------------------------------
3. CONFIGURATION TO CHECK BEFORE A RERUN
--------------------------------------------------------------------------------
  skim_files_comp.py   skim_dir, friend_dir, sample_configs, the cut strings
  bdt_run3.py          skim_dir, BKG_MODE, x_min/x_max, NSIGMA, INPUT_FEATURES,
                       N_SPLITS, SEED, MODEL_STEM, OUT_DIR, SCORED_TUPLE
  bdt_util.py          DEFAULT_XGB_PARAMS
  sb_fit.py            BDT_CUT, x_min2/x_max2, CACHE_FILE
  punzi_scan.py        A_PUNZI (=5), N_GRID, T_START, N_MIN, BR/BRE
  make_fit_tuples.py   OS_REGION (blinding!), BR_3PI, BR_3PIPI0, file paths
  test_data.py         CACHE_FILE path, MASS_LO/MASS_HI per mode, binning
  LFV_model.py         mass ranges, PDF shapes, mass label
 
Paths are hard-coded to /afs/.../adalmia/... in several places. Change skim_dir
in BOTH skim_files_comp.py and bdt_run3.py, and CACHE_FILE in test_data.py, if
you move the working area.
 
 
--------------------------------------------------------------------------------
4. IMPORTANT
--------------------------------------------------------------------------------
  - Blinding. make_fit_tuples.py with OS_REGION="sideband" is the safe default.
    Do not flip it until the selection, window and systematics are frozen.
 
  - The two-pass skim. If pass 2 crashes with a missing cache, you skipped
    sb_fit.py. If punzi_scan.py crashes on the npz, you skipped pass 2.
 
  - Efficiency definition. punzi_scan.signal_efficiency() computes a separate
    efficiency per tau mode (N_pass / n_gen for that mode) and combines them as
    eps = (BR0*eps0 + BR1*eps1) / (BR0+BR1). BR_TOT = 0.0931 + 0.0462 = 0.1393
    is the sum of the two SIMULATED modes only, not the full tau BR. The same
    0.1393 must appear in the branching-fraction formula downstream.
 
  - The mc_weight column in the fit tuple (BR_m / n_gen_m) is for the yield
    side; it is NOT used by the efficiency calculation.
 
  - Seeds. --seed is set in test_data.py and test_toy.py; the BDT seed is SEED
    in bdt_run3.py. Fix all of them to reproduce a number exactly.
 
  - Reruns of bdt_run3.py overwrite the fold models, which invalidates the
    sb_fit cache and everything after it. Retrain -> redo steps 3 to 10.
 
  - Only the OOF-scored events populate mode 0 inside the window in
    make_fit_tuples.py (the "scored rest" frame is sideband-only), so mode 0 and
    mode 1 are scored slightly differently. Worth remembering when comparing
    eps0 and eps1.
 
 
--------------------------------------------------------------------------------
5. MINIMAL PATH TO THE FINAL NUMBER
--------------------------------------------------------------------------------
  skim (pass 1) -> bdt_run3 -> sb_fit -> skim (pass 2) -> punzi_scan
  -> make_fit_tuples -> test_data (norm) -> test_data (signal)
  -> test_toy (validation) -> test_toy --doFeldmanCousinsScan
================================================================================
