import ROOT

SIGNAL_SIM_FILE = "scored_mc_combined.root"
SIGNAL_DAT_FILE = "scored_data_os.root"
SIGNAL_TREE = "tree"

NORM_DIR = " "
SIGNAL_DIR = " "


class FitSamples(object):
    def __init__(self,
                 mode="norm",
                 year="all",
                 Polarity="Both",
                 fitVars=None,
                 setData=True,
                 setMC=True,
                 PIDsel="default",
                 bdt_cut=None):

        # self.cache_folder = "/eos/user/f/fabudine/B2mumupi/cache/angular"
        # self.cache_folder_data = "/eos/user/f/fabudine/B2mumupi/cache/angular"
        self.bdt_cut = bdt_cut

        #self.base_sel = "survive_preselection && survive_full_fiducial && survive_trigger_selection"
        #self.bdt_sel = "xgb_output > 0.984"
        #self.clone_veto = "acos(samesign_costheta) > 1e-3"
        #self.D0veto = "abs(dimuon_M_Kpi - 1865) > 25 && abs(dimuon_M_piK - 1865) > 25"
        #self.jpsiveto = "survive_single_misidveto"
        #self.mctruth = "B_BKGCAT < 51"

        #self.ProbNNpi = 0.5
        #self.ProbNNK = 0.2

        self.modes = {
            "norm": {
                "sim_dir": NORM_DIR,
                "dat_dir": NORM_DIR,
                "sigID": "mumu_norm_sim.root",
                "data": "mumu_norm_dat.root",
                "tree": "DecayTree",
                "mcweight": ""
            },
            "signal": {
                "sim_dir": SIGNAL_DIR,
                "dat_dir": SIGNAL_DIR,
                "sigID": SIGNAL_SIM_FILE,
                "data": SIGNAL_DAT_FILE,
                "tree": SIGNAL_TREE,
                "mcweight": "mc_weight",
            }
        }

        self.dat_selection = ""
        self.sim_selection = ""

        #self.dimuon_M = ROOT.RooRealVar("dimuon_M", "", 0)
        #self.xgb_output = ROOT.RooRealVar("xgb_output", "", 0)
        #self.hadron_ProbNNpi = ROOT.RooRealVar("hadron_ProbNNpi", "", 0)
        #self.hadron_ProbNNk = ROOT.RooRealVar("hadron_ProbNNk", "", 0)
        #self.dimuon_M_Kpi = ROOT.RooRealVar("dimuon_M_Kpi", "", 0)
        #self.dimuon_M_piK = ROOT.RooRealVar("dimuon_M_piK", "", 0)
        #self.survive_single_misidveto = ROOT.RooRealVar(
        # "survive_single_misidveto", "", 0)
        #self.survive_D0_veto = ROOT.RooRealVar("survive_D0_veto", "", 0)
        #self.proportion_and_pid_weight = ROOT.RooRealVar(
        #"proportion_and_pid_weight", "", 0)
        #self.proportion_and_pid_kin_weight = ROOT.RooRealVar(
        #"proportion_and_pid_kin_weight", "", 0)
        #self.proportion_and_pid_kin_dal_weight = ROOT.RooRealVar(
        #"proportion_and_pid_kin_dal_weight", "", 0)
        self.B_BKGCAT = ROOT.RooRealVar("Lb_BKGCAT", "", 0)
        # self.hadron_TRUEID = ROOT.RooRealVar("hadron_TRUEID", "", 0)
        # self.year = ROOT.RooRealVar("year", "", 0)
        self.Polarity = ROOT.RooRealVar("Polarity", "", 0)

        # self.hadron_ProbNNpi = ROOT.RooRealVar("hadron_ProbNNpi", "", 0)
        # self.hadron_ProbNNk = ROOT.RooRealVar("hadron_ProbNNk", "", 0)

        self.bdt = ROOT.RooRealVar("bdt", "", -1.0, 2.0)

        self.fitSet = ROOT.RooArgSet(*fitVars)
        #  self.fitSet.add(self.dimuon_M)
        #  self.fitSet.add(self.xgb_output)
        if mode == "signal":
            self.fitSet.add(self.bdt)

        fitVarNames = [fitVar.GetName() for fitVar in fitVars]

        self.fitSetMC = self.fitSet.Clone()
        #  self.fitSetMC.add(self.proportion_and_pid_weight)
        if mode == "norm":
            self.fitSetMC.add(self.B_BKGCAT)
        # self.fitSet.add(self.year)
        #  self.fitSet.add(self.Polarity)
        self.mcweight_var = None
        mcweight_name = self.modes[mode]["mcweight"]
        if mcweight_name:
            self.mcweight_var = ROOT.RooRealVar(mcweight_name, "", -1.0e9,
                                                1.0e9)
            self.fitSetMC.add(self.mcweight_var)

        self.samples = {"chains": {}, "dataSets": {}}

        yr = "*" if year == "all" else year
        pol = "*" if Polarity == "Both" else Polarity

        if setData:

            self.setDataSample(mode=mode, yr=yr, pol=pol)

        if setMC:

            self.setMCSample(mode=mode, sample="sig", yr=yr, pol=pol)

    def _selection(self, mode, base, pol=None):
        """Build a clean ' && '-joined cut, appending the BDT cut for signal."""
        pieces = []
        if base:
            pieces.append(base)
        if pol == "Up":
            pieces.append("Polarity == 1")
        elif pol == "Down":
            pieces.append("Polarity == -1")
        if mode == "signal" and self.bdt_cut is not None:
            pieces.append(f"bdt > {self.bdt_cut}")
        return " && ".join(pieces)

    def setMCSample(self, mode="norm", sample="sig", yr="*", pol="*"):
        '''
        Sets simulation samples.
        '''

        self.samples["chains"][sample] = ROOT.TChain(self.modes[mode]["tree"])

        file_path = self.modes[mode]["sim_dir"] + "/" + self.modes[mode][
            f"{sample}ID"]

        print(file_path)

        self.samples["chains"][sample].Add(file_path)

        mcweight = self.modes[mode]["mcweight"]
        selection = self._selection(mode, self.sim_selection)

        print(f"\nSelection for {sample} sample",
              self.modes[mode][f"{sample}ID"])
        print("file_path = ", file_path)
        print(f"cut = {selection}", ", weight: ", mcweight)

        kwargs = dict(Import=self.samples["chains"][sample], Cut=selection)
        if mcweight:
            kwargs["WeightVar"] = mcweight
        self.samples["dataSets"][sample] = ROOT.RooDataSet(
            f"{sample}_" + self.modes[mode][f"{sample}ID"], "", self.fitSetMC,
            **kwargs)

        self.samples["dataSets"][sample].Print()

    def setDataSample(self, mode="norm", yr="*", pol="*"):
        '''
        Sets data samples.
        '''

        self.samples["chains"]["dat"] = ROOT.TChain(self.modes[mode]["tree"])

        file_path = self.modes[mode]["dat_dir"] + "/" + self.modes[mode]["data"]

        self.samples["chains"]["dat"].Add(file_path)

        selection = self._selection(mode, self.dat_selection, pol=pol)

        print(f"\nSelection for data sample")
        print("file_path = ", file_path)
        print(f"cut = {selection}")

        self.samples["dataSets"]["dat"] = ROOT.RooDataSet(
            "data_sample",
            "",
            self.fitSet,
            Import=self.samples["chains"]["dat"],
            Cut=selection)

        self.samples["dataSets"]["dat"].Print()

    def getMCSample(self, sample="sig"):
        '''
        Returns the simulation sample.
        '''
        return self.samples["dataSets"][sample]

    def getDataSample(self):
        '''
        Returns the data sample.
        '''
        return self.samples["dataSets"]["dat"]
