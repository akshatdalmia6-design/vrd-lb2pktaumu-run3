from tabulate import tabulate

import ROOT

isBaseCodeCompiled = False


def compileBaseCode():
    '''
    Jits and loads C++ classes.
    Use jitting to load extra functionality.
    '''

    global isBaseCodeCompiled

    if isBaseCodeCompiled:
        return

    rc = ROOT.gSystem.Load('./build/libFitModel.so')
    if rc < 0:
        raise RuntimeError(
            f"failed to load libFitModel.so (gSystem.Load -> {rc})")
    isBaseCodeCompiled = True
    # ROOT.gInterpreter.ProcessLine('#include "./include/FeldmanCousins.h"')

    # ROOT.gSystem.Load('libRooFit')


class LFVModel(object):
    def __init__(self,
                 reparam=False,
                 com_order=2,
                 use_hypatia=False,
                 mode="norm",
                 com_pos=False,
                 mass_range=None):
        self.var = {}
        self.pdf = {}
        self.fnc = {}

        if mode == "signal":
            _mass_label = '#it{m}(#it{K}^{#minus}#it{p}#it{#tau}^{#pm}#it{#mu}^{#mp})'
        else:
            _mass_label = '#it{m}(#it{K}^{#minus}#it{p}#it{#mu}^{#plus}#it{#mu}^{#minus})'

        if mass_range is None:
            mass_range = (4800., 7000.) if mode == "signal" else (5400., 6000.)

        self.mass = ROOT.RooRealVar("Lb_M", _mass_label, mass_range[0],
                                    mass_range[1], "MeV/#it{c}^{2}")

        # PDF scaling and offset
        self.var['mass:scaling'] = ROOT.RooRealVar(
            'mass:scaling', 'mass:scaling', 1.0, 0.5, 2.0)
        self.var['mass:offset'] = ROOT.RooRealVar('mass:offset', 'mass:offset',
                                                  0.0, -10., 10.)

        add_misid = False if "k" in mode else True

        self._width_seed = 60. if mode == "signal" else 20.
        self._width_max = 250. if mode == "signal" else 30.

        self._mode = mode

        self.createPdfShapes(reparam, com_order, use_hypatia, com_pos)

    def createPdfShapes(self,
                        reparam=False,
                        com_order=2,
                        use_hypatia=False,
                        com_pos=False):
        '''
        Create all of the mass/angular PDFs
        '''
        # Create mass shapes
        self.pdf['sig:mass'] = self.createSigMassPdf(use_hypatia)
        self.pdf['com:mass'] = self.createComMassPdf()

        # Create combined PDF
        self.var['sig:yield'] = ROOT.RooRealVar('sig:yield', 'nsig', 100., 0.,
                                                100000000)
        self.var['com:yield'] = ROOT.RooRealVar('com:yield', 'ncom', 100., 0.,
                                                100000000)

        _shape_list_mass = ROOT.RooArgList(self.pdf['sig:mass'],
                                           self.pdf['com:mass'])  # ,
        # self.pdf['mid:mass'] )
        _shape_list = _shape_list_mass  #,
        # self.pdf['mid:total'] )
        _yield_list = ROOT.RooArgList(self.var['sig:yield'],
                                      self.var['com:yield'])  #,
        # self.var['mid:yield'] )
        '''

        if add_misid:
            self.pdf['mid:mass'] = self.createMidMassPdf()
            #self.pdf['mid:angle'] = self.createAnglePolynomial( 'mid:angle', 4 )
            self.pdf['mid:total'] = self.createProduct('mid')
            self.var['mid:yield'] = ROOT.RooRealVar('mid:yield', 'nmid', 100.,
                                                    0., 100000000)
            _shape_list_mass.add(self.pdf['mid:mass'])
            _shape_list.add(self.pdf['mid:total'])
            _yield_list.add(self.var['mid:yield'])

        '''

        _shape_list_mass.Print()

        # This is for fits using only the B mass
        self.pdf['combined_mass'] = ROOT.RooAddPdf(
            'combined_mass', 'combined_mass', _shape_list_mass, _yield_list)

        # This is for mass fits
        self.pdf['combined'] = ROOT.RooAddPdf('combined', 'combined',
                                              _shape_list, _yield_list)
        return

    def createAnglePolynomial(self, name, n):
        '''
        Create a polynomial model with orders up-to n
        '''
        _arg_list = ROOT.RooArgList()
        for i in range(n):
            _arg_name = name + ':p%i' % (i + 1)
            self.var[_arg_name] = ROOT.RooRealVar(_arg_name, _arg_name, 0,
                                                  -1.0, 1.0)
            _arg_list.add(self.var[_arg_name])
        return ROOT.RooChebychev(name, name, self.costheta, _arg_list)

    def createPosAnglePolynomial(self, name):
        """
        Creates the PDF: 1/(2*(1 + b/3)) * (1 + a*(1+b)*cosTheta + b*cosTheta^2)
        """
        # 1. Define Parameters a and b
        # a in [-1, 1]
        self.var[f"{name}:p1"] = ROOT.RooRealVar(f"{name}:p1", "a", 0.0, -1.0,
                                                 1.0)
        # b in [-1, 10] (Setting a reasonable upper bound for 'infinity')
        self.var[f"{name}:p2"] = ROOT.RooRealVar(f"{name}:p2", "b", 0.5, -1.0,
                                                 10.0)

        # 2. Define the Formula String (C++ Syntax)
        # Note: RooGenericPdf will automatically handle normalisation if
        # the formula is simple, but including the explicit normalisation is safer.
        formula = "(1.0 / (2.0 * (1.0 + @1/3.0))) * (1.0 + @0 * (1.0 + @1) * @2 + @1 * @2 * @2)"

        # 3. Create the ArgList for the formula placeholders (@0, @1, @2)
        # Mapping: @0 -> a, @1 -> b, @2 -> costheta
        _arg_list = ROOT.RooArgList(self.var[f"{name}:p1"],
                                    self.var[f"{name}:p2"], self.costheta)

        # 4. Retutest_fit_signal_2024_Both.jsonrn the Generic PDF
        return ROOT.RooGenericPdf(name, name, formula, _arg_list)

    def createCrystalBall(self, name, m0param=[5619.6, 5400., 6000.]):
        '''
        Create a Crystal Ball shape
        '''
        if not (name + ':m0') in self.var:
            self.var[name + ':m0'] = ROOT.RooRealVar(name + ':m0',
                                                     name + ':m0', *m0param)
            self.fnc[name + ':m0:offset'] = self.createOffsetVariable(
                self.var[name + ':m0'], self.var['mass:offset'])

        self.var[name + ':sC'] = ROOT.RooRealVar(
            name + ':sC', name + ':sC', self._width_seed, 0.1, self._width_max)
        self.fnc[name + ':sC:scaling'] = self.createScaledVariable(
            self.var[name + ':sC'], self.var['mass:scaling'])

        self.var[name + ':aL'] = ROOT.RooRealVar(name + ':aL', name + ':aL',
                                                 1., 0., 10.)
        self.var[name + ':aR'] = ROOT.RooRealVar(name + ':aR', name + ':aR',
                                                 1., 0., 10.)
        self.var[name + ':nL'] = ROOT.RooRealVar(name + ':nL', name + ':nL',
                                                 5., 0., 10.)
        self.var[name + ':nR'] = ROOT.RooRealVar(name + ':nR', name + ':nR',
                                                 5., 0., 10.)

        return ROOT.RooCrystalBall(
            name + ':C', name + ':C', self.mass, self.fnc[name + ':m0:offset'],
            self.fnc[name + ':sC:scaling'], self.var[name + ':aL'],
            self.var[name + ':nL'], self.var[name + ':aR'],
            self.var[name + ':nR'])

    def createGaussian(self, name, m0param=[5619.6, 5400., 6000.]):
        '''
        Create a Gaussian shape
        '''
        if not (name + ':m0') in self.var:
            self.var[name + ':m0'] = ROOT.RooRealVar(name + ':m0',
                                                     name + ':m0', *m0param)
            self.fnc[name + ':m0:offset'] = self.createOffsetVariable(
                self.var[name + ':m0'], self.var['mass:offset'])
        self.var[name + ':sG'] = ROOT.RooRealVar(
            name + ':sG', name + ':sG', self._width_seed, 0.1, self._width_max)
        self.fnc[name + ':sG:scaling'] = self.createScaledVariable(
            self.var[name + ':sG'], self.var['mass:scaling'])

        return ROOT.RooGaussian(name + ':G', name + ':G', self.mass,
                                self.fnc[name + ':m0:offset'],
                                self.fnc[name + ':sG:scaling'])

    def createHypatia(self, name):
        '''
        Create a hypatia shape
        '''
        if not (name + ':m0') in self.var:
            self.var[name + ':m0'] = ROOT.RooRealVar(
                name + ':m0', name + ':m0', 5619.6, 5400., 6000.)
            self.fnc[name + ':m0:offset'] = self.createOffsetVariable(
                self.var[name + ':m0'], self.var['mass:offset'])
        self.var[name + ':sH'] = ROOT.RooRealVar(name + ':sH', name + ':sH',
                                                 20., 0., 50.)
        self.fnc[name + ':sH:scaling'] = self.createScaledVariable(
            self.var[name + ':sH'], self.var['mass:scaling'])
        self.var[name + ':lambda'] = ROOT.RooRealVar(
            name + ':lambda', name + ':lambda', -1.0, -10.0, 0.0)
        self.var[name + ':zeta'] = ROOT.RooRealVar(name + ':zeta',
                                                   name + ':zeta', 0.005)
        self.var[name + ':beta'] = ROOT.RooRealVar(name + ':beta',
                                                   name + ':beta', 0.0)
        self.var[name + ':aL'] = ROOT.RooRealVar(name + ':aL', name + ':aL',
                                                 1., 0., 10.)
        self.var[name + ':aR'] = ROOT.RooRealVar(name + ':aR', name + ':aR',
                                                 1., 0., 10.)
        self.var[name + ':nL'] = ROOT.RooRealVar(name + ':nL', name + ':nL',
                                                 5., 0., 10.)
        self.var[name + ':nR'] = ROOT.RooRealVar(name + ':nR', name + ':nR',
                                                 5., 0., 10.)

        return ROOT.RooHypatia2(
            name + ':H', name + ':H', self.mass, self.var[name + ':lambda'],
            self.var[name + ':zeta'], self.var[name + ':beta'],
            self.fnc[name + ':sH:scaling'], self.fnc[name + ':m0:offset'],
            self.var[name + ':aL'], self.var[name + ':nL'],
            self.var[name + ':aR'], self.var[name + ':nR'])

    def createCrystalBallAndGauss(self, name):
        '''
        Create the sum of a Crystal Ball and a Gaussian shape
        '''
        self.var[name + ':m0'] = ROOT.RooRealVar(name + ':m0', name + ':m0',
                                                 5619.6, 5400., 6000.)
        self.fnc[name + ':m0:offset'] = self.createOffsetVariable(
            self.var[name + ':m0'], self.var['mass:offset'])

        self.var[name + ':fG'] = ROOT.RooRealVar(name + ':fG', name + ':fG',
                                                 0.5, 0.0, 1.0)

        self.pdf[name + ':G'] = self.createGaussian(name)
        self.pdf[name + ':C'] = self.createCrystalBall(name)

        return ROOT.RooAddPdf(name, name, self.pdf[name + ':G'],
                              self.pdf[name + ':C'], self.var[name + ':fG'])

    def createOffsetVariable(self, variable, offset):
        '''
        Offset a mass variable by a measured shift
        '''
        formula = '{0} - {1}'.format(variable.GetName(), offset.GetName())
        _arg_list = ROOT.RooArgList(variable, offset)
        return ROOT.RooFormulaVar(variable.GetName() + ':offset', formula,
                                  _arg_list)

    def createScaledVariable(self, variable, scaling):
        '''
        Scale a width by a known scale factor
        '''
        formula = '{0} * {1}'.format(variable.GetName(), scaling.GetName())
        _arg_list = ROOT.RooArgList(variable, scaling)
        return ROOT.RooFormulaVar(variable.GetName() + ':scaling', formula,
                                  _arg_list)

    def createSigMassPdf(self, use_hypatia=False):
        '''
        Create signal mass PDF: Distribution described by the sum of Gaussian 
        function and a double sided Crystal Ball shape
        '''
        if use_hypatia:
            return self.createHypatia('sig:mass')

        if self._mode == "signal":
            return self.createCrystalBall('sig:mass')

        return self.createCrystalBallAndGauss('sig:mass')

    def createSigAngularPdf(self, angle_var="restframe_lepton_costheta"):
        '''
        Create signal angular PDF
        '''
        self.var['FH'] = ROOT.RooRealVar('FH', 'FH', 0.1, 0.0, 1.0)

        if angle_var == "restframe_lepton_costheta":
            self.var['AFB'] = ROOT.RooRealVar('AFB', 'AFB', 0.0, -0.5, 0.5)
            return ROOT.RooPLeptonAngleWithVeto('sig:angle', 'sig:angle',
                                                self.costheta, self.var['FH'],
                                                self.var['AFB'])
        elif angle_var == "abs_restframe_lepton_costheta":

            return ROOT.RooPLeptonAbsAngleWithVeto(
                'sig:angle', 'sig:angle', self.costheta, self.var['FH'])

    def createSigReparamPdf(self):
        '''
        Create signal angular PDF using reparameterised form
        '''
        self.var['AFBp'] = ROOT.RooRealVar('AFBp', 'AFB', 0.0, -1.0, 1.0)
        self.var['FH'] = ROOT.RooRealVar('FH', 'FH', 0.1, 0.0, 1.0)

        _arg_list = ROOT.RooArgList(self.var['FH'], self.var['AFBp'])
        self.fnc['AFB'] = ROOT.RooFormulaVar('reparam:AFB', '0.5*FH*AFBp',
                                             _arg_list)

        return ROOT.RooPLeptonAngleWithVeto('sig:reparam:angle', 'sig:angle',
                                            self.costheta, self.var['FH'],
                                            self.fnc['AFB'])

    def createComMassPdf(self):
        '''
        Create combinatorial background mass PDF: Distribution described by 
        an exponential function 
        '''
        self.var['com:mass:e0'] = ROOT.RooRealVar('com:mass:e0', 'e0',
                                                  -5.0949e-03, -0.02, 0.02)
        return ROOT.RooExponential('com:mass', 'com:mass', self.mass,
                                   self.var['com:mass:e0'])

    def createMidMassPdf(self):
        '''
        Create mis-id mass PDF: Described by the sum of Crystal Ball shape with
        tails on the left- and right-hand side of the distribution and a 
        Gaussian distribution.
        '''
        return self.createCrystalBall('mid:mass')
        #return self.createCrystalBallAndGauss( 'mid:mass' )

    def createMidAngularPdf(self):
        '''
        Create the angular PDF for the misid background: sum of even ordered 
        Chebychev polynomials up-to and including 4th order
        '''
        self.var['mid:angle:p2'] = ROOT.RooRealVar('mid:angle:p2',
                                                   'mid:angle:p2', 0., -2., 2.)
        self.var['mid:angle:p4'] = ROOT.RooRealVar('mid:angle:p4',
                                                   'mid:angle:p4', 0., -2., 2.)

        return ROOT.RooMisidAngleWithVeto(
            'mid:angle', 'mid:angle', self.costheta, self.var['mid:angle:p2'],
            self.var['mid:angle:p4'])

    def fixParam(self, param={}, fixed=True):
        '''
        Fix parameters in the list by name
        '''
        for k, v in param.items():
            if k in self.var:
                self.var[k].setVal(v)
                self.var[k].setConstant(fixed)
        return

    def fixParams(self, keyword="", fixed=True):
        '''
        Fix mass parameters in the list 
        '''
        for k in self.var.keys():
            if keyword in k:
                self.var[k].setConstant(fixed)
        return

    def setParam(self, param={}):
        '''
        Set parameter values from a dictionary 
        '''
        for k, v in param.items():
            if k in self.var:
                self.var[k].setVal(v)
        return

    def getParam(self, param):
        '''
        Get a parameter value
        '''
        if param in self.var:
            return (self.var[param].getVal(), self.var[param].getError())
        elif param in self.fnc:
            return (self.fnc[param].getVal(), 0.)
        return (0, 0)

    def setLegendreEfficiency(self, coeff):
        '''
        Set the signal angular efficiency
        '''
        self.pdf['sig:angle'].setLegEfficiency(coeff)

    def setSignalVetoEfficiency(self, histogram):
        '''
        Set a relative efficiency shape due to a veto
        '''
        self.pdf['sig:angle'].setVeto(histogram)

    def setMisidVetoEfficiency(self, histogram):
        '''
        '''
        if hasattr(self.pdf['mid:angle'], 'setVeto'):
            self.pdf['mid:angle'].setVeto(histogram)

    def setPolynomialEfficiency(self, coeff):
        '''
        Set the signal angular efficiency
        '''
        self.pdf['sig:angle'].setPolEfficiency(coeff)

    def getSignalAngularPdf(self):
        '''
        Get the signal angular PDF
        '''
        return self.pdf['sig:angle']

    def getSignalMassPdf(self):
        '''
        Get the signal mass PDF
        '''
        return self.pdf['sig:mass']

    def getMidAngularPdf(self):
        '''
        Get the mid angular PDF
        '''
        return self.pdf['mid:angle']

    def getCombinedAngularPdf(self):
        '''
        Get the combined angular PDF
        '''
        return self.pdf['combined_angle']

    def getMidMassPdf(self):
        '''
        Get the misid mass PDF
        '''
        return self.pdf['mid:mass']

    def getComAngularPdf(self):
        '''
        Get the com angular PDF
        '''
        return self.pdf['com:angle']

    def getComMassPdf(self):
        '''
        Get the misid mass PDF
        '''
        return self.pdf['com:mass']

    def getCombinedMassPdf(self):
        '''
        Get the total mass PDF
        '''
        return self.pdf['combined_mass']

    def getCombinedPdf(self):
        '''
        Get the total PDF
        '''
        return self.pdf['combined']

    def generateDataset(self, nsig=-1, ncom=-1, seed=42):
        '''
        Generate a new data set from specified yields
        '''

        if nsig > 0: self.var['sig:yield'].setVal(nsig)
        if ncom > 0: self.var['com:yield'].setVal(ncom)

        model = self.getCombinedPdf()

        if nsig == 0:
            print("\n Debug mass shape")
            self.pdf["com:mass"].Print()

            randNumbGenerator = ROOT.TRandom3(seed)
            poisson_ncom = randNumbGenerator.Poisson(ncom)

            print("poisson_ncom = ", poisson_ncom)

            return self.pdf["com:mass"].generate(
                ROOT.RooArgList(self.mass), poisson_ncom)

        else:
            return model.generate(
                ROOT.RooArgList(self.mass), ROOT.RooFit.Extended())

    def printParam(self):
        '''
        Print the used parameters
        '''
        _data = []
        for k, v in self.var.items():
            _data.append(
                [v.GetName(),
                 v.getVal(),
                 v.getError(),
                 v.isConstant()])
        print('\n')
        print(
            tabulate(
                _data,
                headers=['Parameters', 'Value', 'Error', 'Constant'],
                tablefmt='presto'))
        print('\n')
        return _data

    def printSummary(self, param):
        '''
        Print summary table of parameters in "param" 
        '''
        _data = []
        for k, v in param.items():
            _value = self.getParam(k)
            _data.append([k, _value[0] - v, _value[1]])
        print('\n')
        print(
            tabulate(
                _data,
                headers=['Observable', 'Bias', 'Error'],
                tablefmt='presto'))
        print('\n')
        return _data

    def getNsig(self):
        '''
        Return the FH variable
        '''
        return self.var['sig:yield']

    def getValue(self, param):
        '''
        Get the value of a parameter
        '''
        if param in self.fnc:
            return self.fnc[param].getVal()
        if param in self.var:
            return self.var[param].getVal()
        return 0.

    def getError(self, param):
        '''
        Get the uncertainty on a parameter
        '''
        if param in self.var and not self.var[param].isConstant():
            return self.var[param].getError()
        return 0.

    def getSigPdf(self):
        '''
        Get the signal PDF
        '''
        return self.pdf['sig:total']

    def getMidPdf(self):
        '''
        Get the PDF for the misid background
        '''
        return self.pdf['mid:total']

    def getBkgPdf(self):
        '''
        Get the PDF for the combinatorial background
        '''
        return self.pdf['com:total']

    def getComponent(self, name):
        '''
        Get a single PDF component
        '''
        if name in self.pdf:
            return self.pdf[name]
        return None

    def getFitVars(self):
        '''
        Returns list with fit variables
        '''
        return [self.mass]
