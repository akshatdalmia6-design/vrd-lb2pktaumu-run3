#ifndef FELDMANCOUSINS_H
#define FELDMANCOUSINS_H

// std
#include <map>
#include <string>
#include <vector>

// RooFit
#include "RooAbsPdf.h"
#include "RooDataSet.h"
#include "RooRealVar.h"

// ROOT
#include "RooRandom.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"

class FeldmanCousins {
public:
  // Helper structure to store single-toy evaluation results
  struct ToyResult {
    double nsig_hat;  // fitted/measured Nsig in toy
    double delta_nll; // -2 * ln Q (profile likelihood test statistic)
  };

  // 1D Constructor
  FeldmanCousins( RooDataSet* data, RooAbsPdf* pdf, RooRealVar* Nsig, const std::string filename = "" );

  ~FeldmanCousins();

  void setBins( const int numNsig, const double minNsig, const double maxNsig );

  TH1* makeNsigHistogram( const int n = 1000 );

  // NEW: Constructs the 2D Confidence Band (True Nsig vs Measured Nsig)
  TH2D* makeConfidenceBandPlot( const int nToys = 1000, const double cl = 0.95, const int beltMode = 0 );

  TGraph* getBeltGraph() { return m_belt; }

  double getMinNLL() { return m_min_nll; }

  void setMinNLL( double val ) { m_min_nll = val; }

  void setMinNLLNsig( double val ) { m_min_nll_Nsig = val; }

  double getBestNsig() const { return m_best_Nsig; }

  std::pair<double, double> pointInNsig( const double Nsig, const int n = 1000 );

  double scanPointInNsig( RooAbsData* data, const double Nsig, const bool useCache = false );

  void setRandomSeed( int seed ) { RooRandom::randomGenerator()->SetSeed( seed ); }

private:
  int numParameters( RooAbsData* data ) const;

  int numFloatingParameters( RooAbsData* data ) const;

  std::pair<double, double> pvalue( const double vdata, std::vector<double>& vtoys ) const;

  double loglikelihood( RooAbsData* data, const double Nsig );

  void histogram( const std::string name, std::vector<double>& vtoys );

  void resetVariables( const bool useCache = true );

  void cacheVariables();

public:
  // debug flag
  static bool debug;

private:
  std::map<std::string, double> m_params;

  TFile* m_file;

  RooAbsData* m_data;
  RooAbsPdf*  m_pdf;

  RooRealVar* m_Nsig;

  // minimum in data
  double m_min_nll;
  double m_min_nll_Nsig = 0.;
  double m_best_Nsig    = 0.;

  int m_num_AFB  = 10;
  int m_num_Nsig = 10;

  double m_min_Nsig = 0.;
  double m_max_Nsig = 20.;

  TGraph* m_belt = nullptr;
};

#endif
