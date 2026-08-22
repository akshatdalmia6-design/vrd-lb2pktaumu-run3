#include "FeldmanCousins.h"

#include "RooArgList.h"
#include "RooFitResult.h"
#include "RooMsgService.h"

#include "RooAbsCollection.h"
#include "TH1D.h"
#include "TRandom3.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

bool FeldmanCousins::debug = false;

// 1D implementation
FeldmanCousins::FeldmanCousins( RooDataSet* data, RooAbsPdf* pdf, RooRealVar* Nsig, const std::string filename )
    : m_data( data ), m_pdf( pdf ), m_Nsig( Nsig ), m_min_nll( 0. ) {

  // Connect the angular fit parameters
  // Get the list of all variables actually used by the PDF
  RooArgSet* vars = m_pdf->getVariables();

  // Find the variable named after "Nsig" in that set
  RooRealVar* activeNsig = (RooRealVar*)vars->find( Nsig->GetName() );
  if ( activeNsig ) {
    m_Nsig = activeNsig; // Update our pointer to the one the PDF actually uses
  }

  delete vars;

  // Do a first fit to get minimum of Log likelihood
  RooFitResult* result = m_pdf->fitTo( *m_data, RooFit::Save( true ), RooFit::PrintLevel( -1 ) );

  if ( !result ) {
    std::cerr << "FeldmanCousins: initial fit to data failed (null RooFitResult). "
              << "Check that at least one parameter is floating." << std::endl;
    throw std::runtime_error( "FeldmanCousins: initial fit failed" );
  }

  m_min_nll = result->minNll();

  if ( debug ) {
    result->Print();
    std::cout << "==> Likelihood minimum (" << m_Nsig->getVal() << ")" << std::endl;
  }

  delete result;

  if ( !filename.empty() ) {
    m_file = new TFile( filename.c_str(), "RECREATE" );
  } else {
    m_file = nullptr;
  }
}

FeldmanCousins::~FeldmanCousins() {

  if ( m_file ) {
    m_file->Close();
    delete m_file;
  }
}

void FeldmanCousins::resetVariables( const bool useCache ) {

  // reset all parameters to their cached values

  RooArgSet* args = m_pdf->getParameters( m_data );

  for ( auto* arg : *args ) {

    auto* v = dynamic_cast<RooRealVar*>( arg );

    if ( !v ) continue;

    if ( !v->isConstant() ) {

      if ( useCache ) {
        // Reset non-constant model parameters to cached values
        auto it = m_params.find( v->GetName() );
        if ( it != m_params.end() ) { v->setVal( it->second ); }
        // v->setVal( m_params[v->GetName()] );
      } else {
        // Reset non-constant model parameters to 0
        // if (v==m_FH){
        //   m_FH->setVal( 0.1 );
        // } else {
        // v->setVal( 0.5 * ( v->getMin() + v->getMax() ) );
        if ( std::string( v->GetName() ).find( "yield" ) != std::string::npos ) {
          v->setVal( 0.5 * double( m_data->numEntries() ) );
        } else {
          v->setVal( 0 );
        }
        // }
      }
      v->removeError();

      if ( debug ) { std::cout << "==> Setting " << v->GetName() << " to " << v->getVal() << std::endl; }
    }
  }

  delete args;

  return;
}

void FeldmanCousins::cacheVariables() {

  // cache parameter values

  RooArgSet* args = m_pdf->getParameters( m_data );

  if ( debug ) { std::cout << "==> Caching variables" << std::endl; }

  for ( auto* arg : *args ) {

    auto* v = dynamic_cast<RooRealVar*>( arg );
    if ( !v ) continue;
    m_params[v->GetName()] = v->getVal();
  }

  delete args;

  return;
}

std::pair<double, double> FeldmanCousins::pvalue( const double vdata, std::vector<double>& vtoys ) const {

  if ( vtoys.empty() ) {
    std::cerr << "FeldmanCousins::pvalue: no surviving toys -- returning p = 1" << std::endl;
    return std::make_pair( 1., 0. );
  }

  std::sort( vtoys.begin(), vtoys.end() );

  auto isearch = std::lower_bound( vtoys.begin(), vtoys.end(), vdata );

  unsigned int n = vtoys.size();
  unsigned int i = std::distance( isearch, vtoys.end() );

  double eff = double( i ) / double( n );
  double err = std::sqrt( eff * ( 1. - eff ) / double( n ) );

  // if ( vdata < 0 ) {
  //   std::cout << "vdata = " << vdata << ", "
  //             << "eff = " << eff << std::endl;
  // };

  return std::make_pair( eff, err );
}

int FeldmanCousins::numParameters( RooAbsData* data ) const {

  RooArgSet* params = m_pdf->getParameters( data );
  int        result = params->getSize();

  delete params;
  return result;
}

int FeldmanCousins::numFloatingParameters( RooAbsData* data ) const {

  RooArgSet* params = m_pdf->getParameters( data );
  int        result = 0;
  for ( auto* arg : *params ) {
    auto* v = dynamic_cast<RooRealVar*>( arg );
    if ( v && !v->isConstant() ) ++result;
  }

  delete params;
  return result;
}

double FeldmanCousins::loglikelihood( RooAbsData* data, const double Nsig ) {

  double value = 0;

  m_Nsig->setVal( Nsig );
  m_Nsig->setConstant( true );

  const bool noNuisance = ( 0 == numFloatingParameters( data ) );

  if ( noNuisance ) {
    // RooNLLVar nll( "nll", "nll", *m_pdf, *data ); // This is old and deprecated
    RooAbsReal* nll = m_pdf->createNLL( *data );
    value           = nll->getVal();
    delete nll;
  } else {
    RooFitResult* result = m_pdf->fitTo( *data, RooFit::Save( true ), RooFit::PrintLevel( -1 ) );

    if ( result ) {
      value = result->minNll();
      delete result;
    } else {
      value = std::numeric_limits<double>::quiet_NaN();
    }
  }

  m_Nsig->setConstant( false );

  if ( debug ) { std::cout << "==> Evaluating likelihood, value = " << value << std::endl; }

  return value;
}

void FeldmanCousins::histogram( const std::string name, std::vector<double>& vtoys ) {

  // If no file was provided, skip the writing process
  if ( !m_file ) return;

  TH1D* hist = new TH1D( name.c_str(), name.c_str(), 101, -0.5, 100.5 );

  if ( debug ) { std::cout << "==> Creating histogram : " << name << std::endl; }

  for ( auto v : vtoys ) { hist->Fill( v ); }

  m_file->cd();
  hist->Write();

  return;
}

std::pair<double, double> FeldmanCousins::pointInNsig( const double Nsig, const int n ) {

  if ( debug ) { std::cout << "==> Point at Nsig  : " << Nsig << std::endl; }

  // Reset variables without using cached values so it starts at the same default point
  resetVariables( false );

  double vfloat = 0;
  double vfixed = 0;
  double vdata  = scanPointInNsig( m_data, Nsig, false );

  cacheVariables();

  std::vector<double> vtoys;

  int n_rejected = 0;

  for ( int i = 0; i < n; i++ ) {

    // Reset variables using cached values obtained from fit to data
    resetVariables( true );

    gRandom->SetSeed( i + 1 );

    const double n_exp = m_pdf->expectedEvents( m_data->get() );
    const int    nData = gRandom->Poisson( n_exp );
    // generate toy data
    RooDataSet* toy = m_pdf->generate( *m_data->get(), nData );

    if ( debug ) {
      std::cout << "==> Created new dataset:" << std::endl;
      toy->Print();
    }

    vfixed = scanPointInNsig( toy, Nsig, false );

    // floating
    resetVariables( false );
    // m_AFB->setVal( 0.0 );
    // m_Nsig->setVal( 0.1 );

    RooFitResult* rfloat = m_pdf->fitTo( *toy, RooFit::Save( true ), RooFit::PrintLevel( -1 ) );
    vfloat               = rfloat ? rfloat->minNll() : std::numeric_limits<double>::quiet_NaN();
    // store the result
    if ( std::isnormal( vfloat ) && vfloat <= vfixed ) {
      vtoys.push_back( vfixed - vfloat );
    } else {
      ++n_rejected;
    }

    // cleanup
    if ( rfloat ) delete rfloat;
    delete toy;
  }

  if ( n_rejected > n / 20 ) {
    std::cout << "==> WARNING: at Nsig = " << Nsig << ", " << n_rejected << " / " << n
              << " toys rejected by the fit-quality cut" << std::endl;
  }

  histogram( std::string( "Nsig" ) + std::to_string( Nsig ), vtoys );

  double diff = vdata - m_min_nll_Nsig;

  if ( diff < 0 ) {

    std::cout << "\n==> Negative dif! ==> " << std::endl;
    std::cout << "==> Point at Nsig  : " << Nsig << std::endl;
    std::cout << "==> Likelihood diff = " << diff << std::endl;
    std::cout << "==> m_min_nll_Nsig = " << m_min_nll_Nsig << ", vdata = " << vdata << std::endl;
    std::cout << "==> (the grid minimum missed the true conditional minimum: "
              << "widen or refine the scan range)" << std::endl;
  }

  return pvalue( vdata - m_min_nll_Nsig, vtoys );
}

double FeldmanCousins::scanPointInNsig( RooAbsData* data, const double Nsig, const bool useCache ) {

  resetVariables( useCache );

  return loglikelihood( data, Nsig );
}

void FeldmanCousins::setBins( const int numNsig, const double minNsig, const double maxNsig ) {
  m_num_Nsig = numNsig;
  m_min_Nsig = minNsig;
  m_max_Nsig = maxNsig;
}

TH1* FeldmanCousins::makeNsigHistogram( const int n ) {
  RooMsgService::instance().setGlobalKillBelow( RooFit::FATAL );

  TH1D* hist = new TH1D( "hist_interval_Nsig", "", m_num_Nsig, m_min_Nsig, m_max_Nsig );
  hist->SetXTitle( "#it{N}_{sig}" );
  hist->SetDirectory( nullptr );

  TH1D* nllhist = (TH1D*)hist->Clone( "hist_nll_Nsig" );
  nllhist->Reset();
  nllhist->SetDirectory( nullptr );

  for ( int i = 0; i < nllhist->GetNbinsX(); i++ ) {
    double vNsig = nllhist->GetXaxis()->GetBinCenter( i + 1 );
    double vdata = scanPointInNsig( m_data, vNsig, false );
    nllhist->SetBinContent( i + 1, vdata );
    // std::cout << "\n ==> Point at vNsig  : " << std::endl;
    // std::cout << "Old min = " << m_min_nll << ", New min = " << vdata << std::endl;
  }

  // Find bin containing optimal value
  int min_bin    = nllhist->GetMinimumBin();
  m_min_nll_Nsig = nllhist->GetBinContent( min_bin );
  m_best_Nsig    = nllhist->GetXaxis()->GetBinCenter( min_bin );

  if ( min_bin == 1 || min_bin == nllhist->GetNbinsX() ) {
    std::cout << "==> WARNING: conditional NLL minimum is at the edge of the scan range (" << m_best_Nsig
              << "). Widen setBins()." << std::endl;
  }

  delete nllhist;

  for ( int i = 0; i < hist->GetNbinsX(); i++ ) {
    double vNsig      = hist->GetXaxis()->GetBinCenter( i + 1 );
    auto [pval, perr] = pointInNsig( vNsig, n );
    hist->SetBinContent( i + 1, pval );
    hist->SetBinError( i + 1, perr );
  }

  RooMsgService::instance().reset();

  return hist;
}

static std::pair<double, double> shortestInterval( const std::vector<double>& v, size_t nKeep ) {

  if ( v.empty() ) return {0., 0.};
  if ( nKeep >= v.size() ) return {v.front(), v.back()};
  if ( nKeep < 1 ) nKeep = 1;

  size_t best_i = 0;
  double best_w = v[nKeep - 1] - v[0];

  for ( size_t i = 1; i + nKeep <= v.size(); ++i ) {
    const double w = v[i + nKeep - 1] - v[i];
    if ( w < best_w ) {
      best_w = w;
      best_i = i;
    }
  }

  return {v[best_i], v[best_i + nKeep - 1]};
}

TH2D* FeldmanCousins::makeConfidenceBandPlot( const int nToys, const double cl, const int beltMode ) {
  RooMsgService::instance().setGlobalKillBelow( RooFit::FATAL );

  // Create 2D plot: X-axis = Measured Nsig (Nsig_hat), Y-axis = True Nsig
  std::string hName    = "h_FC_ConfidenceBand_" + std::to_string( int( cl * 100 ) );
  TH2D*       bandHist = new TH2D( hName.c_str(), ";#hat{#it{N}}_{sig} (Measured);#it{N}_{sig} (True)", m_num_Nsig,
                             m_min_Nsig, m_max_Nsig, m_num_Nsig, m_min_Nsig, m_max_Nsig );
  bandHist->SetDirectory( nullptr );

  double stepNsig = ( m_max_Nsig - m_min_Nsig ) / double( m_num_Nsig );

  std::vector<double> belt_y, belt_lo, belt_hi;

  for ( int i = 0; i < m_num_Nsig; ++i ) {
    double trueNsig = m_min_Nsig + ( i + 0.5 ) * stepNsig;
    if ( debug ) { std::cout << "==> Generating FC Belt slice for True Nsig = " << trueNsig << std::endl; }

    std::vector<ToyResult> toys;
    toys.reserve( nToys );

    resetVariables( false );

    scanPointInNsig( m_data, trueNsig, false );
    cacheVariables();

    std::vector<double> all_hats;
    int                 n_bad_fit = 0;

    // Generate Toys at trueNsig
    for ( int t = 0; t < nToys; ++t ) {
      resetVariables( true );
      gRandom->SetSeed( t + 1 );

      const double n_exp = m_pdf->expectedEvents( m_data->get() );
      const int    nData = gRandom->Poisson( n_exp );
      RooDataSet*  toy   = m_pdf->generate( *m_data->get(), nData );

      // 1. Evaluate fixed Nll (given trueNsig)
      double vfixed = scanPointInNsig( toy, trueNsig, false );

      // 2. Evaluate floating fit to find best-fit Nsig (Nsig_hat) and minimum NLL
      resetVariables( false );
      RooFitResult* rfloat = m_pdf->fitTo( *toy, RooFit::Save( true ), RooFit::PrintLevel( -1 ) );

      double vfloat   = rfloat ? rfloat->minNll() : std::numeric_limits<double>::quiet_NaN();
      double nsig_hat = m_Nsig->getVal();

      all_hats.push_back( nsig_hat );
      if ( !( std::isnormal( vfloat ) && vfloat <= vfixed ) ) ++n_bad_fit;

      if ( std::isnormal( vfloat ) && vfloat <= vfixed ) { toys.push_back( {nsig_hat, vfixed - vfloat} ); }

      if ( rfloat ) delete rfloat;
      delete toy;
    }

    {
      std::vector<double> h = all_hats;
      std::sort( h.begin(), h.end() );
      auto q = [&]( double f ) { return h[std::min( h.size() - 1, size_t( f * double( h.size() ) ) )]; };
      std::cout << "  slice true=" << trueNsig << "  nsig_hat p2.5=" << q( 0.025 ) << " median=" << q( 0.50 )
                << " p97.5=" << q( 0.975 ) << "  |  ncom_last=" << m_params["com:yield"] << "  dropped=" << n_bad_fit
                << std::endl;
    }

    if ( toys.empty() ) {
      std::cout << "==> WARNING: no surviving toys at true Nsig = " << trueNsig << std::endl;
      continue;
    }

    // Sort toys by Delta NLL (Feldman-Cousins rank ordering according to likelihood ratio R)
    std::sort( toys.begin(), toys.end(),
               []( const ToyResult& a, const ToyResult& b ) { return a.delta_nll < b.delta_nll; } );

    // Accept top (cl * N_toys) into the confidence interval slice
    size_t nAccept = static_cast<size_t>( cl * toys.size() );
    if ( nAccept < 1 ) nAccept = 1;

    for ( size_t k = 0; k < nAccept; ++k ) { bandHist->Fill( toys[k].nsig_hat, trueNsig, 1.0 / double( nAccept ) ); }

    double lo = 0.;
    double hi = 0.;

    if ( beltMode == 0 ) {

      lo = toys[0].nsig_hat;
      hi = toys[0].nsig_hat;
      for ( size_t k = 0; k < nAccept; ++k ) {
        lo = std::min( lo, toys[k].nsig_hat );
        hi = std::max( hi, toys[k].nsig_hat );
      }

    } else {

      std::vector<double> hats;
      hats.reserve( toys.size() );
      for ( const auto& tr : toys ) { hats.push_back( tr.nsig_hat ); }
      std::sort( hats.begin(), hats.end() );

      const size_t nKeep = static_cast<size_t>( std::ceil( cl * double( hats.size() ) ) );
      auto         iv    = shortestInterval( hats, nKeep );
      lo                 = iv.first;
      hi                 = iv.second;
    }

    belt_y.push_back( trueNsig );
    belt_lo.push_back( lo );
    belt_hi.push_back( hi );
  }

  if ( !belt_y.empty() ) {
    const size_t N = belt_y.size();
    m_belt         = new TGraph( 2 * N + 1 );
    m_belt->SetName( "g_FC_belt" );
    for ( size_t k = 0; k < N; ++k ) { m_belt->SetPoint( k, belt_lo[k], belt_y[k] ); }
    for ( size_t k = 0; k < N; ++k ) { m_belt->SetPoint( N + k, belt_hi[N - 1 - k], belt_y[N - 1 - k] ); }
    m_belt->SetPoint( 2 * N, belt_lo[0], belt_y[0] ); // close it
  }

  if ( m_file ) {
    m_file->cd();
    bandHist->Write();
    if ( m_belt ) m_belt->Write();
  }

  RooMsgService::instance().reset();
  return bandHist;
}
