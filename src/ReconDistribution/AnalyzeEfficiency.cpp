#include <iostream>
#include <vector>
#include <string>

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp> 
#include <boost/filesystem.hpp>

#include <TFile.h>
#include <TGraphErrors.h>
#include <TH2D.h>

#include <B2Reader.hh>

#include <McsClass.hpp>

#include "AnalyzeEfficiency.hpp"

namespace logging = boost::log;
namespace fs = boost::filesystem;

void AnalyzeEfficiency(std::string b2filename,
		       std::string momchfilename,
		       std::string outputfilename) {

  BOOST_LOG_TRIVIAL(info) << "==========Efficiency mode==========";

  // input B2 file
  B2Reader reader(b2filename);

  // input momch file
  if ( !fs::exists(momchfilename) ) {
    throw std::runtime_error("File not found : " + momchfilename);   
  }
  auto ev_vec = Momentum_recon::ReadEventInfomationBin(momchfilename);

  // output file
  TFile *outputfile = new TFile((TString)outputfilename, "recreate");
  BOOST_LOG_TRIVIAL(info) << "Output filename : " << outputfilename;

  double mom_bin[10];
  double ang_bin[10];

  for ( auto ev : ev_vec ) {
    
  }

  TGraphErrors *g_eff_muon_mom;
  TGraphErrors *g_eff_muon_ang;
  TH2D *hist_eff_muon_mom_ang;
  

}
