#include <iostream>
#include <vector>
#include <cmath>
#include <string>

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/filesystem.hpp>

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>

#include <B2Reader.hh>

#include <McsClass.hpp>

namespace logging = boost::log;
namespace fs = boost::filesystem;

// program to evaluate the number of sand muons detected in ECC5

int main (int argc, char* argv[]) {

  logging::core::get()->set_filter
    (
     logging::trivial::severity >= logging::trivial::info
     );

  BOOST_LOG_TRIVIAL(info) << "==========Wall Check==========";

  if ( argc != 3 ) {
    BOOST_LOG_TRIVIAL(error) << "Usage : " << argv[0]
			     << " <input momch file> <output file>";
    std::exit(1);
  }

  std::string momchfilename = argv[1];
  if ( !fs::exists(momchfilename) ) std::exit(1);
  auto ev_vec = Momentum_recon::ReadEventInformationBin(momchfilename);

  TFile *ofile = new TFile((TString)argv[2], "recreate");
  TH1D *hist_multi = new TH1D("hist_multi", ";# of tracks;Entries", 5, 0.5, 5.5);

  for ( auto ev : ev_vec ) {
    
    if ( ev.chains.empty() ) continue;

    if ( ev.vertex_material == -2 ) {
      hist_multi->Fill(ev.chains.size(), ev.weight);
    }

  }

  ofile->cd();
  hist_multi->Write();
  ofile->Close();
  std::exit(0);

}
