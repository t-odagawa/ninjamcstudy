// system includes
#include <sstream>

// boost includes
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

// ROOT includes
#include <TFile.h>
#include <TH2D.h>

#include <B2Reader.hh>
#include <McsClass.hpp>

namespace logging = boost::log;

int main (int argc, char* argv[]) {

  logging::core::get()->set_filter
    (
     logging::trivial::severity >= logging::trivial::info
     );

  
  BOOST_LOG_TRIVIAL(info) << "==========Start==========";

  try {

    // input B2 file
    std::stringstream b2filename_ss;
    b2filename_ss << argv[1];
    B2Reader reader(b2filename_ss.str());
    
    // input momch file
    auto ev_vec = Momentum_recon::ReadEventInformationBin(argv[2]);
    
    // output file
    std::stringstream ofilename_ss;
    ofilename_ss << argv[3];
    TFile *ofile = new TFile(ofilename_ss.str().c_str(), "recreate");
    
    TH2D *hist_muon_mom_recon_true_mcs_stop = new TH2D("hist_muon_mom_recon_true_mcs_stop",
						       "Muon MCS momentum vs true (stop only);p_{true} [MeV/c];p_{MCS} [MeV/c]",
						       50, 0., 2000., 50, 0., 2000.);
    TH2D *hist_muon_mom_recon_true_mcs_all = new TH2D("hist_muon_mom_recon_true_mcs_all",
						      "Muon MCS momentum vs true (no stop cut);p_{true} [MeV/c];p_{MCS} [MeV/c]",
						      50, 0., 2000., 50, 0., 2000.);
    TH2D *hist_muon_mom_recon_range_mcs_stop = new TH2D("hist_muon_mom_recon_range_mcs_stop",
							"Muon MCS vs rane momentum (stop only);p_{range} [MeV/c];p_{MCS} [MeV/c]",
							50, 0., 2000., 50, 0., 2000.);

    for ( auto ev : ev_vec ) {
      
      reader.ReadSpill(ev.groupid);
      
      if ( ev.vertex_material != -2 ) continue; // penetrate only

      if ( ev.chains.empty() ) continue;

      bool muon_correct_id_flag = false;

      for ( auto chain : ev.chains ) {
	if ( chain.particle_flag % 10000 == 13 &&
	     chain.particle_flag / 10000 == 13 ) muon_correct_id_flag = true;
      }

      if ( !muon_correct_id_flag ) continue;

      for ( auto chain : ev.chains ) {
	if ( chain.particle_flag% 10000 != 13 ) continue;

	double true_momentum = -1.;
	for ( auto true_chain : ev.true_chains ) {
	  if ( chain.chainid == true_chain.chainid ) {
	    true_momentum = true_chain.bm_range_mom;
	    break;
	  }
	}

	if ( chain.stop_flag == 1 ) {
	  hist_muon_mom_recon_true_mcs_stop->Fill(true_momentum, chain.ecc_mcs_mom[0], ev.weight);
	  hist_muon_mom_recon_range_mcs_stop->Fill(chain.bm_range_mom, chain.ecc_mcs_mom[0], ev.weight);
	}
	hist_muon_mom_recon_true_mcs_all->Fill(true_momentum, chain.ecc_mcs_mom[0], ev.weight);
	
      } // chain

    } // ev

    ofile->cd();
    hist_muon_mom_recon_true_mcs_stop->Write();
    hist_muon_mom_recon_range_mcs_stop->Write();
    hist_muon_mom_recon_true_mcs_all->Write();
    ofile->Close();

  } catch ( const std::runtime_error &error ) {
    BOOST_LOG_TRIVIAL(fatal) << "Runtime error : " << error.what();
    std::exit(1);
  } catch ( const std::invalid_argument &error ) {
    BOOST_LOG_TRIVIAL(fatal) << "Invalid argument error : " << error.what();
    std::exit(1);
  }

  BOOST_LOG_TRIVIAL(info) << "==========Finish==========";
  std::exit(0);

}
