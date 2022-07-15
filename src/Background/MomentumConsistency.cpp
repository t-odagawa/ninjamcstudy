#include <iostream>
#include <vector>
#include <cmath>
#include <string>

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/filesystem.hpp>

#include <TFile.h>
#include <TTree.h>

#include <McsClass.hpp>

namespace fs = boost::filesystem;

int main(int argc, char* argv[]) {

  if ( argc != 3 ) {
    std::cerr << "Usage : " << argv[0]
	      << " <input momch file w/ PID> <output root file>" << std::endl;
    std::exit(1);
  }

  if ( !fs::exists((std::string)argv[1]) ) {
    std::cerr << "File does not exists : " << argv[1] << std::endl;
    std::exit(1);
  }

  auto ev_vec = Momentum_recon::ReadEventInformationBin((std::string)argv[1]);

  TFile *ofile = new TFile((TString)argv[2], "recreate");
  TTree *otree = new TTree("tree", "tree");

  double weight;
  int true_particle_id;
  int recon_particle_id;
  int stop_flag;
  double angle;
  double true_momentum;
  double recon_mom_bm_range;
  double recon_mom_ecc_mcs_muon;
  double recon_mom_ecc_range_proton;
  double recon_mom_ecc_mcs_proton;
  double recon_mom_ecc_range_pion;
  double recon_mom_ecc_mcs_pion;
  double recon_mom_bm_range_err[2];
  double recon_mom_ecc_mcs_muon_err[2];
  double recon_mom_ecc_range_proton_err[2];
  double recon_mom_ecc_mcs_proton_err[2];
  double recon_mom_ecc_range_pion_err[2];
  double recon_mom_ecc_mcs_pion_err[2];

  otree->Branch("weight", &weight, "weight/D");
  otree->Branch("true_particle_id", &true_particle_id, "true_particle_id/I");
  otree->Branch("recon_particle_id", &recon_particle_id, "recon_particle_id/I");
  otree->Branch("stop_flag", &stop_flag, "stop_flag/I");
  otree->Branch("angle", &angle, "angle/D");
  otree->Branch("true_momentum", &true_momentum, "true_momentum/D");
  otree->Branch("recon_mom_bm_range", &recon_mom_bm_range, "recon_mom_bm_range/D");
  otree->Branch("recon_mom_ecc_mcs_muon", &recon_mom_ecc_mcs_muon, "recon_mom_ecc_mcs_muon/D");
  otree->Branch("recon_mom_ecc_range_proton", &recon_mom_ecc_range_proton, "recon_mom_ecc_range_proton/D");
  otree->Branch("recon_mom_ecc_mcs_proton", &recon_mom_ecc_mcs_proton, "recon_mom_ecc_mcs_proton/D");
  otree->Branch("recon_mom_ecc_range_pion", &recon_mom_ecc_range_pion, "recon_mom_ecc_range_pion/D");
  otree->Branch("recon_mom_ecc_mcs_pion", &recon_mom_ecc_mcs_pion, "recon_mom_ecc_mcs_pion/D");
  otree->Branch("recon_mom_bm_range_err", recon_mom_bm_range_err, "recon_mom_bm_range_err[2]/D");
  otree->Branch("recon_mom_ecc_mcs_muon_err", recon_mom_ecc_mcs_muon_err, "recon_mom_ecc_mcs_muon_err[2]/D");
  otree->Branch("recon_mom_ecc_range_proton_err", recon_mom_ecc_range_proton_err, "recon_mom_ecc_range_proton_err[2]/D");
  otree->Branch("recon_mom_ecc_mcs_proton_err", recon_mom_ecc_mcs_proton_err, "recon_mom_ecc_mcs_proton_err[2]/D");
  otree->Branch("recon_mom_ecc_range_pion_err", recon_mom_ecc_range_pion_err, "recon_mom_ecc_range_pion_err[2]/D");
  otree->Branch("recon_mom_ecc_mcs_pion_err", recon_mom_ecc_mcs_pion_err, "recon_mom_ecc_mcs_pion_err[2]/D");
  
  for ( auto ev : ev_vec ) {
    
    if ( ev.chains.empty() ) continue;

    for ( auto chain : ev.chains ) {
      weight = ev.weight;
      true_particle_id = chain.particle_flag / 10000;
      recon_particle_id = chain.particle_flag % 10000;
      stop_flag = chain.stop_flag;
      
      for ( auto true_chain : ev.true_chains ) {
	if ( true_chain.chainid == chain.chainid ) {
	  true_momentum = true_chain.bm_range_mom;
	  if ( true_chain.direction == 1 ) {
	    angle = std::hypot(true_chain.base.back().ax,
			       true_chain.base.back().ay);
	  }
	  else if ( true_chain.direction == -1 ) {
	    angle = std::hypot(true_chain.base.front().ax,
			       true_chain.base.front().ay);
	  }
	  break;
	}
      }
      
      recon_mom_bm_range = chain.bm_range_mom;
      recon_mom_ecc_mcs_muon = chain.ecc_mcs_mom[0];
      recon_mom_ecc_range_proton = chain.ecc_range_mom[1];
      recon_mom_ecc_mcs_proton = chain.ecc_mcs_mom[1];
      recon_mom_ecc_range_pion = chain.ecc_range_mom[0];
      recon_mom_ecc_mcs_pion = chain.bm_curvature_mom;
      for (int i = 0; i < 2; i++) {
	recon_mom_bm_range_err[i] = chain.bm_range_mom_error[i];
	recon_mom_ecc_mcs_muon_err[i] = chain.ecc_mcs_mom_error[0][i];
	recon_mom_ecc_range_proton_err[i] = chain.ecc_range_mom_error[1][i];
	recon_mom_ecc_mcs_proton_err[i] = chain.ecc_mcs_mom_error[1][i];
	recon_mom_ecc_range_pion_err[i] = chain.ecc_range_mom_error[0][i];
	recon_mom_ecc_mcs_pion_err[i] = chain.bm_curvature_mom_error[i];
      }

      otree->Fill();
    }
  }

  ofile->cd();
  otree->Write();
  ofile->Close();

  std::exit(0);

}
