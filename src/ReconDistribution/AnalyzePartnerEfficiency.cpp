#include <iostream>
#include <vector>
#include <string>
#include <cmath>

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/filesystem.hpp>

#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>
#include <TMath.h>

#include <B2Reader.hh>

#include <McsClass.hpp>

#include "AnalyzePartnerEfficiency.hpp"

namespace logging = boost::log;
namespace fs = boost::filesystem;

void AnalyzePartnerEfficiency(std::string b2filename,
			      std::string momchfilename,
			      std::string outputfilename) {
  
  BOOST_LOG_TRIVIAL(info) << "==========Partner efficiency mode==========";

  // input B2 file
  B2Reader reader(b2filename);

  // input momch file
  if ( !fs::exists(momchfilename) ) {
    throw std::runtime_error("File not found : " + momchfilename);
  }
  auto ev_vec = Momentum_recon::ReadEventInformationBin(momchfilename);

  // output file
  TFile *outputfile = new TFile((TString)outputfilename, "recreate");
  TH1D *true_proton_mom_hist = new TH1D("true_proton_mom_hist", "", 100, 0., 1500.);
  TH1D *true_proton_ang_hist = new TH1D("true_proton_ang_hist", "", 60, 0., 180.);
  TH1D *true_pion_mom_hist = new TH1D("true_pion_mom_hist", "", 100, 0., 1500.);
  TH1D *true_pion_ang_hist = new TH1D("true_pion_ang_hist", "", 60, 0., 180.);
  TTree *otree = new TTree("tree", "tree");

  BOOST_LOG_TRIVIAL(info) << "Output filename : " << outputfilename;

  double mom_bin_edge[12] = {0., 100., 200., 300., 400., 500., 600., 700., 800., 900., 1000., 1500.};
  double ang_bin_edge[37] = {0., 5., 10., 15., 20., 25., 
			     30., 35., 40., 45., 50., 55.,
			     60., 65., 70., 75., 80., 85.,
			     90., 95., 100., 105., 110., 115.,
			     120., 125., 130., 135., 140., 145.,
			     150., 155., 160., 165., 170., 175.,
			     180.};
  
  double partner_mom_entry_total[11] = {};
  double partner_ang_entry_total[36] = {};
  double partner_mom_ang_entry_total[11][36] = {{}};
  double partner_mom_total[11] = {};
  double partner_ang_total[36] = {};
  double partner_mom_ang_total[11][36] = {{}};
  double partner_mom_eff[11] = {};
  double partner_ang_eff[36] = {};
  double partner_mom_ang_eff[11][36] = {{}};

  double proton_mom_entry_total[11] = {};
  double proton_ang_entry_total[36] = {};
  double proton_mom_ang_entry_total[11][36] = {{}};
  double proton_mom_total[11] = {};
  double proton_ang_total[36] = {};
  double proton_mom_ang_total[11][36] = {{}};
  double proton_mom_eff[11] = {};
  double proton_ang_eff[36] = {};
  double proton_mom_ang_eff[11][36] = {{}};

  double pion_mom_entry_total[11] = {};
  double pion_ang_entry_total[36] = {};
  double pion_mom_ang_entry_total[11][36] = {{}};
  double pion_mom_total[11] = {};
  double pion_ang_total[36] = {};
  double pion_mom_ang_total[11][36] = {{}};
  double pion_mom_eff[11] = {};
  double pion_ang_eff[36] = {};
  double pion_mom_ang_eff[11][36] = {{}};

  otree->Branch("partner_mom_entry_total", partner_mom_entry_total, "partner_mom_entry_total[11]/D");
  otree->Branch("partner_ang_entry_total", partner_ang_entry_total, "partner_ang_entry_total[36]/D");
  otree->Branch("partner_mom_ang_entry_total", partner_mom_ang_entry_total, "partner_mom_ang_entry_total[11][36]/D");
  otree->Branch("partner_mom_total", partner_mom_total, "partner_mom_total[11]/D");
  otree->Branch("partner_ang_total", partner_ang_total, "partner_ang_total[36]/D");
  otree->Branch("partner_mom_ang_total", partner_mom_ang_total, "partner_mom_ang_total[11][36]/D");
  otree->Branch("partner_mom_eff", partner_mom_eff, "partner_mom_eff[11]/D");
  otree->Branch("partner_ang_eff", partner_ang_eff, "partner_ang_eff[36]/D");
  otree->Branch("partner_mom_ang_eff", partner_mom_ang_eff, "partner_mom_ang_eff[11][36]/D");
  otree->Branch("proton_mom_entry_total", proton_mom_entry_total, "proton_mom_entry_total[11]/D");
  otree->Branch("proton_ang_entry_total", proton_ang_entry_total, "proton_ang_entry_total[36]/D");
  otree->Branch("proton_mom_ang_entry_total", proton_mom_ang_entry_total, "proton_mom_ang_entry_total[11][36]/D");
  otree->Branch("proton_mom_total", proton_mom_total, "proton_mom_total[11]/D");
  otree->Branch("proton_ang_total", proton_ang_total, "proton_ang_total[36]/D");
  otree->Branch("proton_mom_ang_total", proton_mom_ang_total, "proton_mom_ang_total[11][36]/D");
  otree->Branch("proton_mom_eff", proton_mom_eff, "proton_mom_eff[11]/D");
  otree->Branch("proton_ang_eff", proton_ang_eff, "proton_ang_eff[36]/D");
  otree->Branch("proton_mom_ang_eff", proton_mom_ang_eff, "proton_mom_ang_eff[11][36]/D");
  otree->Branch("pion_mom_entry_total", pion_mom_entry_total, "pion_mom_entry_total[11]/D");
  otree->Branch("pion_ang_entry_total", pion_ang_entry_total, "pion_ang_entry_total[36]/D");
  otree->Branch("pion_mom_ang_entry_total", pion_mom_ang_entry_total, "pion_mom_ang_entry_total[11][36]/D");
  otree->Branch("pion_mom_total", pion_mom_total, "pion_mom_total[11]/D");
  otree->Branch("pion_ang_total", pion_ang_total, "pion_ang_total[36]/D");
  otree->Branch("pion_mom_ang_total", pion_mom_ang_total, "pion_mom_ang_total[11][36]/D");
  otree->Branch("pion_mom_eff", pion_mom_eff, "pion_mom_eff[11]/D");
  otree->Branch("pion_ang_eff", pion_ang_eff, "pion_ang_eff[36]/D");
  otree->Branch("pion_mom_ang_eff", pion_mom_ang_eff, "pion_mom_ang_eff[11][36]/D");



  for ( auto ev : ev_vec ) {
    
    // good muon id and good vertex id
    int true_mu_id = -1;
    bool mu_detect_flag = false;
    for ( auto true_chain : ev.true_chains ) {
      if ( true_chain.particle_flag == 13 &&
	   true_chain.direction == 1 ) {
	true_mu_id = true_chain.chainid;
      }
    }

    for ( auto chain : ev.chains ) {
      int recon_particle_id = chain.particle_flag % 10000;
      if ( chain.chainid == true_mu_id &&
	   recon_particle_id == 13 &&
	   ev.vertex_pl / 1000 == ev.vertex_pl % 1000 ) {
	mu_detect_flag = true;
	break;
      }
    }

    if ( !mu_detect_flag ) continue;

    int true_id = -1;
    double true_mom = -1.;
    double true_ang = -1.;
    double ax, ay;
    int imom = -1;
    int iang = -1;
    for ( auto true_chain : ev.true_chains ) {
      
      if ( !(true_chain.particle_flag == 2212) &&
	   !(std::abs(true_chain.particle_flag) == 211) ) continue;
      
      true_mom = true_chain.bm_range_mom;
      imom = true_mom / 100.;
      if ( imom >= 10 ) imom = 10;
      if ( true_chain.direction == 1 ) {
	ax = true_chain.base.back().ax;
	ay = true_chain.base.back().ay;
      }
      else if ( true_chain.direction == -1 ) {
	ax = true_chain.base.front().ax;
	ay = true_chain.base.front().ay;
      }
      true_ang = std::atan(std::hypot(ax, ay)) * TMath::RadToDeg();
      if ( true_chain.direction == -1 ) {
	true_ang = -true_ang + 180.;
      }
      iang = true_ang / 5.;

      true_id = true_chain.chainid;

      partner_mom_entry_total[imom] += 1.;
      partner_ang_entry_total[iang] += 1.;
      partner_mom_ang_entry_total[iang][imom] += 1.;
      partner_mom_total[imom] += ev.weight;
      partner_ang_total[iang] += ev.weight;
      partner_mom_ang_total[imom][iang] += ev.weight;

      if ( true_chain.particle_flag == 2212 ) { // proton
	true_proton_mom_hist->Fill(true_mom, ev.weight);
	true_proton_ang_hist->Fill(true_ang, ev.weight);

	proton_mom_entry_total[imom] += 1.;
	proton_ang_entry_total[iang] += 1.;
	proton_mom_ang_entry_total[imom][iang] += 1.;
	proton_mom_total[imom] += ev.weight;
	proton_ang_total[iang] += ev.weight;
	proton_mom_ang_total[imom][iang] += ev.weight;

      } else if ( std::abs(true_chain.particle_flag) == 211 ) { // pion
	true_pion_mom_hist->Fill(true_mom, ev.weight);
	true_pion_ang_hist->Fill(true_ang, ev.weight);

	pion_mom_entry_total[imom] += 1.;
	pion_ang_entry_total[iang] += 1.;
	pion_mom_ang_entry_total[imom][iang] += 1.;
	pion_mom_total[imom] += ev.weight;
	pion_ang_total[iang] += ev.weight;
	pion_mom_ang_total[imom][iang] += ev.weight;

      }

      for ( auto chain : ev.chains ) {
	int recon_particle_id = chain.particle_flag % 10000;
	if ( chain.chainid == true_id ) { // same track is reconstructed
	  partner_mom_eff[imom] += ev.weight;
	  partner_ang_eff[iang] += ev.weight;
	  partner_mom_ang_eff[imom][iang] += ev.weight;

	  double recon_mom;
	  if ( recon_particle_id == 2212 &&
	       chain.stop_flag == 2 ) {
	    recon_mom = chain.ecc_range_mom[1];
	  }
	  else if ( recon_particle_id == 2212 &&
		    chain.stop_flag == 0 ) {
	    recon_mom = chain.ecc_mcs_mom[1];
	  }
	  else if ( recon_particle_id == 211 &&
		    chain.stop_flag == 2 ) {
	    recon_mom = chain.ecc_range_mom[0];	   
	  }
	  else if ( recon_particle_id == 211 &&
		    chain.stop_flag == 0 ) {
	    recon_mom = chain.bm_curvature_mom;
	  }

	  if ( recon_particle_id == true_chain.particle_flag &&
	       recon_particle_id == 2212 ) {
	    proton_mom_eff[imom] += ev.weight;
	    proton_ang_eff[iang] += ev.weight;
	    proton_mom_ang_eff[imom][iang] += ev.weight;
	  }
	  else if ( recon_particle_id == true_chain.particle_flag &&
		    recon_particle_id == 211 ) {
	    pion_mom_eff[imom] += ev.weight;
	    pion_ang_eff[iang] += ev.weight;
	    pion_mom_ang_eff[imom][iang] += ev.weight;
	  }      
	}
      }
    }
  }

  outputfile->cd();
  true_proton_mom_hist->Write();
  true_proton_ang_hist->Write();
  true_pion_mom_hist->Write();
  true_pion_ang_hist->Write();
  otree->Fill();
  otree->Write();
  outputfile->Close();
  
}
