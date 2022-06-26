#include <iostream>
#include <vector>
#include <string>

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp> 
#include <boost/filesystem.hpp>

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
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
  auto ev_vec = Momentum_recon::ReadEventInformationBin(momchfilename);

  // output file
  TFile *outputfile = new TFile((TString)outputfilename, "recreate");
  TTree *otree = new TTree("tree", "tree");
  
  BOOST_LOG_TRIVIAL(info) << "Output filename : " << outputfilename;

  double mom_entry_total[10] = {};
  double ang_entry_total[10] = {};
  double mom_ang_entry_total[10][10] = {};
  double mom_total[10] = {};
  double ang_total[10] = {};
  double mom_ang_total[10][10] = {};
  double mom_eff[10] = {};
  double ang_eff[10] = {};
  double mom_ang_eff[10][10] = {};

  otree->Branch("mom_entry_total", mom_entry_total, "mom_entry_total[10]/D");
  otree->Branch("ang_entry_total", ang_entry_total, "ang_entry_total[10]/D");
  otree->Branch("mom_ang_entry_total", mom_ang_entry_total, "mom_ang_entry_total[10][10]/D");
  otree->Branch("mom_total", mom_total, "mom_total[10]/D");
  otree->Branch("ang_total", ang_total, "ang_total[10]/D");
  otree->Branch("mom_ang_total", mom_ang_total, "mom_ang_total[10][10]/D");
  otree->Branch("mom_eff", mom_eff, "mom_eff[10]/D");
  otree->Branch("ang_eff", ang_eff, "ang_eff[10]/D");
  otree->Branch("mom_ang_eff", mom_ang_eff, "mom_ang_eff[10][10]/D");

  double true_mu_mom = -1;
  double true_mu_ang = -1;
  int true_mu_id = -1;

  int imom = -1;
  int iang = -1;

  for ( auto ev : ev_vec ) {

    BOOST_LOG_TRIVIAL(debug) << "Entry : " << ev.groupid;
    bool mu_detect_flag = false;

    for ( auto chain : ev.true_chains ) {
      if ( chain.particle_flag == 13 &&
	   chain.direction == 1 ) {
	true_mu_mom = chain.bm_range_mom;
	imom = true_mu_mom / 100.;
	if ( imom > 9 ) imom = 9;
	double ax = chain.base.back().ax;
	double ay = chain.base.back().ay;
	true_mu_ang = std::atan(std::hypot(ax, ay)) * TMath::RadToDeg();
	iang = true_mu_ang / 5.;
	if ( iang > 9 ) iang = 9;
	true_mu_id = chain.chainid;
	BOOST_LOG_TRIVIAL(debug) << "True momentm : " << true_mu_mom << " [MeV/c], "
				 << "True angle : " << true_mu_ang << " [deg] , "
				 << "True mu id : " << true_mu_id << ", "
				 << "Bin id : " << imom << ", " << iang;

	break;
      }
    }


    for ( auto chain : ev.chains ) {
      int recon_particle_id = chain.particle_flag % 10000;
      if ( chain.chainid == true_mu_id && // same track
	   recon_particle_id == 13 && // muon
	   ev.vertex_pl / 1000 == ev.vertex_pl % 1000 ) { // same vertex
	mu_detect_flag = true;
	break;
      }
    }


    if ( mu_detect_flag ) {
      BOOST_LOG_TRIVIAL(debug) << "Muon detect flag true";
      mom_eff[imom] += ev.weight;
      ang_eff[iang] += ev.weight;
      mom_ang_eff[imom][iang] += ev.weight;
    }
    else {
      BOOST_LOG_TRIVIAL(debug) << "Muon detect flag false";
    }
    mom_total[imom] += ev.weight;
    ang_total[iang] += ev.weight;
    mom_ang_total[imom][iang] += ev.weight;
    mom_entry_total[imom] += 1.;
    ang_entry_total[iang] += 1.;
    mom_ang_entry_total[imom][iang] += 1.;
    
  }

  double mom_bin_edge[11] = {0., 100., 200., 300., 400., 500., 600., 700., 800., 900., 1500.};
  double ang_bin_edge[11] = {0., 5., 10., 15., 20., 25., 30., 35., 40., 45., 90.};
  TH1D *hist_eff_muon_mom = new TH1D("hist_eff_muon_mom", "Muon efficiency;Momentum [MeV/c];Efficiency",
				  10, mom_bin_edge);
  TH1D *hist_eff_muon_ang = new TH1D("hist_eff_muon_ang", "Muon efficiency;Angle [deg];Efficiency",
				     10, ang_bin_edge);
  TH2D *hist_eff_muon_mom_ang = new TH2D("hist_eff_muon_mom_ang", "Muon efficiency;Momentum [MeV/c];Angle[deg]",
					 10, mom_bin_edge, 10, ang_bin_edge);
  hist_eff_muon_mom->GetYaxis()->SetRangeUser(0., 1.);
  hist_eff_muon_ang->GetYaxis()->SetRangeUser(0., 1.);
  hist_eff_muon_mom_ang->GetZaxis()->SetRangeUser(0., 1.);

  for ( int i = 0; i < 10; i++ ) {
    double imom = (mom_bin_edge[i] + mom_bin_edge[i+1] / 2.);
    double iang = (ang_bin_edge[i] + ang_bin_edge[i+1] / 2.);
    hist_eff_muon_mom->Fill(imom, mom_eff[i]/mom_total[i]);
    hist_eff_muon_ang->Fill(iang, ang_eff[i]/ang_total[i]);
    for ( int j = 0; j < 10; j++ ) {
      double jang = (ang_bin_edge[j] + ang_bin_edge[j+1] / 2.);
      hist_eff_muon_mom_ang->Fill(imom, jang, mom_ang_eff[i][j] / mom_ang_total[i][j]);
    }
  }
  
  outputfile->cd();
  hist_eff_muon_mom->Write();
  hist_eff_muon_ang->Write();
  hist_eff_muon_mom_ang->Write();
  otree->Fill();
  otree->Write();
  outputfile->Close();
}
