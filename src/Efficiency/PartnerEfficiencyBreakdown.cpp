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
#include <B2Enum.hh>
#include <B2TrackSummary.hh>
#include <B2HitSummary.hh>
#include <B2EventSummary.hh>
#include <B2VertexSummary.hh>

#include <McsClass.hpp>

namespace logging = boost::log;
namespace fs = boost::filesystem;

int main (int argc, char* argv[]) {

  logging::core::get()->set_filter
    (
     logging::trivial::severity >= logging::trivial::info
     );

  if ( argc != 4 ) {
    BOOST_LOG_TRIVIAL(error) << "Usage : " << argv[0]
			     << " <B2 file> <momch file w/ BM> <outputfile>";
    std::exit(1);
  }

  BOOST_LOG_TRIVIAL(info) << "==========Partner efficiency breakdown plots creation start==========";

  // input B2 file
  std::string inputfilename = argv[1];
  B2Reader reader(inputfilename);

  // input momch file
  std::string momchfilename = argv[2];
  auto ev_vec = Momentum_recon::ReadEventInformationBin(momchfilename);

  // output file
  std::string outputfilename = argv[3];
  TFile *ofile = new TFile((TString)outputfilename, "recreate");
  TTree *otree = new TTree("tree", "tree");

  TH1D *hist_true_proton_mom;
  TH1D *hist_true_proton_mom_cut1; // One iron plate
  TH1D *hist_true_proton_mom_cut2; // Partner attach (MD)
  TH1D *hist_true_proton_mom_cut3; // PID

  TH1D *hist_true_proton_ang;
  TH1D *hist_true_proton_ang_cut1;
  TH1D *hist_true_proton_ang_cut2;
  TH1D *hist_true_proton_ang_cut3;

  TH1D *hist_true_pion_mom;
  TH1D *hist_true_pion_mom_cut1;
  TH1D *hist_true_pion_mom_cut2;
  TH1D *hist_true_pion_mom_cut3;

  TH1D *hist_true_pion_ang;
  TH1D *hist_true_pion_ang_cut1;
  TH1D *hist_true_pion_ang_cut2;
  TH1D *hist_true_pion_ang_cut3;

  double proton_mom_entry_total[11] = {};
  double proton_ang_entry_total[36] = {};
  double proton_mom_ang_entry_total[11][36] = {};
  double proton_mom_cut0[11] = {};
  double proton_mom_cut1[11] = {};
  double proton_mom_cut2[11] = {};
  double proton_mom_cut3[11] = {};
  double proton_ang_cut0[36] = {};
  double proton_ang_cut1[36] = {};
  double proton_ang_cut2[36] = {};
  double proton_ang_cut3[36] = {};
  double proton_mom_ang_cut0[11][36] = {{}};
  double proton_mom_ang_cut1[11][36] = {{}};
  double proton_mom_ang_cut2[11][36] = {{}};
  double proton_mom_ang_cut3[11][36] = {{}};

  double pion_mom_entry_total[11] = {};
  double pion_ang_entry_total[36] = {};
  double pion_mom_ang_entry_total[11][36] = {};
  double pion_mom_cut0[11] = {};
  double pion_mom_cut1[11] = {};
  double pion_mom_cut2[11] = {};
  double pion_mom_cut3[11] = {};
  double pion_ang_cut0[36] = {};
  double pion_ang_cut1[36] = {};
  double pion_ang_cut2[36] = {};
  double pion_ang_cut3[36] = {};
  double pion_mom_ang_cut0[11][36] = {{}};
  double pion_mom_ang_cut1[11][36] = {{}};
  double pion_mom_ang_cut2[11][36] = {{}};
  double pion_mom_ang_cut3[11][36] = {{}};

  otree->Branch("proton_mom_entry_total", proton_mom_entry_total, "proton_mom_entry_total[11]/D");
  otree->Branch("proton_ang_entry_total", proton_ang_entry_total, "proton_ang_entry_total[36]/D");
  otree->Branch("proton_mom_ang_entry_total", proton_mom_ang_entry_total, "proton_mom_ang_entry_total[11][36]/D");
  otree->Branch("proton_mom_cut0", proton_mom_cut0, "proton_mom_cut0[11]/D");

  for ( auto ev : ev_vec ) {
    
    BOOST_LOG_TRIVIAL(debug) << "Entry : " << ev.groupid;

    // When the event is correctly detected as CC event
    int true_mu_id = -1;
    bool mu_detect_flag = false;
    for ( auto true_chain : ev.true_chains ) {
      if ( true_chain.particle_flag == 13 &&
	   true_chain.direction == 1 ) {
	true_mu_id = true_chain.chainid;
	break;
      }
    }

    for ( auto chain : ev.chains ) {
      int recon_particle_id = chain.particle_flag % 10000;
      if ( chain.chainid == true_mu_id &&
	   recon_particle_id == 13 && // true mu is id-ed as muon
	   ev.vertex_pl / 1000 == ev.vertex_pl % 1000 && // vertex is correct
	   ev.vertex_material >= 0 ) { // start from FV
	mu_detect_flag = true;
	break;
      }
    }
    
    if ( !mu_detect_flag ) continue;

    for ( auto true_chain : ev.true_chains ) {
      
      if ( !(true_chain.particle_flag == 2212) &&
	   !(std::abs(true_chain.particle_flag) == 221) ) continue;

      bool detect_flag_cut1 = false;
      bool detect_flag_cut2 = false;
      bool detect_flag_cut3 = false;

      double true_mom = -1.;
      double true_ang = -1.;
      double ax, ay;
      int true_id = -1;
      
      int imom = -1;
      int iang = -1;

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

      // Penetrate at least one iron plate
      if ( true_chain.base.size() == 2 &&
	   true_chain.base.front().pl <= 15 ) {
	detect_flag_cut1 = true;
      }
      else if ( true_chain.base.size() == 2 &&
		true_chain.base.front().pl % 2 == 0 ) {
	detect_flag_cut1 = true;
      }
      else if ( true_chain.base.size() > 2 ) {
	detect_flag_cut1 = true;
      }


      if ( true_chain.particle_flag == 2212 ) { // proton
	proton_mom_entry_total[imom] += 1.;
	proton_ang_entry_total[iang] += 1.;
	proton_mom_ang_entry_total[imom][iang] += 1.;
	proton_mom_cut0[imom] += ev.weight;
	proton_ang_cut0[iang] += ev.weight;
	proton_mom_ang_cut0[imom][iang] += ev.weight;
	if ( detect_flag_cut1 ) {
	  proton_mom_cut1[imom] += ev.weight;
	  proton_ang_cut1[iang] += ev.weight;
	  proton_mom_ang_cut1[imom][iang] += ev.weight;	  
	}
      }
      else if ( std::abs(true_chain.particle_flag) == 211 ) { // pion
	pion_mom_entry_total[imom] += 1.;
	pion_ang_entry_total[iang] += 1.;
	pion_mom_ang_entry_total[imom][iang] += 1.;
	pion_mom_cut0[imom] += ev.weight;
	pion_ang_cut0[iang] += ev.weight;
	pion_mom_ang_cut0[imom][iang] += ev.weight;
	if ( detect_flag_cut1 ) {
	  pion_mom_cut1[imom] += ev.weight;
	  pion_ang_cut1[iang] += ev.weight;
	  pion_mom_ang_cut1[imom][iang] += ev.weight;	  
	}
      }

      // Attach as the partner
      for ( auto chain : ev.chains ) {
	if ( chain.chainid == true_id ) {
	  detect_flag_cut2 = true;
	  // PID is correct
	  int recon_particle_id = chain.particle_flag % 10000;
	  if ( recon_particle_id == true_chain.particle_flag ) {
	    detect_flag_cut3 = true;
	  }
	  break;
	}
      }
      
      // ???
      if ( true_chain.particle_flag == 2212 ) {
	if ( detect_flag_cut1 && detect_flag_cut2 ) {
	  proton_mom_cut2[imom] += ev.weight;
	  proton_ang_cut2[iang] += ev.weight;
	  proton_mom_ang_cut2[imom][iang] += ev.weight;
	  if ( detect_flag_cut3 ) {
	    proton_mom_cut3[imom] += ev.weight;
	    proton_ang_cut3[iang] += ev.weight;
	    proton_mom_ang_cut3[imom][iang] += ev.weight;
	  }
	}
      }
      else if ( std::abs(true_chain.particle_flag) == 211 ) {
	if ( detect_flag_cut1 && detect_flag_cut2 ) {
	  pion_mom_cut2[imom] += ev.weight;
	  pion_ang_cut2[iang] += ev.weight;
	  pion_mom_ang_cut2[imom][iang] += ev.weight;
	  if ( detect_flag_cut3 ) {
	    pion_mom_cut3[imom] += ev.weight;
	    pion_ang_cut3[iang] += ev.weight;
	    pion_mom_ang_cut3[imom][iang] += ev.weight;
	  }
	}
      }
    }

  }

  ofile->cd();
  hist_true_proton_mom->Write();
  hist_true_proton_mom_cut1->Write();
  hist_true_proton_mom_cut2->Write();
  hist_true_proton_mom_cut3->Write();
  hist_true_proton_ang->Write();
  hist_true_proton_ang_cut1->Write();
  hist_true_proton_ang_cut2->Write();
  hist_true_proton_ang_cut3->Write();
  hist_true_proton_mom_ang->Write();
  hist_true_proton_mom_ang_cut1->Write();
  hist_true_proton_mom_ang_cut2->Write();
  hist_true_proton_mom_ang_cut3->Write();
  hist_true_pion_mom->Write();
  hist_true_pion_mom_cut1->Write();
  hist_true_pion_mom_cut2->Write();
  hist_true_pion_mom_cut3->Write();
  hist_true_pion_ang->Write();
  hist_true_pion_ang_cut1->Write();
  hist_true_pion_ang_cut2->Write();
  hist_true_pion_ang_cut3->Write();
  hist_true_pion_mom_ang->Write();
  hist_true_pion_mom_ang_cut1->Write();
  hist_true_pion_mom_ang_cut2->Write();
  hist_true_pion_mom_ang_cut3->Write();

  otree->Fill();
  otree->Write();
  ofile->Close();

}
