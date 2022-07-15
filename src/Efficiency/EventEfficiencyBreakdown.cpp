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
     //logging::trivial::severity >= logging::trivial::debug
     );


  if ( argc != 5 ) {
    BOOST_LOG_TRIVIAL(error) << "Usage : " << argv[0]
			     << " <B2 file> <momch file w/o BM > <momch file w/ BM> <output file>";
    std::exit(1);
  }

  BOOST_LOG_TRIVIAL(info) << "==========Efficiency breakdown plots creation start==========";

  // input B2 file
  std::string inputfilename = argv[1];
  B2Reader reader(inputfilename);

  // input momch file
  std::string momchfilename = argv[2];
  auto ev_vec = Momentum_recon::ReadEventInformationBin(momchfilename);

  // input momch file w/ BM
  std::string momchbmfilename = argv[3];
  auto ev_bm_vec = Momentum_recon::ReadEventInformationBin(momchbmfilename);
  std::map<int, Momentum_recon::Event_information > ev_bm_map;
  for ( auto ev_bm : ev_bm_vec ) {
    ev_bm_map.insert(std::make_pair(ev_bm.groupid, ev_bm));
  }

  // output file
  std::string outputfilename = argv[4];
  TFile *ofile = new TFile((TString)outputfilename, "recreate");
  TTree *otree = new TTree("tree", "tree");

  TH1D *hist_true_muon_mom = new TH1D("hist_true_muon_mom", "", 100, 0., 2000.);
  TH1D *hist_true_muon_mom_cut1 = new TH1D("hist_true_muon_mom_cut1", "", 100, 0., 2000.); // ISS acceptance
  TH1D *hist_true_muon_mom_cut2 = new TH1D("hist_true_muon_mom_cut2", "", 100, 0., 2000.); // BM detection
  TH1D *hist_true_muon_mom_cut3 = new TH1D("hist_true_muon_mom_cut3", "", 100, 0., 2000.); // track matching
  TH1D *hist_true_muon_mom_cut4 = new TH1D("hist_true_muon_mom_cut4", "", 100, 0., 2000.); // FV
  TH1D *hist_true_muon_mom_cut5 = new TH1D("hist_true_muon_mom_cut5", "", 100, 0., 2000.); // good vertexing
  TH1D *hist_true_muon_mom_cut6 = new TH1D("hist_true_muon_mom_cut6", "", 100, 0., 2000.); // range measureble

  TH1D *hist_true_muon_ang = new TH1D("hist_true_muon_ang", "", 30, 0., 90.);
  TH1D *hist_true_muon_ang_cut1 = new TH1D("hist_true_muon_ang_cut1", "", 30, 0., 90.);
  TH1D *hist_true_muon_ang_cut2 = new TH1D("hist_true_muon_ang_cut2", "", 30, 0., 90.);
  TH1D *hist_true_muon_ang_cut3 = new TH1D("hist_true_muon_ang_cut3", "", 30, 0., 90.);
  TH1D *hist_true_muon_ang_cut4 = new TH1D("hist_true_muon_ang_cut4", "", 30, 0., 90.);
  TH1D *hist_true_muon_ang_cut5 = new TH1D("hist_true_muon_ang_cut5", "", 30, 0., 90.);
  TH1D *hist_true_muon_ang_cut6 = new TH1D("hist_true_muon_ang_cut6", "", 30, 0., 90.);

  TH2D *hist_true_muon_mom_ang = new TH2D("hist_true_muon_mom_ang", "", 100, 0., 2000., 30, 0., 90.);
  TH2D *hist_true_muon_mom_ang_cut1 = new TH2D("hist_true_muon_mom_ang_cut1", "", 100, 0., 2000., 30, 0., 90.);
  TH2D *hist_true_muon_mom_ang_cut2 = new TH2D("hist_true_muon_mom_ang_cut2", "", 100, 0., 2000., 30, 0., 90.);
  TH2D *hist_true_muon_mom_ang_cut3 = new TH2D("hist_true_muon_mom_ang_cut3", "", 100, 0., 2000., 30, 0., 90.);
  TH2D *hist_true_muon_mom_ang_cut4 = new TH2D("hist_true_muon_mom_ang_cut4", "", 100, 0., 2000., 30, 0., 90.);
  TH2D *hist_true_muon_mom_ang_cut5 = new TH2D("hist_true_muon_mom_ang_cut5", "", 100, 0., 2000., 30, 0., 90.);
  TH2D *hist_true_muon_mom_ang_cut6 = new TH2D("hist_true_muon_mom_ang_cut6", "", 100, 0., 2000., 30, 0., 90.);

  double mom_bin_edge[13] = {0., 100., 200., 300., 400., 500.,
			     600., 700., 800., 900., 1000., 1500., 2000.};
  double ang_bin_edge[11] = {0., 5., 10., 15., 20., 25., 30., 35., 40., 45., 90.};

  double mom_entry_total[12] = {};
  double ang_entry_total[10] = {};
  double mom_ang_entry_total[12][10] = {{}};
  double mom_cut0[12] = {};
  double mom_cut1[12] = {};
  double mom_cut2[12] = {};
  double mom_cut3[12] = {};
  double mom_cut4[12] = {};
  double mom_cut5[12] = {};
  double mom_cut6[12] = {};
  double mom_cut7[12] = {};
  double ang_cut0[10] = {};
  double ang_cut1[10] = {};
  double ang_cut2[10] = {};
  double ang_cut3[10] = {};
  double ang_cut4[10] = {};
  double ang_cut5[10] = {};
  double ang_cut6[10] = {};
  double ang_cut7[10] = {};
  double mom_ang_cut0[12][10] = {{}};
  double mom_ang_cut1[12][10] = {{}};
  double mom_ang_cut2[12][10] = {{}};
  double mom_ang_cut3[12][10] = {{}};
  double mom_ang_cut4[12][10] = {{}};
  double mom_ang_cut5[12][10] = {{}};
  double mom_ang_cut6[12][10] = {{}};
  double mom_ang_cut7[12][10] = {{}};

  otree->Branch("mom_entry_total", mom_entry_total, "mom_entry_total[12]/D");
  otree->Branch("ang_entry_total", ang_entry_total, "ang_entry_total[10]/D");
  otree->Branch("mom_ang_entry_total", mom_ang_entry_total, "mom_ang_entry_total[12][10]/D");
  otree->Branch("mom_cut0", mom_cut0, "mom_cut0[12]/D");
  otree->Branch("mom_cut1", mom_cut1, "mom_cut1[12]/D");
  otree->Branch("mom_cut2", mom_cut2, "mom_cut2[12]/D");
  otree->Branch("mom_cut3", mom_cut3, "mom_cut3[12]/D");
  otree->Branch("mom_cut4", mom_cut4, "mom_cut4[12]/D");
  otree->Branch("mom_cut5", mom_cut5, "mom_cut5[12]/D");
  otree->Branch("mom_cut6", mom_cut6, "mom_cut6[12]/D");
  otree->Branch("mom_cut7", mom_cut7, "mom_cut7[12]/D");
  otree->Branch("ang_cut0", ang_cut0, "ang_cut0[10]/D");
  otree->Branch("ang_cut1", ang_cut1, "ang_cut1[10]/D");
  otree->Branch("ang_cut2", ang_cut2, "ang_cut2[10]/D");
  otree->Branch("ang_cut3", ang_cut3, "ang_cut3[10]/D");
  otree->Branch("ang_cut4", ang_cut4, "ang_cut4[10]/D");
  otree->Branch("ang_cut5", ang_cut5, "ang_cut5[10]/D");
  otree->Branch("ang_cut6", ang_cut6, "ang_cut6[10]/D");
  otree->Branch("ang_cut7", ang_cut7, "ang_cut7[10]/D");
  otree->Branch("mom_ang_cut0", mom_ang_cut0, "mom_ang_cut0[12][10]/D");
  otree->Branch("mom_ang_cut1", mom_ang_cut1, "mom_ang_cut1[12][10]/D");
  otree->Branch("mom_ang_cut2", mom_ang_cut2, "mom_ang_cut2[12][10]/D");
  otree->Branch("mom_ang_cut3", mom_ang_cut3, "mom_ang_cut3[12][10]/D");
  otree->Branch("mom_ang_cut4", mom_ang_cut4, "mom_ang_cut4[12][10]/D");
  otree->Branch("mom_ang_cut5", mom_ang_cut5, "mom_ang_cut5[12][10]/D");
  otree->Branch("mom_ang_cut6", mom_ang_cut6, "mom_ang_cut6[12][10]/D");
  otree->Branch("mom_ang_cut7", mom_ang_cut7, "mom_ang_cut7[12][10]/D");

  for ( auto ev : ev_vec ) {

    Momentum_recon::Event_information ev_bm;
    if ( ev_bm_map.find(ev.groupid) != ev_bm_map.end() ) {
      ev_bm = ev_bm_map.at(ev.groupid);
    }
    else 
      ev_bm.weight = 0.;

    BOOST_LOG_TRIVIAL(debug) << "Entry : " << ev.groupid << ", " << ev_bm.groupid;
    BOOST_LOG_TRIVIAL(debug) << "weight : " << ev.weight << ", " << ev_bm.weight;
    reader.ReadSpill(ev.groupid);
    auto &spill_summary = reader.GetSpillSummary();

    bool mu_detect_flag_cut1 = false;
    bool mu_detect_flag_cut2 = false;
    bool mu_detect_flag_cut3 = false;
    bool mu_detect_flag_cut4 = false;
    bool mu_detect_flag_cut5 = false;
    bool mu_detect_flag_cut6 = false;

    double true_mu_mom = -1.;
    double true_mu_ang = -1.;
    int true_mu_id = -1;
    
    int imom = -1;
    int iang = -1;

    if ( ev.true_chains.empty() ) continue;

    for ( auto chain : ev.true_chains ) {
      if ( chain.particle_flag == 13 &&
	   chain.direction == 1 ) {
	true_mu_mom = chain.bm_range_mom;
	imom = true_mu_mom / 100;
	if ( imom >= 10 && imom < 15 ) imom = 10;
	else if ( imom >= 15 ) imom = 15;
	double ax = chain.base.back().ax;
	double ay = chain.base.back().ay;
	true_mu_ang = std::atan(std::hypot(ax, ay)) * TMath::RadToDeg();
	iang = true_mu_ang / 5.;
	if ( iang >= 9 ) iang = 9;
	true_mu_id = chain.chainid;
	hist_true_muon_mom->Fill(true_mu_mom, ev.weight);
	hist_true_muon_ang->Fill(true_mu_ang, ev.weight);
	break;
      }
    }

    if ( imom < 0 || iang < 0 ) continue;

    // ECC acceptance and reconstruction
    for ( auto chain : ev.chains ) {
      if ( (chain.chainid == true_mu_id) &&
	   (chain.base.front().pl == 3 || chain.base.front().pl == 4) ) {
	mu_detect_flag_cut1 = true;
	break;
      }
    }

    // Baby MIND detection
    double track_length = -1.;
    std::vector<B2TrackSummary* > tracks;
    if ( mu_detect_flag_cut1 ) {
      auto it_recon_vertex = spill_summary.BeginReconVertex();
      while ( auto *vertex = it_recon_vertex.Next() ) {
	auto it_track = vertex->BeginTrack();
	while ( auto *track = it_track.Next() ) {
	  if ( !track->HasDetector(B2Detector::kBabyMind)) continue;
	  if ( track->GetTrackType() != B2TrackType::kPrimaryTrack ) continue;
	  if ( track->GetPrimaryTrackType() == B2PrimaryTrackType::kBabyMind3DTrack ) {
	    if ( track->GetTrackLengthTotal() > track_length ) {
	      track_length = track->GetTrackLengthTotal();
	      tracks.push_back(track);
	    }
	  }
	  else if ( track->GetPrimaryTrackType() == B2PrimaryTrackType::kMatchingTrack ) {
	    if ( track->HasDetector(B2Detector::kProtonModule) ||
		 track->HasDetector(B2Detector::kWagasciUpstream) ||
		 track->HasDetector(B2Detector::kWagasciDownstream)) {
	      if ( track->GetTrackLengthTotal() > track_length ) {
		track_length = track->GetTrackLengthTotal();
		tracks.push_back(track);
	      }
	    }
	  }
	}
      }
    }

    if ( !tracks.empty() ) {
      if ( tracks.back()->GetParticlePdg() == PDG_t::kMuonMinus ) {
	int n_num_hit = 0;
	auto it_hit = tracks.back()->BeginHit();
	while ( const auto *hit = it_hit.Next() ) {
	  if ( hit->GetParentTrack().GetParticlePdg() == PDG_t::kMuonMinus )
	    n_num_hit++;
	}
	if ( n_num_hit > 2 )
	  mu_detect_flag_cut2 = true;
      }
    }

    // Muon track matching
    int muon_stop_flag = -1;
    if ( mu_detect_flag_cut2 ) {
      for ( auto chain : ev_bm.chains ) {
	int recon_particle_id = chain.particle_flag % 10000;
	int true_particle_id = chain.particle_flag / 10000;
	if ( chain.chainid == true_mu_id &&
	     recon_particle_id == 13 ) {
	  muon_stop_flag = chain.stop_flag;
	  mu_detect_flag_cut3 = true;
	  break;
	}
      }
    }

    // Fiducial volume
    if ( mu_detect_flag_cut3 &&
	 ev_bm.vertex_material >= 0 ) {
      mu_detect_flag_cut4 = true;
    }

    // Vertex plate is consistent with the true information
    if ( mu_detect_flag_cut4 && 
	 ev_bm.vertex_pl / 1000 == ev_bm.vertex_pl % 1000 ) {
      mu_detect_flag_cut5 = true;
    }
    
    // Baby MIND range measurable
    if ( mu_detect_flag_cut5 &&
	 muon_stop_flag == 1 ) {
      mu_detect_flag_cut6 = true;
    }

    mom_entry_total[imom] += 1.;
    ang_entry_total[iang] += 1.;
    mom_ang_entry_total[imom][iang] += 1.;
    mom_cut0[imom] += ev.weight;
    ang_cut0[iang] += ev.weight;
    mom_ang_cut0[imom][iang] += ev.weight;

    if ( mu_detect_flag_cut1 ) {
      hist_true_muon_mom_cut1->Fill(true_mu_mom, ev.weight);
      hist_true_muon_ang_cut1->Fill(true_mu_ang, ev.weight);
      hist_true_muon_mom_ang_cut1->Fill(true_mu_mom, true_mu_ang, ev.weight);
      mom_cut1[imom] += ev.weight;
      ang_cut1[iang] += ev.weight;
      mom_ang_cut1[imom][iang] += ev.weight;
    }
    if ( mu_detect_flag_cut2 ) {
      hist_true_muon_mom_cut2->Fill(true_mu_mom, ev.weight);
      hist_true_muon_ang_cut2->Fill(true_mu_ang, ev.weight);
      hist_true_muon_mom_ang_cut2->Fill(true_mu_mom, true_mu_ang, ev.weight);
      mom_cut2[imom] += ev.weight;
      ang_cut2[iang] += ev.weight;
      mom_ang_cut2[imom][iang] += ev.weight;
    }
    if ( mu_detect_flag_cut3 ) {
      hist_true_muon_mom_cut3->Fill(true_mu_mom, ev_bm.weight);
      hist_true_muon_ang_cut3->Fill(true_mu_ang, ev_bm.weight);
      hist_true_muon_mom_ang_cut3->Fill(true_mu_mom, true_mu_ang, ev_bm.weight);
      mom_cut3[imom] += ev_bm.weight;
      ang_cut3[iang] += ev_bm.weight;
      mom_ang_cut3[imom][iang] += ev_bm.weight;
    }
    if ( mu_detect_flag_cut4 ) {
      hist_true_muon_mom_cut4->Fill(true_mu_mom, ev_bm.weight);
      hist_true_muon_ang_cut4->Fill(true_mu_ang, ev_bm.weight);
      hist_true_muon_mom_ang_cut4->Fill(true_mu_mom, true_mu_ang, ev_bm.weight);
      mom_cut4[imom] += ev_bm.weight;
      ang_cut4[iang] += ev_bm.weight;
      mom_ang_cut4[imom][iang] += ev_bm.weight;
    }
    if ( mu_detect_flag_cut5 ) {
      hist_true_muon_mom_cut5->Fill(true_mu_mom, ev_bm.weight);
      hist_true_muon_ang_cut5->Fill(true_mu_ang, ev_bm.weight);
      hist_true_muon_mom_ang_cut5->Fill(true_mu_mom, true_mu_ang, ev_bm.weight);
      mom_cut5[imom] += ev_bm.weight;
      ang_cut5[iang] += ev_bm.weight;
      mom_ang_cut5[imom][iang] += ev_bm.weight;
    }
    if ( mu_detect_flag_cut6 ) {
      hist_true_muon_mom_cut6->Fill(true_mu_mom, ev_bm.weight);
      hist_true_muon_ang_cut6->Fill(true_mu_ang, ev_bm.weight);
      hist_true_muon_mom_ang_cut6->Fill(true_mu_mom, true_mu_ang, ev_bm.weight);
      mom_cut6[imom] += ev_bm.weight;
      ang_cut6[iang] += ev_bm.weight;
      mom_ang_cut6[imom][iang] += ev_bm.weight;
    }

  }

  ofile->cd();
  hist_true_muon_mom->Write();
  hist_true_muon_mom_cut1->Write();
  hist_true_muon_mom_cut2->Write();
  hist_true_muon_mom_cut3->Write();
  hist_true_muon_mom_cut4->Write();
  hist_true_muon_mom_cut5->Write();
  hist_true_muon_mom_cut6->Write();
  hist_true_muon_ang->Write();
  hist_true_muon_ang_cut1->Write();
  hist_true_muon_ang_cut2->Write();
  hist_true_muon_ang_cut3->Write();
  hist_true_muon_ang_cut4->Write();
  hist_true_muon_ang_cut5->Write();
  hist_true_muon_ang_cut6->Write();
  hist_true_muon_mom_ang->Write();
  hist_true_muon_mom_ang_cut1->Write();
  hist_true_muon_mom_ang_cut2->Write();
  hist_true_muon_mom_ang_cut3->Write();
  hist_true_muon_mom_ang_cut4->Write();
  hist_true_muon_mom_ang_cut5->Write();
  hist_true_muon_mom_ang_cut6->Write();

  otree->Fill();
  otree->Write();
  ofile->Close();

}
