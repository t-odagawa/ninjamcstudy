// system includes
#include <vector>
#include <algorithm>

// boost includes
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

// B2 includes
#include <B2Reader.hh>
#include <B2Enum.hh>
#include <B2Const.hh>
#include <B2SpillSummary.hh>
#include <B2TrackSummary.hh>
#include <B2EventSummary.hh>
#include <B2EmulsionSummary.hh>

// NTBM includes
#include <NTBMConst.hh>
#include <NTBMSummary.hh>

// root includes
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TVector3.h>

namespace logging = boost::log;

const TString modename[5] = {"CCQE", "2p2h", "CC 1#pi", "CC Multi#pi", "CC Other", "NC"};
const Int_t color[5] = {424, 624, 394, 395, 401, 408};
const Int_t style[5] = {1001, 1001, 3006, 3005, 1001, 1001};

int main (int argc, char *argv[]) {

  logging::core::get()->set_filter
    (
     //logging::trivial::severity >= logging::trivial::debug
     logging::trivial::severity >= logging::trivial::info
     );

  if (argc != 3) {
    BOOST_LOG_TRIVIAL(error) << "Usage : " << argv[0]
			     << " <fileid> <output file name>";
    std::exit(1);
  }

  BOOST_LOG_TRIVIAL(info) << "==========Start==========";

  try {

    TString fileid_string = argv[1];
    TString matchfilename = "/home/t2k/odagawa/data/mc_data/pra1/trackmatch/ninja_mc_h2o_ninjamatch_";
    matchfilename += fileid_string;
    matchfilename += ".root";
    BOOST_LOG_TRIVIAL(debug) << "NTBM filename : " << matchfilename;
    TFile *matchfile = new TFile(matchfilename, "read");
    TTree *matchtree = (TTree*)matchfile->Get("tree");
    NTBMSummary *ntbm = nullptr;
    matchtree->SetBranchAddress("NTBMSummary", &ntbm);

    TString b2filename = "/home/t2k/odagawa/data/mc_data/pra1/ninja_mc_h2o_";
    b2filename += fileid_string;
    b2filename += ".root";
    BOOST_LOG_TRIVIAL(debug) << "B2MC filename : " << b2filename;
    B2Reader reader(b2filename);

    TFile *outputfile = new TFile(argv[2], "recreate");
    THStack *hs = new THStack("hs", "");
    TH1D *h_multi = new TH1D("h_multi", "Charged Track Multiplicty;# of particles;Entries", 10, 0, 10);
    TH1D *h_pro_mom = new TH1D("h_pro_mom", "Proton Momentum;Momentum [GeV/c];Entries", 8, 0, 1.6);
    TH1D *h_pro_ang = new TH1D("h_pro_ang", "Proton Angle;Angle [deg];Entries", 10, 0, 180);
    TH2D *h_pro_mom_ang = new TH2D("h_pro_mom_ang", ";Momentum [GeV/c];Angle [deg]", 32, 0, 1.6, 45, 0, 180);

    for (Int_t imatchentry = 0; imatchentry < matchtree->GetEntries(); imatchentry++) {
      matchtree->GetEntry(imatchentry);
      Double_t weight = ntbm->GetNormalization() * ntbm->GetTotalCrossSection();
      weight *= std::pow(10, -38) * 6.02 * std::pow(10, 23) * 1. * 13.34 * 0.292 / 960.;
      Bool_t analyze_event = false;

      if (ntbm->GetNumberOfTracks() < 1) continue; // No muon ID
      for (Int_t imatchcluster = 0; imatchcluster < ntbm->GetNumberOfNinjaClusters(); imatchcluster++) {
	Int_t babymind_track_id = ntbm->GetBabyMindTrackId(imatchcluster);
	if (babymind_track_id != -1) {
	  analyze_event = true;
	  break;
	}	
      } // imatchcluster
      
      if (analyze_event && reader.ReadSpill(imatchentry) > 0) {
	auto &spill_summary = reader.GetSpillSummary();
	const Int_t num_true_tracks = spill_summary.GetNumTrueTracks();
	BOOST_LOG_TRIVIAL(debug) << "Muon matched, entry : " << imatchentry
				 << " # of true tracks : " << num_true_tracks;
	Int_t emulsion_hits[num_true_tracks] = {};
        Double_t proton_momentum[num_true_tracks] = {};
	Double_t proton_tan_x[num_true_tracks] = {};
	Double_t proton_tan_y[num_true_tracks] = {};
	Double_t proton_angle[num_true_tracks] = {};

	// Analyze all emulsion hits
	BOOST_LOG_TRIVIAL(debug) << "Analyze all emulsion hits";
	auto it_emulsion = spill_summary.BeginEmulsion();
	while (const auto *emulsion = it_emulsion.Next()) {
	  BOOST_LOG_TRIVIAL(debug) << *emulsion;
	  if (emulsion->GetParentTrackId() < 1 ||
	      emulsion->GetParentTrackId() > num_true_tracks) continue;
	  auto &track_summary = emulsion->GetParentTrack();
	  BOOST_LOG_TRIVIAL(debug) << "True track id : " << emulsion->GetParentTrackId();
	  if (emulsion->GetTangent().GetValue().X() > 4.0 ||
	      emulsion->GetTangent().GetValue().Y() > 4.0 ||
	      emulsion->GetEcc() != 4 ||
	      emulsion->GetFilmType() != B2EmulsionType::kECC) continue;
	  emulsion_hits[emulsion->GetParentTrackId()]++;
	}

	// Then loop for all true tracks
	BOOST_LOG_TRIVIAL(debug) << "Analyze all true tracks";
	auto it_true_track = spill_summary.BeginTrueTrack();
	while (const auto *true_track = it_true_track.Next()) {
	  if (true_track->GetParticlePdg() != 2212) continue;
	  UInt_t true_track_id = true_track->GetTrackId();
	  proton_momentum[true_track_id] = std::sqrt(std::pow(true_track->GetTotalEnergy().GetValue(), 2.0)
						     - std::pow(true_track->GetMass().GetValue(), 2.0)) / 1000.;
	  proton_tan_x[true_track_id] = TMath::Tan(true_track->GetViewAngle().GetValue().X() * TMath::DegToRad());
	  proton_tan_y[true_track_id] = TMath::Tan(true_track->GetViewAngle().GetValue().Y() * TMath::DegToRad());
	  TVector3 initial_direction = true_track->GetInitialDirection().GetValue();
	  Double_t internal_product = initial_direction.Z();
	  Double_t distance = initial_direction.Mag();
	  proton_angle[true_track_id] = TMath::ACos(internal_product / distance) * TMath::RadToDeg();
	  BOOST_LOG_TRIVIAL(debug) << "Proton angle" << proton_angle[true_track_id] << " [deg]";
	}
	
	// Fill histograms
	BOOST_LOG_TRIVIAL(debug) << "Fill histograms";
	Int_t num_detected_tracks = 0;
	for (Int_t itrue_track = 0; itrue_track < num_true_tracks; itrue_track++) {
	  if (emulsion_hits[itrue_track] < 2) continue;
	  num_detected_tracks++;
	  if (proton_momentum[itrue_track] == 0. || proton_angle[itrue_track] == 0.) continue;
	  if (fabs(proton_tan_x[itrue_track]) > 4. || fabs(proton_tan_y[itrue_track]) > 4.) continue;
	  h_pro_mom->Fill(proton_momentum[itrue_track], weight);
	  h_pro_ang->Fill(proton_angle[itrue_track], weight);
	  h_pro_mom_ang->Fill(proton_momentum[itrue_track], proton_angle[itrue_track], weight);
	}

	if (num_detected_tracks > 0)
	  h_multi->Fill((Double_t)num_detected_tracks, weight);

      }

    } // imatchentry

    BOOST_LOG_TRIVIAL(info) << "Create output file : " << argv[2];
    outputfile->cd();
    h_multi->Write();
    h_pro_mom->Write();
    h_pro_ang->Write();
    h_pro_mom_ang->Write();
    outputfile->Close();

  } catch (const std::runtime_error &error) {
    BOOST_LOG_TRIVIAL(fatal) << "Runtime error : " << error.what();
    std::exit(1);
  } catch (const std::invalid_argument &error) {
    BOOST_LOG_TRIVIAL(fatal) << "Invalid argument error : " << error.what();
    std::exit(1);
  }

  BOOST_LOG_TRIVIAL(info) << "==========Finish==========";
  std::exit(0);

}
