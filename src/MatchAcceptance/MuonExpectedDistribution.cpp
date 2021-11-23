// boost includes
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

// B2 includes
#include <B2Reader.hh>
#include <B2SpillSummary.hh>
#include <B2TrackSummary.hh>
#include <B2EventSummary.hh>
#include <B2VertexSummary.hh>
#include <B2EmulsionSummary.hh>
#include <B2Pdg.hh>

// NTBM includes
#include <NTBMSummary.hh>
#include <NTBMConst.hh>

// root includes
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TVector3.h>

namespace logging = boost::log;
namespace fs = boost::filesystem;

Int_t GetInteractionEcc(const B2VertexSummary &primary_vertex_summary) {
  if ( primary_vertex_summary.GetDetector() != B2Detector::kNinja ) 
    return -1;
  TVector3 vertex_position = primary_vertex_summary.GetAbsolutePosition().GetValue();
  vertex_position.SetX(vertex_position.X() - NINJA_POS_X - NINJA_ECC_POS_X);
  vertex_position.SetY(vertex_position.Y() - NINJA_POS_Y - NINJA_ECC_POS_Y);
  Int_t top_id, side_id;

  if ( std::fabs(vertex_position.X() + NINJA_ECC_GAP_X) < 0.5 * NINJA_DESIC_WIDTH )
    top_id = 0;
  else if ( std::fabs(vertex_position.X()) < 0.5 * NINJA_DESIC_WIDTH )
    top_id = 1;
  else if ( std::fabs(vertex_position.X() - NINJA_ECC_GAP_X) < 0.5 * NINJA_DESIC_WIDTH )
    top_id = 2;
  else return -1;

  if ( std::fabs(vertex_position.Y() + NINJA_ECC_GAP_Y) < 0.5 * NINJA_DESIC_HEIGHT )
    side_id = 0;
  else if ( std::fabs(vertex_position.Y()) < 0.5 * NINJA_DESIC_HEIGHT )
    side_id = 1;
  else if ( std::fabs(vertex_position.Y() - NINJA_ECC_GAP_Y) < 0.5 * NINJA_DESIC_HEIGHT )
    side_id = 2;
  else return -1;

  return 3 * side_id + top_id;

}

int main (int argc, char *argv[]) {

  logging::core::get()->set_filter
    (
     logging::trivial::severity >= logging::trivial::info
     );

    if ( argc != 4 ) {
      BOOST_LOG_TRIVIAL(error) << "Usage : " << argv[0]
			       << " <input B2 file name> <input NTBM file name> <output file name>";
      std::exit(1);
    }

  BOOST_LOG_TRIVIAL(info) << "==========Start==========";

  const fs::path input_file(argv[1]);

  try {

    B2Reader reader(input_file);

    TFile *nt_file = new TFile(argv[2], "read");
    TTree *nt_tree = (TTree*)nt_file->Get("tree");
    NTBMSummary *ntbm = nullptr;
    nt_tree->SetBranchAddress("NTBMSummary", &ntbm);
    Int_t nt_entry = 0;

    TFile *ofile = new TFile(argv[3], "recreate");
    TH1D *hist_muon_mom = new TH1D("hist_muon_mom", "Muon momentum;p_{#mu} [GeV/c];Entries/0.1 GeV/c", 40, 0, 4);
    TH1D *hist_muon_cos = new TH1D("hist_muon_cos", "Muon angle;cos#theta_{#mu};Entries/0.01", 200, -1, 1); // cos theta distribution
    TH1D *hist_muon_deg = new TH1D("hist_muon_deg", "Muon angle;#theta_{#mu} [degree];Entries/2 deg", 90, 0, 180);// angle [degree] distribution
    TH2D *hist_muon_mom_cos = new TH2D("hist_muon_mom_cos", "Muon 2d kinematics;p_{#mu} [GeV/c];cos#theta_{#mu}", 40, 0, 4, 200, -1, 1);
    TH2D *hist_muon_mom_deg = new TH2D("hist_muon_mom_deg", "Muon 2d kinematics;p_{#mu} [GeV/c];#theta_{#mu} [degree]", 40, 0, 4, 90, 0, 180);
    TH1D *hist_muon_bm_deg = new TH1D("hist_muon_bm_deg", ";Baby MIND reconstructed angle (#theta_{#mu, BM}) [degree];Entries/2 deg", 90, 0, 180);
    TH2D *hist_muon_deg_bm_deg = new TH2D("hist_muon_deg_bm_deg",";Baby MIND reconstructed angle (#theta_{#mu, BM}) [degree];#theta_{#mu} [degree]",
					  90, 180, 0, 90, 180, 0);

    while ( reader.ReadNextSpill() > 0 ) {
      nt_tree->GetEntry(nt_entry);

      auto &input_spill_summary = reader.GetSpillSummary();

      auto it_event = input_spill_summary.BeginTrueEvent();
      const auto *event = it_event.Next(); // get true event summary
      const Double_t norm = event->GetNormalization();

      auto &primary_vertex_summary = event->GetPrimaryVertex();

      Int_t interaction_ecc_id = GetInteractionEcc(primary_vertex_summary);
      const Double_t total_cross_section = primary_vertex_summary.GetTotalCrossSection();
      const Double_t weight = norm * total_cross_section * 1.e-38 * 6.02e23;

      Bool_t detect_muon = false;
      Double_t muon_mom = -1; 
      Double_t muon_ang = -1;
      Double_t muon_cos = -2;

      // Check if a muon is generated from the intreaction
      auto it_true_track = primary_vertex_summary.BeginTrack();
      while ( const auto *track = it_true_track.Next() ) {
	const Int_t particle_id = track->GetParticlePdg();
	if (particle_id == PDG_t::kMuonMinus) {
	  BOOST_LOG_TRIVIAL(debug) << "Found muon!";
	  detect_muon = true;
	  muon_mom = track->GetInitialAbsoluteMomentum().GetValue() / 1.e3;
	  TVector3 direction = track->GetInitialDirection().GetValue();
	  muon_ang = TMath::ATan(TMath::Hypot(direction.X()/direction.Z(), direction.Y()/direction.Z())) * TMath::RadToDeg();
	  muon_cos = TMath::Cos(muon_ang * TMath::DegToRad());
	  break;
	}
      }

      // Check if the muon is w/i acceptance of the ISS/OSS
      Bool_t detect_muon_iss = false;
      Bool_t detect_muon_oss = false;
      Int_t iss_ecc_id = -1;
      Int_t oss_ecc_id = -1;
      TVector3 tss_position;
      TVector3 tss_tangent;
      auto it_emulsion = input_spill_summary.BeginEmulsion();
      while ( const auto *emulsion = it_emulsion.Next() ) {
	if ( emulsion->GetParentTrackId() == 0 ||
	     emulsion->GetParentTrackId() >= primary_vertex_summary.GetNumOutgoingTracks() ) continue;
	if ( emulsion->GetParentTrack().GetParticlePdg() != PDG_t::kMuonMinus ) 
	  continue;
	if ( (emulsion->GetFilmType() == B2EmulsionType::kECC) &&
	     (emulsion->GetPlate() == 2 || emulsion->GetPlate() == 3) ) {
	  detect_muon_iss = true;
	  iss_ecc_id = emulsion->GetEcc();
	} // track detected in Inside SS
	else if ( (emulsion->GetFilmType() == B2EmulsionType::kShifter) &&
		  (emulsion->GetPlate() < 4 )) {
	  detect_muon_oss = true;
	  oss_ecc_id = emulsion->GetEcc();
	} // track detected in Outside SS
	else if ( (emulsion->GetFilmType() == B2EmulsionType::kShifter) &&
		  (emulsion->GetPlate() == 15)) {
	  tss_position = emulsion->GetAbsolutePosition().GetValue();
	  tss_tangent = emulsion->GetAbsolutePosition().GetValue();
	}
	else continue;
      }

      if ( !detect_muon_iss || !detect_muon_oss )
	detect_muon = false;
      if ( interaction_ecc_id != iss_ecc_id ||
	   interaction_ecc_id != oss_ecc_id )
	detect_muon = false; // interaction ECC and ISS/OSS ECC should be the same

      // Check if the track penetrates the tracker fiducial area
      TVector3 extrapolated_position;
      extrapolated_position.SetX(tss_position.X() + tss_tangent.X() * 30.);
      extrapolated_position.SetY(tss_position.Y() + tss_tangent.Y() * 10.);
      extrapolated_position.SetX(extrapolated_position.X() - NINJA_POS_X - NINJA_TRACKER_POS_X);
      extrapolated_position.SetY(extrapolated_position.Y() - NINJA_POS_Y - NINJA_TRACKER_POS_Y);

      if ( extrapolated_position.Y() < -448. ||
	   extrapolated_position.Y() > 600. ||
	   extrapolated_position.X() < -600. ||
	   extrapolated_position.X() > 448. ) continue;

      // Check if Baby MIND has track made by the muon
      if (detect_muon) {
	Bool_t detect_muon_bm = false;
	std::vector<const B2TrackSummary*> recon_track_vector;
	recon_track_vector.reserve(input_spill_summary.GetNumReconTracks());
	auto it_recon_track = input_spill_summary.BeginReconTrack();
	while ( const auto *recon_track = it_recon_track.Next() ) {
	  if ( recon_track->GetTrackType() == B2TrackType::kPrimaryTrack ) {
	    if ( recon_track->GetPrimaryTrackType() == B2PrimaryTrackType::kBabyMind3DTrack ) {
	      recon_track_vector.push_back(recon_track);
	    }
	    else if ( recon_track->GetPrimaryTrackType() == B2PrimaryTrackType::kMatchingTrack &&
		      recon_track->HasDetector(B2Detector::kBabyMind) ) {
              if ( recon_track->HasDetector(B2Detector::kProtonModule) ||
		   recon_track->HasDetector(B2Detector::kWagasciUpstream) ||
		   recon_track->HasDetector(B2Detector::kWagasciDownstream) ) {
		recon_track_vector.push_back(recon_track);
              }
	    }
	  }
	}

	Int_t matched_babymind_track_id = -1;
	for ( int i = 0; i < recon_track_vector.size(); i++ ) {
	  if ( recon_track_vector.at(i)->GetParticlePdg() == PDG_t::kMuonMinus ) {
	    detect_muon_bm = true;
	    matched_babymind_track_id = i;
	  }	    
	}

	// Search corresponding muon track from NTBMSummary and get 3d pre-reconstructed angle
	if ( detect_muon_bm ) {
	  std::vector<Double_t> tangent = ntbm->GetBabyMindTangent(matched_babymind_track_id);
	  Double_t theta_bm = TMath::ATan(TMath::Hypot(tangent.at(0), tangent.at(1))) * TMath::RadToDeg();

	  hist_muon_mom->Fill(muon_mom, weight);
	  hist_muon_cos->Fill(muon_cos, weight);
	  hist_muon_deg->Fill(muon_ang, weight);
	  hist_muon_mom_cos->Fill(muon_mom, muon_cos, weight);
	  hist_muon_mom_deg->Fill(muon_mom, muon_ang, weight);
	  hist_muon_bm_deg->Fill(theta_bm, weight);
	  hist_muon_deg_bm_deg->Fill(muon_ang, theta_bm, weight);
	}
      }
  
      nt_entry++;
    }

    ofile->cd();
    hist_muon_mom->Write();
    hist_muon_cos->Write();
    hist_muon_deg->Write();
    hist_muon_mom_cos->Write();
    hist_muon_mom_deg->Write();
    hist_muon_bm_deg->Write();
    hist_muon_deg_bm_deg->Write();
    ofile->Close();

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
