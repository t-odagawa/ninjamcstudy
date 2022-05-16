// boost includes
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

// B2 includes
#include <B2Reader.hh>
#include <B2SpillSummary.hh>
#include <B2TrackSummary.hh>
#include <B2EventSummary.hh>
#include <B2ClusterSummary.hh>
#include <B2HitSummary.hh>
#include <B2VertexSummary.hh>
#include <B2EmulsionSummary.hh>
#include <B2Enum.hh>
#include <B2Pdg.hh>

// root includes
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TMath.h>
#include <TVector3.h>

#include "MyFunction.hpp"

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

int main ( int argc, char *argv[] ) {
  
  logging::core::get()->set_filter
    (
     logging::trivial::severity >= logging::trivial::debug
     );

  if ( argc != 3 ) {
    BOOST_LOG_TRIVIAL(error) << "Usage : " << argv[0]
			     << " <input B2 filename> <output file name>";
    std::exit(1);
  }

  const fs::path input_file(argv[1]);

  try {

    B2Reader reader(input_file);

    TFile *ofile = new TFile(argv[2], "recreate");
    TH1D *hist_muon_tan_x = new TH1D("hist_muon_tan_x", ";Shifter angle (tan#theta_{#mu,x});Entries/0.04", 30, 0, 1.2);
    TH1D *hist_muon_tan_y = new TH1D("hist_muon_tan_y", ";Shifter angle (tan#theta_{#mu,y});Entries/0.04", 30, 0, 1.2);
    TH1D *hist_ang_3d = new TH1D("hist_ang_3d", "", 100, 0, 50);

    while ( reader.ReadNextSpill() > 0 ) {
      auto &input_spill_summary = reader.GetSpillSummary();

      auto it_event = input_spill_summary.BeginTrueEvent();
      const auto *event = it_event.Next();
      const Double_t norm = event->GetNormalization();

      auto &primary_vertex_summary = event->GetPrimaryVertex();

      Int_t interaction_ecc_id = GetInteractionEcc(primary_vertex_summary);
      const Double_t total_cross_section = primary_vertex_summary.GetTotalCrossSection();
      const Double_t weight = norm * total_cross_section * 1.e-38 * 6.02e23;     

      std::vector<Int_t > bm_hit_id;
      auto it_true_track = primary_vertex_summary.BeginTrack();
      while ( const auto *track = it_true_track.Next() ) {

	if ( !track->HasDetector(B2Detector::kBabyMind) ) continue;

	if ( B2Pdg::IsMuonPlusOrMinus(track->GetParticlePdg()) ) {
	  BOOST_LOG_TRIVIAL(debug) << "Found muon";
	  auto it_cluster = track->BeginCluster();
	  const auto *cluster = it_cluster.Next();
	  auto it_hit = cluster->BeginHit();
	  while ( const auto *hit = it_hit.Next() ) {
	    if ( hit->GetDetectorId() == B2Detector::kBabyMind ) {
	      bm_hit_id.push_back(hit->GetHitId());
	    }
	  }
	  break;
	}
      }

      if ( bm_hit_id.empty() ) continue;

      // Check if the muon is w/i acceptance of the ISS/OSS
      Bool_t detect_muon_iss = false;
      Bool_t detect_muon_oss = false;
      Int_t iss_ecc_id = -1;
      Int_t oss_ecc_id = -1;
      TVector3 tss_position;
      TVector3 tss_tangent;
      Double_t muon_tan_x = -1;
      Double_t muon_tan_y = -1;
      Double_t ang_3d = -1;
      auto it_emulsion = input_spill_summary.BeginEmulsion();
      while ( const auto *emulsion = it_emulsion.Next() ) {
	if ( emulsion->GetParentTrackId() == 0 ||
	     emulsion->GetParentTrackId() >= primary_vertex_summary.GetNumOutgoingTracks() )
	  continue;
	if ( (emulsion->GetFilmType() == B2EmulsionType::kECC) &&
	     (emulsion->GetPlate() == 2 || emulsion->GetPlate() == 3) ) {
	  detect_muon_iss = true;
	  iss_ecc_id = emulsion->GetEcc();
	} // track detected in ISS
	else if ( (emulsion->GetFilmType() == B2EmulsionType::kShifter) &&
		  (emulsion->GetPlate() < 4) ) {
	  detect_muon_oss = true;
	  oss_ecc_id = emulsion->GetEcc();
	} // track detected in OSS
	else if ( (emulsion->GetFilmType() == B2EmulsionType::kShifter) &&
		  (emulsion->GetPlate() == 15) ) {
	  tss_position = emulsion->GetAbsolutePosition().GetValue();
	  tss_tangent = emulsion->GetTangent().GetValue();
	  muon_tan_x = std::fabs(tss_tangent.X());
	  muon_tan_y = std::fabs(tss_tangent.Y());
	  ang_3d = std::atan(std::hypot(tss_tangent.X(), tss_tangent.Y())) * TMath::RadToDeg();
	}
	else continue;
      }

      if ( !detect_muon_iss || !detect_muon_oss )
	continue;
      if ( interaction_ecc_id != iss_ecc_id ||
	   interaction_ecc_id != oss_ecc_id )
	continue;

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

      BOOST_LOG_TRIVIAL(debug) << "Target muon";

      // Check if Baby MIND has track made by the muon
      Bool_t detect_muon_bm = false;
      Int_t num_bm_hit_in_recon_track = 0;
      auto it_recon_vertex = input_spill_summary.BeginReconVertex();
      while ( auto *vertex = it_recon_vertex.Next() ) {
	auto it_outgoing_track = vertex->BeginTrack();
	while ( auto *track = it_outgoing_track.Next() ) {
	  auto it_hit = track->BeginHit();
	  while ( auto *hit = it_hit.Next() ) {
	    if ( hit->GetDetectorId() != B2Detector::kBabyMind ) continue;
	    if ( std::find(bm_hit_id.begin(),
			   bm_hit_id.end(),
			   hit->GetHitId()) != bm_hit_id.end() ) {
	      num_bm_hit_in_recon_track++;
	    }
	  }
	  if ( num_bm_hit_in_recon_track > 3 ) break;
	  else num_bm_hit_in_recon_track = 0;
	}
      }

      if ( num_bm_hit_in_recon_track < 4 ) continue;

      hist_muon_tan_x->Fill(muon_tan_x, weight);
      hist_muon_tan_y->Fill(muon_tan_y, weight);
      hist_ang_3d->Fill(ang_3d, weight);
    }

    ofile->cd();
    hist_muon_tan_x->Write();
    hist_muon_tan_y->Write();
    hist_ang_3d->Write();
    ofile->Close();

  } catch (const std::runtime_error &error) {
    BOOST_LOG_TRIVIAL(fatal) << "Runtime error : " << error.what();
    std::exit(1);
  } catch (const std::invalid_argument &error) {
    BOOST_LOG_TRIVIAL(fatal) << "Invalid argument error : " << error.what();
    std::exit(1);
  }

  std::exit(0);

}
