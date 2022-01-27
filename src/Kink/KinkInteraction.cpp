// system includes
#include <fstream>

// boost includes
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

// B2 includes
#include <B2Reader.hh>
#include <B2SpillSummary.hh>
#include <B2ClusterSummary.hh>
#include <B2TrackSummary.hh>
#include <B2EventSummary.hh>
#include <B2VertexSummary.hh>
#include <B2HitSummary.hh>
#include <B2Enum.hh>
#include <B2Pdg.hh>

// MCs includes
#include <McsConst.hpp>
#include <McsClass.hpp>
#include <McsFunction.hpp>

// root includes
#include <TVector3.h>

// my includes
#include "MyFunction.hpp"

namespace logging = boost::log;
namespace fs = boost::filesytem;

int main (int argc, char *argv[]) {

  logging::core::get()->set_filter
    (
     logging::trivial::severity >= logging::trivial::debug
     );

  if ( argc != 3 ) {
    BOOST_LOG_TRIVIAL(error) << "Usage : " << argv[0]
			     << " <input B2 file name> <output file name>";
    std::exit(1);
  }

  BOOST_LOG_TRIVIAL(info) << "==========Start==========";

  const fs::path input_file(argv[1]);

  try {

    B2Reader reader(input_file);

    std::ofstream ofs(argv[2], std::ios::binary);
    Momentum_recon::Mom_chain mom_chain;
    Momentum_recon::Mom_basetrack mom_basetrack;
    std::pair<Momentum_recon::Mom_basetrack, Momentum_recon::Mom_basetrack> mom_basetrack_pair;

    while ( reader.ReadNextSpill() > 0 ) {
      
      auto &spill_summary = reader.GetSpillSummary();

      // Get muon track id
      int muon_track_id = -2;
      std::vector<int> ninja_hit_id;
      std::vector<int> bm_hit_id;
      auto it_event = spill_summary.BeginTrueEvent();
      const auto *event = it_event.Next();
      double norm = event->GetNormalization();
      auto &primary_vertex_summary = event->GetPrimaryVertex();
      int vertex_ecc = GetInteractionEcc(primary_vertex_summary);
      double total_cross_section = primary_vertex_summary.GetTotalCrossSection();
      auto it_outgoing_track = primary_vertex_summary.BeginTrack();

      while ( const auto *track = it_outgoing_track.Next() ) {
	
	if ( !track->HasDetector(B2Detector::kNinja) &&
	     !track->HasDetector(B2Detector::kBabyMind) ) continue;

	if ( B2Pdg::IsMuonPlusOrMinus(track->GetParticlePdg()) ) {
	  muon_track_id = track->GetTrackId();
	  // Get scintillator hit id of true muon track
	  auto it_cluster = track->BeginCluster();
	  const auto *cluster = it_cluster.Next();
	  auto it_hit = cluster->BeginHit();
	  while ( const auto *hit = it_hit.Next() ) {
	    if ( hit->GetDetectorId() == B2Detector::kBabyMind )
	      bm_hit_id.push_back(hit->GetHitId());
	    if ( hit->GetDetectorId() == B2Detector::kNinja )
	      ninja_hit_id.push_back(hit->GetHitId());
	  }
	  break;
	}
      }

      if ( muon_track_id < 0 ) continue;

      // Check muon id & collect basetracks
      Bool_t muon_detect_flag = false;
      Bool_t iss_detect_flag = false;
      Bool_t oss_detect_flag = false;      
      Bool_t bm_nt_connect_flag = false;

      std::vector<const B2EmulsionSummary*> emulsions;
      auto it_emulsion = spill_summary.BeginEmusion();
      while ( const auto *emulsion = it_emulsion.Next() ) {
	if ( emulsion->GetParentTrackId() == 0 ) continue;
	if ( emulsion->GetEcc() != vertex_ecc ) continue;
	if ( emulsion->GetFilmType() == B2EmulsionType::kECC ) {
	  
	}
	else if ( emulsion->GetFilmType() == B2EmulsionType::kShifter ) {
	  if ( emulsion->GetPlate() == 0 ) {
	    oss_detect_flag = true;
	  }
	}
      }
      
      // Check vertex ECC id

      std::vector<std::vector<const B2EmulsionSummary>> particle_emulsions;

      for ( auto one_particle_emulsions : particle_emulsions ) {
	int num_base = one_particle_emulsions.size();
	int num_link = num_base - 1;

	mom_chain.base.clear();
	mom_chain.base_pair.clear();
	mom_chain.base.reserve(num_base);
	mom_chain.base_pair.reserve(num_link);

	mom_chain.groupid = reader.GetEntryNumber();
	mom_chain.chainid = one_particle_emulsions.at(0)->GetParentTrackId();
	mom_chain.unixtime = -1;
	mom_chain.tracker_track_id = muon_track_id;
	mom_chain.entry_in_daily_file = reader.GetEntryNumber();
	mom_chain.stop_flag = 0;
	mom_chain.particle_flag = one_particle_emulsions.at(0)->GetParentTrack().GetParticlePdg();
	mom_chain.bm_range_mom = one_particle_emulsions.at(0)->GetParentTrakc().GetInitialAbsoluteMomentum().GetValue();
	mom_chain.bm_curvature_mom = norm * total_cross_section * 1e-38 * 6.02e23;

	Double_t downstream_position_z = 0;
	
	for ( auto emulsion : one_particle_emulsions ) {



	} // emulsion

	Momentum_recon::WriteMomChainHeader(ofs, mom_chain);
	for ( int ibase = 0; ibase < num_base; ibase++ ) {
	  ofs.write((char*)& mom_chain.base.at(ibase), sizeof(Momentum_recon::Mom_basetrack));
	}
	for ( int ilink = 0; ilink < num_link; ilink++ ) {
	  ofs.write((char*)& mom_chain.base_pair.at(ilink).first, sizeof(Momentum_recon::Mom_basetrack));
	  ofs.write((char*)& mom_chain.base_pair.at(ilink).second, sizeof(Momentum_recon::Mom_basetrack));
	}

      } // one particle emulsions

    } // reader

    ofs.close();

  } catch (const std::runtime_error &error) {
    BOOST_LOG_TRIVIAL(fatal) << "Runtime error : " << error.what();
    std::exit(1);
  } catch (const std::invalid_argument &error) {
    BOOST_LOG_TRIVIAL(fatal) << "Invalid argument error : " << error.what();
    std::exit(1);
  } catch (const std::out_of_range &error) {
    BOOST_LOG_TRIVIAL(fatal) << "Out of range error : " << error.what();
    std::exit(1);
  }

  BOOST_LOG_TRIVIAL(info) < "==========Finish==========";
  std::exit(0);

}
