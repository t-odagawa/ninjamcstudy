// boost include
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

// B2includes
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

Int_t GetVertexPlate(Double_t z_pos /*um*/) {
  // um -> mm
  z_pos /= 1.e3;

}

void CalcPosInEccCoordinate(TVector3 &vertex_position, Int_t ecc_id) {
  // center of ECC5 dessicator
  vertex_position.SetX(vertex_position.X() - NINJA_POS_X - NINJA_ECC_POS_X);
  vertex_position.SetY(vertex_position.Y() - NINJA_POS_Y - NINJA_ECC_POS_Y);
  vertex_position.SetZ(vertex_position.Z() - NINJA_POS_Z - NINJA_ECC_POS_Z);
  // film coordinate
  vertex_position.SetX(vertex_position.X()
		       + 0.5 * NINJA_ECC_FILM_XY);
  vertex_position.SetY(vertex_position.Y()
		       + 0.5 * NINJA_DESIC_HEIGHT
		       - NINJA_DESIC_THICK
		       - NINJA_ENV_THICK);
  vertex_position.SetZ(vertex_position.Z()
		       - 0.5 * NINJA_DESIC_DEPTH
		       + NINJA_DESIC_THICK
		       + NINJA_ENV_THICK
		       + NINJA_EMULSION_LAYER_THICK
		       + NINJA_BASE_LAYER_THICK);
  // move to each ECC
  vertex_position.SetX(vertex_position.X()
		       + NINJA_ECC_GAP_X * (1 - ecc_id % 3));
  vertex_position.SetY(vertex_position.Y()
		       + NINJA_ECC_GAP_Y * (ecc_id / 3 - 1));
  // mm -> um
  vertex_position.SetX(vertex_position.X() * 1.e3);
  vertex_position.SetY(vertex_position.Y() * 1.e3);
  vertex_position.SetZ(vertex_position.Z() * 1.e3);
}

TVector3 SmearPosition(TVector3 true_position /*um*/ ) {

}

TVector3 SmearTangent(TVector3 true_tangent) {

}

Double_t GetMinimumDistance(TVector3 parent_pos, TVector3 daughter_pos, TVector3 parent_dir, TVector3 daughter_dir,
			    std::array<Double_t,2> z_range, std::array<Double_t,2> &extrapolate_z) {
  std::array<Double_t,2> extrapolate_distance;
  TVector3 position_difference = parent_pos - daughter_pos;
  // Almost parallel
  if ( TMath::ACos((parent_dir * daughter_dir) / (parent_dir.Mag() * daughter_dir.Mag())) < 1.e-4) {
    extrapolate_distance.at(0) = (parent_pos.Z() + daughter_pos.Z()) / 2. - parent_pos.Z();
    extrapolate_distance.at(1) = (parent_pos.Z() + daughter_pos.Z()) / 2. - daughter_pos.Z();
  }
  else {
    Double_t delta = parent_dir.Mag2() * daughter_dir.Mag2() - (parent_dir * daughter_dir) * (parent_dir * daughter_dir);
    extrapolate_distance.at(0) = ( 1 * (position_difference * parent_dir) * daughter_dir.Mag2()
				   - (parent_dir * daughter_dir) * (position_difference * daughter_dir)) / delta;
    extrapolate_distance.at(1) = (-1 * (position_difference * daughter_dir) * parent_dir.Mag2()
				   + (parent_dir * daughter_dir) * (position_difference * parent_dir)) / delta;
  }
  // z_range.at(0) : small, z_range.at(1) : large
  if ( z_range.at(0) > z_range.at(1) ) {
    std::swap(z_range.at(0), z_range.at(1));
  }
  if ( parent_pos.Z() + extrapolate_distance.at(0) < z_range.at(0) ||
       daughter_pos.Z() + extrapolate_distnace.at(1) < z_range.at(0)) {
    extrapolate_distance.at(0) = z_range.at(0) - parent_pos.Z();
    extrapolate_distance.at(1) = z_range.at(0) - daughter_pos.Z();
  }
  else if ( parent_pos.Z() + extrapolate_distance.at(0) > z_range.at(1) ||
	    daughter_pos.Z() + extrapolate_distance.at(1) > z_range.at(1)) {
    extrapolate_distance.at(0) = z_range.at(1) - parent_pos.Z();
    extrapolate_distance.at(1) - z_range.at(1) - daughter_pos.Z();
  }

  extrapolate_z.at(0) = extrapolate_distance.at(0);
  extrapolate_z.at(1) = extrapolate_distance.at(1);

  TVector3 calculate_parent_position = parent_pos + extrapolate_distance.at(0) * parnet_dir;
  TVector3 calculate_daughter_position = daughter_pos + extrapolate_distance.at(1) * daughter_dir;

  TVector3 distance_vec = calculate_parent_position - calculate_daughter_position;

  return distance_vec.Mag();

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
    TTree *otree = new TTree("tree", "tree");
    Double_t weight;
    std::vector<Double_t > true_vertex_position(3);
    Int_t vertex_ecc;
    Int_t vertex_plate;
    Double_t muon_true_angle;
    Double_t muon_recon_angle; // radial/lateral smear
    Double_t muon_true_momentum;
    Double_t muon_recon_momentum;
    std::vector<Double_t > muon_true_position(3);
    std::vector<Double_t > muon_recon_position(3); // accuracy
    Int_t number_of_partners;
    std::vector<Int_t > partner_pdg_code;
    std::vector<Bool_t > partner_one_seg_flag;
    std::vector<Int_t > partner_direction_sign;
    std::vector<Int_t > partner_plate;
    std::vector<Int_t > partner_second_plate;
    std::vector<Double_t > partner_true_angle;
    std::vector<Double_t > partner_recon_angle;
    std::vector<std::vector<Double_t > > partner_true_position;
    std::vector<std::vector<Double_t > > partner_recon_position;
    std::vector<Double_t > partner_true_second_angle;
    std::vector<Double_t > partner_recon_second_angle;
    std::vector<std::vector<Double_t > > partner_true_second_position;
    std::vector<std::vector<Double_t > > partner_recon_second_position;    
    std::vector<Double_t > partner_true_momentum;
    std::vector<Double_t > partner_recon_momentum;
    std::vector<Double_t > extrapolate_distance;
    std::vector<Double_t > second_extrapolate_distance;
    std::vector<Double_t > minimum_distance;
    std::vector<Double_t > second_minimum_distance;
    std::vector<std::vector<Double_t > > partner_recon_vertex; // middle point of MD

    otree->Branch("weight", &weight, "weight/D");
    otree->Branch("true_vertex_position", &true_vertex_position);
    otree->Branch("vertex_ecc", &vertex_ecc, "vertex_ecc/I");
    otree->Branch("vertex_plate", &vertex_plate, "vertex_plate/I");
    otree->Branch("muon_true_angle", &muon_true_angle, "muon_true_angle/D");
    otree->Branch("muon_recon_angle", &muon_recon_angle, "muon_recon_angle/D");
    otree->Branch("muon_true_momentum", &muon_true_momentum, "muon_true_momentum/D");
    otree->Branch("muon_recon_momentum", &muon_recon_momentum, "muon_recon_momentum/D");
    otree->Branch("muon_true_position", &muon_true_position);
    otree->Branch("muon_recon_position", &muon_recon_position);
    otree->Branch("number_of_partners", &number_of_partners, "number_of_partners/I");
    otree->Branch("partner_pdg_code", &partner_pdg_code);
    otree->Branch("partner_one_seg_flag", &partner_one_seg_flag);
    otree->Branch("partner_direction_sign", &partner_direction_sign);
    otree->Branch("partner_plate", &partner_plate);
    otree->Branch("partner_second_plate", &partner_second_plate);
    otree->Branch("partner_true_angle", &partner_true_angle);
    otree->Branch("partner_recon_angle", &partner_recon_angle);
    otree->Branch("partner_true_position", &partner_true_position);
    otree->Branch("partner_recon_position", &partner_recon_position);
    otree->Branch("partner_true_second_angle", &partner_true_second_angle);
    otree->Branch("partner_recon_second_angle", &partner_recon_second_angle);
    otree->Branch("partner_true_second_position", &partner_true_second_position);
    otree->Branch("partner_recon_second_position", &partner_recon_second_position);
    otree->Branch("partner_true_momentum", &partner_true_momantum);
    otree->Branch("partner_recon_momentum", &partner_recon_momentum);
    otree->Branch("extrapolate_distance", &extrapolate_distance);
    otree->Branch("second_extrapolate_distance", &seconde_extrapolate_distance);
    otree->Branch("minimum_distance", &minimum_distance);
    otree->Branch("second_minimum_distance", &second_minimum_distance);
    otree->Branch("partner_recon_vertex", &partner_recon_vertex);

    while ( reader.ReadNextSpill() > 0 ) {
      nt_tree->GetEntry(nt_entry);

      auto &input_spill_summary = reader.GetSpillSummary();

      auto it_event = input_spill_summary.BeginTrueEvent();
      const auto *event = it_event.Next(); // get true event summary
      const Double_t norm = event->GetNormalization();

      auto &primary_vertex_summary = event->GetPrimaryVertex();
      vertex_ecc = GetInteractionEcc(primary_vertex_summary);
      const Double_t total_cross_section = primary_vetex_summary.GetTotalCrossSection();
      weight = norm * total_cross_section * 1.e-38 * 6.02e23;

      TVector3 vertex_position = primary_vertex_summary.GetAbsolutePosition().GetValue();
      CalcPosInEccCoordinate(vertex_position, vertex_ecc);
      vertex_plate = GetVertexPlateId(vertex_position.Z());

      // true information
      Bool_t cc_interaction_flag = false;
      Int_t muon_track_id = -1;
      std::vector<Int_t > partner_track_id_vec;
      std::vector<Int_t > partner_pdg_code_vec;
      std::vector<Double_t > partner_true_momentum_vec;
      auto it_track = primary_vertex_summary.BeginTrack();
      while ( *track = it_track.Next() ) {
	if ( B2Pdg::IsMuonPlusOrMinus(track->GetParticlePdg()) ) {
	  muon_track_id = (Int_t)track->GetTrackId();
	  muon_true_momentum = track->GetInitialMomentum().GetValue().Mag();
	  cc_interaction_flag = true;
	} else if ( B2Pdg::IsChargedPion(track->GetParticlePdg()) ||
		    track->GetParticlePdg() == PDG_t::kProton ) {
	  partner_track_id_vec.push_back((Int_t)track->GetTrackId());
	  partner_pdg_code_vec.push_back((Int_t)track->GetParticlePdg());
	  partner_true_momentum_vec.push_back(track->GetInitialMomentum().GetValue().Mag());
	}
      } // track

      // muon matching check
      if (!cc_interaction_flag) continue;
      for ( Int_t icluster = 0; icluster< ntbm->GetNumberOfNinjaClusters(); icluster++ ) {
	if ( ntbm->GetNumberOfHits(icluster, B2View::kSideView) > 0 &&
	     ntbm->GetNumberOfHits(icluster, B2View::kTopView) > 0 ) {
	  break;
	} // 2d cluster
      } // icluster
      // ISS/OSS/Shifter acceptance

      // search for relevant basetracks
      TVector3 muon_position_smeared;
      TVector3 muon_second_position_smeared;
      TVector3 muon_tangent_smeared;
      TVector3 muon_second_tangent_smeared;
      number_of_partners = 0;
      std::vector<TVector3 > partner_position_smeared;
      std::vector<TVector3 > partner_second_position_smeared;
      std::vector<TVector3 > partner_tangent_smeared;
      std::vector<TVector3 > partner_second_tangent_smeared;
      auto it_emulsion = input_spill_summary.BeginEmulsion();
      while ( const auto *emulsion = it_emulsion.Next() ) {
	if ( emulsion->GetFilmType() != B2EmulsionType::kECC ) continue;
	if ( std::fabs(emulsion->GetPlate() - vertex_plate) > 2 ) continue;
	// muon basetrack
	if ( emulsion->GetParentTrackId() == muon_track_id ) {
	  TVector3 position = emulsion->GetAbsolutePosition().GetValue();
	  CalcPosInEccCoordinate(position, vertex_ecc);
	  TVector3 tangent = emulsion->GetTangent().GetValue();
	  if ( emulsion->GetPlate() == vertex_plate + 1 ) {
	    muon_position_smeared = SmearPosition(position);
	    muon_tangent_smeared = SmearTangent(tangent);
	  }
	  else if ( emulsion->GetPlate() == vertex_plate + 2 ) {
	    muon_second_position_smeared = SmearPosition(position);
	    muon_second_tangent_smeared = SmearTangent(tangent);
	  }
	}
	// partner basetrack
	auto partner_emulsion = std::find(partner_track_id_vec.begin(),
					  partner_track_id_vec.end(),
					  emulsion->GetParentTrackId());
	else if ( partner_emulsion != partner_track_id_vec.end() ) {
	  number_of_partners++;
	  TVector3 position = emulsion->GetAbsolutePosition().GetValue();
	  CalcPosInEccCoordinate(position, vertex_ecc);
	}

      }

      // minimum distance check
      std::array<Double_t, 2 > z_range = {};
      if ( material == B2Material::kWater ) {
	z_range.at(0) = - NINJA_WATER_LAYER_THICK - NINJA_FILM_THICK - 2 * NINJA_ENV_THICK - 1000.;
	z_range.at(1) = + NINJA_BASE_LAYER_THICK + NINJA_EMULSION_LAYER_THICK + 30.;
      }
      else if ( material == B2Material::kIron ) {
	z_range.at(0) = - NINJA_IRON_LAYER_THICK - NINJA_FILM_THICK - 50.;
	z_range.at(1) = + NINJA_BASE_LAYER_THICK + NINJA_EMULSION_LAYER_THICK + 30.;
      }
      z_range.at(0) += muon_position_smeared.Z();
      z_range.at(1) += muon_position_smeared.Z();
      std::array<Double_t, 2 > extrapolate_z;
      for ( Int_t ipartner = 0; ipartner < number_of_partners; ipartner++) {
	minimum_distance.push_back(GetMinimumDistance(muon_position_smeared,
						      partner_position_smeared.at(ipartner),
						      muon_tangent_smeared,
						      partner_tangent_smeared.at(ipartner),
						      z_range, extrapolate_z));
	if (!one_seg_flag) {
	  second_minimum_distance.push_back(GetMinimumDistance(muon_position_smeared,
							       partner_second_position_smeared.at(ipartner),
							       muon_tangent_smeared,
							       partner_second_tangent_smeared.at(ipartner),
							       z_range, extrapolate_z));
	}
      }

      true_vertex_position.at(0) = vertex_position.X();
      true_vertex_position.at(1) = vertex_position.Y();
      true_vertex_position.at(2) = vertex_position.Z();
      vertex_ecc = vertex_ecc + 1;
      vertex_plate = vertex_plate + 1;

      otree->Fill();
      nt_entry++;
    }

    ofile->cd();
    otree->Write();
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
