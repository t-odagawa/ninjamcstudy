// system includes
#include <vector>
#include <array>
#include <algorithm>
#include <iostream>
#include <cmath>

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
#include <TRandom.h>
#include <TMath.h>
#include <TVector3.h>

#include <McsFunction.hpp>

#include "MyFunction.hpp"

namespace logging = boost::log;
namespace fs = boost::filesystem;

int main (int argc, char *argv[]) {

  logging::core::get()->set_filter
    (
      logging::trivial::severity >= logging::trivial::info
     //logging::trivial::severity >= logging::trivial::debug
     );

  if ( argc != 5 ) {
    BOOST_LOG_TRIVIAL(error) << "Usage : " << argv[0]
			     << " <input B2 file name> <input NTBM file name> <output file name> <material>";
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

    Int_t material = std::stoi(argv[4]);

    TFile *ofile = new TFile(argv[3], "recreate");
    TTree *otree = new TTree("tree", "tree");
    Int_t entry_in_original_file = 0;
    Double_t weight = 0.;
    Int_t muon_detection_flag = 0;
    std::vector<Double_t > true_vertex_position(3);
    Int_t vertex_ecc = -1;
    Int_t vertex_plate = -1;
    Double_t muon_true_angle = -1.;
    Double_t muon_recon_angle = -1.; // radial/lateral smear
    Double_t muon_true_momentum = -1.;
    Double_t muon_recon_momentum = -1.;
    std::vector<Double_t > muon_true_position(3);
    std::vector<Double_t > muon_recon_position(3); // accuracy
    std::vector<Double_t > muon_true_tangent(3);
    std::vector<Double_t > muon_recon_tangent(3);
    Int_t number_of_partners = 0;
    std::vector<Int_t > partner_pdg_code = {};
    std::vector<Int_t > partner_detect_flag = {};
    std::vector<Int_t > partner_one_seg_flag = {};
    std::vector<Int_t > partner_direction_sign = {};
    std::vector<Int_t > partner_plate = {};
    std::vector<Int_t > partner_second_plate = {};
    std::vector<Double_t > partner_true_angle = {};
    std::vector<Double_t > partner_recon_angle = {};
    std::vector<std::vector<Double_t > > partner_true_position = {};
    std::vector<std::vector<Double_t > > partner_recon_position = {};
    std::vector<std::vector<Double_t > > partner_true_tangent = {};
    std::vector<std::vector<Double_t > > partner_recon_tangent = {};
    std::vector<Double_t > partner_true_second_angle = {};
    std::vector<Double_t > partner_recon_second_angle = {};
    std::vector<std::vector<Double_t > > partner_true_second_position = {};
    std::vector<std::vector<Double_t > > partner_recon_second_position = {};
    std::vector<std::vector<Double_t > > partner_true_second_tangent = {};
    std::vector<std::vector<Double_t > > partner_recon_second_tangent = {};
    std::vector<Double_t > partner_true_momentum = {};
    std::vector<Double_t > partner_recon_momentum = {};
    std::vector<Double_t > extrapolate_distance = {};
    std::vector<Double_t > second_extrapolate_distance = {};
    std::vector<Double_t > minimum_distance = {};
    std::vector<Double_t > second_minimum_distance = {};
    std::vector<std::vector<Double_t > > partner_recon_vertex = {}; // middle point of MD
    std::vector<std::vector<Double_t > > partner_recon_second_vertex = {}; // middle point of MD

    otree->Branch("entry_in_original_file", &entry_in_original_file, "entry_in_original_file/I");
    otree->Branch("weight", &weight, "weight/D");
    otree->Branch("muon_detection_flag", &muon_detection_flag, "muon_detection_flag/B");
    otree->Branch("true_vertex_position", &true_vertex_position);
    otree->Branch("vertex_ecc", &vertex_ecc, "vertex_ecc/I");
    otree->Branch("vertex_plate", &vertex_plate, "vertex_plate/I");
    otree->Branch("muon_true_angle", &muon_true_angle, "muon_true_angle/D");
    otree->Branch("muon_recon_angle", &muon_recon_angle, "muon_recon_angle/D");
    otree->Branch("muon_true_momentum", &muon_true_momentum, "muon_true_momentum/D");
    otree->Branch("muon_recon_momentum", &muon_recon_momentum, "muon_recon_momentum/D");
    otree->Branch("muon_true_position", &muon_true_position);
    otree->Branch("muon_recon_position", &muon_recon_position);
    otree->Branch("muon_true_tangent", &muon_true_tangent);
    otree->Branch("muon_recon_tangent", &muon_recon_tangent);
    otree->Branch("number_of_partners", &number_of_partners, "number_of_partners/I");
    otree->Branch("partner_pdg_code", &partner_pdg_code);
    otree->Branch("partner_detect_flag", &partner_detect_flag);
    otree->Branch("partner_one_seg_flag", &partner_one_seg_flag);
    otree->Branch("partner_direction_sign", &partner_direction_sign);
    otree->Branch("partner_plate", &partner_plate);
    otree->Branch("partner_second_plate", &partner_second_plate);
    otree->Branch("partner_true_angle", &partner_true_angle);
    otree->Branch("partner_recon_angle", &partner_recon_angle);
    otree->Branch("partner_true_position", &partner_true_position);
    otree->Branch("partner_recon_position", &partner_recon_position);
    otree->Branch("partner_true_tangent", &partner_true_tangent);
    otree->Branch("partner_recon_tangent", &partner_recon_tangent);
    otree->Branch("partner_true_second_angle", &partner_true_second_angle);
    otree->Branch("partner_recon_second_angle", &partner_recon_second_angle);
    otree->Branch("partner_true_second_position", &partner_true_second_position);
    otree->Branch("partner_recon_second_position", &partner_recon_second_position);
    otree->Branch("partner_true_second_tangent", &partner_true_second_tangent);
    otree->Branch("partner_recon_second_tangent", &partner_recon_second_tangent);
    otree->Branch("partner_true_momentum", &partner_true_momentum);
    otree->Branch("partner_recon_momentum", &partner_recon_momentum);
    otree->Branch("extrapolate_distance", &extrapolate_distance);
    otree->Branch("second_extrapolate_distance", &second_extrapolate_distance);
    otree->Branch("minimum_distance", &minimum_distance);
    otree->Branch("second_minimum_distance", &second_minimum_distance);
    otree->Branch("partner_recon_vertex", &partner_recon_vertex);
    otree->Branch("partner_recon_second_vertex", &partner_recon_second_vertex);

    while ( reader.ReadNextSpill() > 0 ) {
      nt_tree->GetEntry(nt_entry);
      entry_in_original_file = nt_entry;
      nt_entry++;

      auto &input_spill_summary = reader.GetSpillSummary();

      auto it_event = input_spill_summary.BeginTrueEvent();
      const auto *event = it_event.Next(); // get true event summary
      const Double_t norm = event->GetNormalization();

      auto &primary_vertex_summary = event->GetPrimaryVertex();
      vertex_ecc = GetInteractionEcc(primary_vertex_summary);
      const Double_t total_cross_section = primary_vertex_summary.GetTotalCrossSection();
      weight = norm * total_cross_section * 1.e-38 * 6.02e23;

      TVector3 vertex_position = primary_vertex_summary.GetAbsolutePosition().GetValue();
      CalcPosInEccCoordinate(vertex_position, vertex_ecc, kTRUE);
      vertex_plate = GetVertexPlate(vertex_position.Z());

      // true information
      muon_detection_flag = 0;
      Int_t muon_track_id = -1;
      number_of_partners = 0;
      std::vector<Int_t > partner_track_id_vec;
      std::vector<Int_t > partner_pdg_code_vec;
      std::vector<Double_t > partner_true_momentum_vec;
      auto it_track = primary_vertex_summary.BeginTrack();
      while ( const auto *track = it_track.Next() ) {
	if ( B2Pdg::IsMuonPlusOrMinus(track->GetParticlePdg()) ) {
	  muon_track_id = (Int_t)track->GetTrackId();
	  muon_true_momentum = track->GetInitialMomentum().GetValue().Mag();
	  muon_detection_flag = 1;
	} else if ( B2Pdg::IsChargedPion(track->GetParticlePdg()) ||
		    track->GetParticlePdg() == PDG_t::kProton ) {
	  partner_track_id_vec.push_back((Int_t)track->GetTrackId());
	  partner_pdg_code_vec.push_back((Int_t)track->GetParticlePdg());
	  partner_true_momentum_vec.push_back(track->GetInitialMomentum().GetValue().Mag());
	  number_of_partners++;
	}
      } // track

      if ( number_of_partners <= 0 ) continue;
      partner_pdg_code.resize(number_of_partners);
      partner_detect_flag.resize(number_of_partners);
      partner_one_seg_flag.resize(number_of_partners);
      partner_direction_sign.resize(number_of_partners);
      partner_plate.resize(number_of_partners);
      partner_second_plate.resize(number_of_partners);
      partner_true_angle.resize(number_of_partners);
      partner_recon_angle.resize(number_of_partners);
      partner_true_position.resize(number_of_partners);
      partner_recon_position.resize(number_of_partners);
      partner_true_tangent.resize(number_of_partners);
      partner_recon_tangent.resize(number_of_partners);
      partner_true_second_angle.resize(number_of_partners);
      partner_recon_second_angle.resize(number_of_partners);
      partner_true_second_position.resize(number_of_partners);
      partner_recon_second_position.resize(number_of_partners);
      partner_true_second_tangent.resize(number_of_partners);
      partner_recon_second_tangent.resize(number_of_partners);
      partner_true_momentum.resize(number_of_partners);
      partner_recon_momentum.resize(number_of_partners);
      extrapolate_distance.resize(number_of_partners);
      second_extrapolate_distance.resize(number_of_partners);
      minimum_distance.resize(number_of_partners);
      second_minimum_distance.resize(number_of_partners);
      partner_recon_vertex.resize(number_of_partners);
      partner_recon_second_vertex.resize(number_of_partners);
      for ( Int_t ipartner = 0; ipartner < number_of_partners; ipartner++ ) {
	partner_detect_flag.at(ipartner) = 0;
	partner_true_position.at(ipartner).resize(3);
	partner_recon_position.at(ipartner).resize(3);
	partner_true_tangent.at(ipartner).resize(3);
	partner_recon_tangent.at(ipartner).resize(3);
	partner_true_second_position.at(ipartner).resize(3);
	partner_recon_second_position.at(ipartner).resize(3);
	partner_true_second_tangent.at(ipartner).resize(3);
	partner_recon_second_tangent.at(ipartner).resize(3);
	partner_recon_vertex.at(ipartner).resize(3);
	partner_recon_second_vertex.at(ipartner).resize(3);
      }
      if ( muon_detection_flag == 0 ) continue;
      BOOST_LOG_TRIVIAL(debug) << "CC interaction";

      // ISS/OSS/Shifter acceptance
      muon_detection_flag = 0;
      Bool_t iss_flag = false;
      Bool_t oss_flag = false;
      Bool_t shifter_flag = false;
      auto it_emulsion_acceptance = input_spill_summary.BeginEmulsion();
      while (const auto *emulsion = it_emulsion_acceptance.Next() ) {
	if ( emulsion->GetFilmType() == B2EmulsionType::kECC &&
	     emulsion->GetPlate() == 0 &&
	     emulsion->GetEcc() == vertex_ecc &&
	     emulsion->GetParentTrackId() == muon_track_id ) {
	  iss_flag = true;
	}
	else if ( emulsion->GetFilmType() == B2EmulsionType::kShifter &&
		  emulsion->GetPlate() == 0 &&
		  emulsion->GetEcc() == vertex_ecc &&
		  emulsion->GetParentTrackId() == muon_track_id ) 
	  oss_flag = true;
	else if ( emulsion->GetFilmType() == B2EmulsionType::kShifter &&
		  emulsion->GetPlate() == 17 &&
		  emulsion->GetParentTrackId() == muon_track_id ) {
	  shifter_flag = true;
	}
	else continue;
      }
      if ( iss_flag && oss_flag && shifter_flag )
	muon_detection_flag = 1;
      if ( muon_detection_flag == 0) continue;
      BOOST_LOG_TRIVIAL(debug) << "ISS/OSS/Shifter matching";

      // Tracker matching
      muon_detection_flag = 0;
      for ( Int_t icluster = 0; icluster< ntbm->GetNumberOfNinjaClusters(); icluster++ ) {
	if ( ntbm->GetNumberOfHits(icluster, B2View::kSideView) > 0 &&
	     ntbm->GetNumberOfHits(icluster, B2View::kTopView) > 0 ) {
	  muon_detection_flag = 1;
	  break;
	} // 2d cluster
      } // icluster
      if ( muon_detection_flag == 0) continue;
      BOOST_LOG_TRIVIAL(debug) << "Tracker matching";

      // search for relevant basetracks
      TVector3 muon_position_smeared; // used in minimum distance calculation so more global
      TVector3 muon_tangent_smeared;
      std::vector<TVector3 > partner_position_smeared(number_of_partners);
      std::vector<TVector3 > partner_second_position_smeared(number_of_partners);
      std::vector<TVector3 > partner_tangent_smeared(number_of_partners);
      std::vector<TVector3 > partner_second_tangent_smeared(number_of_partners);
      auto it_emulsion = input_spill_summary.BeginEmulsion();
      while ( const auto *emulsion = it_emulsion.Next() ) {
	if ( emulsion->GetFilmType() != B2EmulsionType::kECC ) continue;
	if ( emulsion->GetEcc() != vertex_ecc ) continue;
	if ( emulsion->GetPlate() - vertex_plate != -1 &&
	     emulsion->GetPlate() - vertex_plate != 0 &&
	     emulsion->GetPlate() - vertex_plate != 1 &&
	     emulsion->GetPlate() - vertex_plate != 2 ) continue;
	// muon basetrack
	if ( emulsion->GetParentTrackId() == muon_track_id &&
	     emulsion->GetPlate() == vertex_plate ) {
	  TVector3 position = emulsion->GetAbsolutePosition().GetValue();
	  CalcPosInEccCoordinate(position, vertex_ecc, kTRUE);
	  muon_true_position.at(0) = position.X();
	  muon_true_position.at(1) = position.Y();
	  muon_true_position.at(2) = position.Z();
	  TVector3 tangent = emulsion->GetTangent().GetValue();
	  muon_true_angle = std::hypot(tangent.X(), tangent.Y());
	  muon_true_tangent.at(0) = tangent.X();
	  muon_true_tangent.at(1) = tangent.Y();
	  muon_true_tangent.at(2) = tangent.Z();
	  // smear
	  muon_position_smeared = position;
	  SmearPosition(muon_position_smeared);
	  muon_recon_position.at(0) = muon_position_smeared.X();
	  muon_recon_position.at(1) = muon_position_smeared.Y();
	  muon_recon_position.at(2) = muon_position_smeared.Z();
	  muon_tangent_smeared = tangent;
	  SmearTangentVector(tangent);
	  muon_recon_angle = std::hypot(muon_tangent_smeared.X(), muon_tangent_smeared.Y());
	  muon_recon_tangent.at(0) = muon_tangent_smeared.X();
	  muon_recon_tangent.at(1) = muon_tangent_smeared.Y();
	  muon_recon_tangent.at(2) = muon_tangent_smeared.Z();
	  continue; // go to next basetrack after muon process
	}

	// partner basetrack
	auto partner_emulsion_track_id = std::find(partner_track_id_vec.begin(),
						   partner_track_id_vec.end(),
						   emulsion->GetParentTrackId());
	if ( partner_emulsion_track_id != partner_track_id_vec.end() ) {
	  Int_t partner_index = partner_emulsion_track_id - partner_track_id_vec.begin();
	  if ( partner_detect_flag.at(partner_index) == 0 ) {
	    partner_pdg_code.at(partner_index) = emulsion->GetParentTrack().GetParticlePdg();
	    partner_one_seg_flag.at(partner_index) = 1;
	    partner_direction_sign.at(partner_index) = 0;
	    partner_detect_flag.at(partner_index) = 1;
	  }	  
	  TVector3 position = emulsion->GetAbsolutePosition().GetValue();
	  CalcPosInEccCoordinate(position, vertex_ecc, kTRUE);
	  TVector3 tangent = emulsion->GetTangent().GetValue();
	  if ( emulsion->GetPlate() - vertex_plate == 0 ) { 
	    partner_plate.at(partner_index) = emulsion->GetPlate() + 1;
	    std::vector<Double_t > true_position_tmp; true_position_tmp.resize(3);
	    true_position_tmp.at(0) = position.X();
	    true_position_tmp.at(1) = position.Y();
	    true_position_tmp.at(2) = position.Z();
	    partner_true_position.at(partner_index) = true_position_tmp;
	    partner_true_angle.at(partner_index) = std::hypot(tangent.X(), tangent.Y());
	    partner_true_tangent.at(partner_index).at(0) = tangent.X();
	    partner_true_tangent.at(partner_index).at(1) = tangent.Y();
	    partner_true_tangent.at(partner_index).at(2) = tangent.Z();
	    // smear
	    partner_position_smeared.at(partner_index) = position;
	    SmearPosition(partner_position_smeared.at(partner_index));
	    std::vector<Double_t > recon_position_tmp; recon_position_tmp.resize(3);
	    recon_position_tmp.at(0) = partner_position_smeared.at(partner_index).X();
	    recon_position_tmp.at(1) = partner_position_smeared.at(partner_index).Y();
	    recon_position_tmp.at(2) = partner_position_smeared.at(partner_index).Z();
	    partner_recon_position.at(partner_index) = recon_position_tmp;
	    partner_tangent_smeared.at(partner_index) = tangent;
	    SmearTangentVector(partner_tangent_smeared.at(partner_index));
	    partner_recon_angle.at(partner_index) = std::hypot(partner_tangent_smeared.at(partner_index).X(),
							       partner_tangent_smeared.at(partner_index).Y());
	    partner_recon_tangent.at(partner_index).at(0) = partner_tangent_smeared.at(partner_index).X();
	    partner_recon_tangent.at(partner_index).at(1) = partner_tangent_smeared.at(partner_index).Y();
	    partner_recon_tangent.at(partner_index).at(2) = partner_tangent_smeared.at(partner_index).Z();
	    // momentum
	    partner_true_momentum.at(partner_index) = partner_true_momentum_vec.at(partner_index);
	    partner_direction_sign.at(partner_index) = 1;
	  }
	  else if ( emulsion->GetPlate() - vertex_plate == 1 ){
	    partner_plate.at(partner_index) = emulsion->GetPlate() + 1;
	    std::vector<Double_t > true_position_tmp; true_position_tmp.resize(3);
	    true_position_tmp.at(0) = position.X();
	    true_position_tmp.at(1) = position.Y();
	    true_position_tmp.at(2) = position.Z();
	    partner_true_position.at(partner_index) = true_position_tmp;
	    partner_true_angle.at(partner_index) = std::hypot(tangent.X(), tangent.Y());
	    partner_true_tangent.at(partner_index).at(0) = tangent.X();
	    partner_true_tangent.at(partner_index).at(1) = tangent.Y();
	    partner_true_tangent.at(partner_index).at(2) = tangent.Z();
	    // smear
	    partner_position_smeared.at(partner_index) = position;
	    SmearPosition(partner_position_smeared.at(partner_index));
	    std::vector<Double_t > recon_position_tmp; recon_position_tmp.resize(3);
	    recon_position_tmp.at(0) = partner_position_smeared.at(partner_index).X();
	    recon_position_tmp.at(1) = partner_position_smeared.at(partner_index).Y();
	    recon_position_tmp.at(2) = partner_position_smeared.at(partner_index).Z();
	    partner_recon_position.at(partner_index) = recon_position_tmp;
	    partner_tangent_smeared.at(partner_index) = tangent;
	    SmearTangentVector(partner_tangent_smeared.at(partner_index));
	    partner_recon_angle.at(partner_index) = std::hypot(partner_tangent_smeared.at(partner_index).X(),
							       partner_tangent_smeared.at(partner_index).Y());
	    partner_recon_tangent.at(partner_index).at(0) = partner_tangent_smeared.at(partner_index).X();
	    partner_recon_tangent.at(partner_index).at(1) = partner_tangent_smeared.at(partner_index).Y();
	    partner_recon_tangent.at(partner_index).at(2) = partner_tangent_smeared.at(partner_index).Z();
	    // momentum
	    partner_true_momentum.at(partner_index) = partner_true_momentum_vec.at(partner_index);
	    partner_direction_sign.at(partner_index) = -1;
	  }	    
	  else if ( emulsion->GetPlate() - vertex_plate == -1 ) {
	    partner_second_plate.at(partner_index) = emulsion->GetPlate() + 1;
	    std::vector<Double_t > true_position_tmp; true_position_tmp.resize(3);
	    true_position_tmp.at(0) = position.X();
	    true_position_tmp.at(1) = position.Y();
	    true_position_tmp.at(2) = position.Z();
	    partner_true_second_position.at(partner_index) = true_position_tmp;
	    partner_true_second_angle.at(partner_index) = std::hypot(tangent.X(), tangent.Y());
	    partner_true_second_tangent.at(partner_index).at(0) = tangent.X();
	    partner_true_second_tangent.at(partner_index).at(1) = tangent.Y();
	    partner_true_second_tangent.at(partner_index).at(2) = tangent.Z();
	    // smear
	    partner_second_position_smeared.at(partner_index) = position;
	    SmearPosition(partner_second_position_smeared.at(partner_index));
	    std::vector<Double_t > recon_position_tmp; recon_position_tmp.resize(3);
	    recon_position_tmp.at(0) = partner_second_position_smeared.at(partner_index).X();
	    recon_position_tmp.at(1) = partner_second_position_smeared.at(partner_index).Y();
	    recon_position_tmp.at(2) = partner_second_position_smeared.at(partner_index).Z();
	    partner_recon_second_position.at(partner_index) = recon_position_tmp;
	    partner_second_tangent_smeared.at(partner_index) = tangent;
	    SmearTangentVector(partner_second_tangent_smeared.at(partner_index));
	    partner_recon_second_angle.at(partner_index) = std::hypot(partner_second_tangent_smeared.at(partner_index).X(),
								      partner_second_tangent_smeared.at(partner_index).Y());
	    partner_recon_second_tangent.at(partner_index).at(0) = partner_second_tangent_smeared.at(partner_index).X();
	    partner_recon_second_tangent.at(partner_index).at(1) = partner_second_tangent_smeared.at(partner_index).Y();
	    partner_recon_second_tangent.at(partner_index).at(2) = partner_second_tangent_smeared.at(partner_index).Z();
	    partner_one_seg_flag.at(partner_index) = 0;
	  }
	  else if ( emulsion->GetPlate() - vertex_plate == 2 ) {
	    partner_second_plate.at(partner_index) = emulsion->GetPlate() + 1;
	    std::vector<Double_t > true_position_tmp; true_position_tmp.resize(3);
	    true_position_tmp.at(0) = position.X();
	    true_position_tmp.at(1) = position.Y();
	    true_position_tmp.at(2) = position.Z();
	    partner_true_second_position.at(partner_index) = true_position_tmp;
	    partner_true_second_angle.at(partner_index) = std::hypot(tangent.X(), tangent.Y());
	    partner_true_second_tangent.at(partner_index).at(0) = tangent.X();
	    partner_true_second_tangent.at(partner_index).at(1) = tangent.Y();
	    partner_true_second_tangent.at(partner_index).at(2) = tangent.Z();
	    // smear
	    partner_second_position_smeared.at(partner_index) = position;
	    SmearPosition(partner_second_position_smeared.at(partner_index));
	    std::vector<Double_t > recon_position_tmp; recon_position_tmp.resize(3);
	    recon_position_tmp.at(0) = partner_second_position_smeared.at(partner_index).X();
	    recon_position_tmp.at(1) = partner_second_position_smeared.at(partner_index).Y();
	    recon_position_tmp.at(2) = partner_second_position_smeared.at(partner_index).Z();
	    partner_recon_second_position.at(partner_index) = recon_position_tmp;
	    partner_second_tangent_smeared.at(partner_index) = tangent;
	    SmearTangentVector(partner_second_tangent_smeared.at(partner_index));
	    partner_recon_second_angle.at(partner_index) = std::hypot(partner_second_tangent_smeared.at(partner_index).X(),
								      partner_second_tangent_smeared.at(partner_index).Y());
	    partner_recon_second_tangent.at(partner_index).at(0) = partner_second_tangent_smeared.at(partner_index).X();
	    partner_recon_second_tangent.at(partner_index).at(1) = partner_second_tangent_smeared.at(partner_index).Y();
	    partner_recon_second_tangent.at(partner_index).at(2) = partner_second_tangent_smeared.at(partner_index).Z();
	    partner_one_seg_flag.at(partner_index) = 0;
	  }
	}
	else continue;
      }


      // minimum distance check
      std::array<Double_t, 2 > z_range = {};
      if ( material == B2Material::kWater ) {
	z_range.at(0) = (- NINJA_WATER_LAYER_THICK - NINJA_FILM_THICK - 2 * NINJA_ENV_THICK) * 1.e3 - 1000.;
	z_range.at(1) = (+ NINJA_BASE_LAYER_THICK + NINJA_EMULSION_LAYER_THICK) * 1.e3 + 30.;
      }
      else if ( material == B2Material::kIron ) {
	z_range.at(0) = (- NINJA_IRON_LAYER_THICK - NINJA_FILM_THICK) * 1.e3 - 50.;
	z_range.at(1) = (+ NINJA_BASE_LAYER_THICK + NINJA_EMULSION_LAYER_THICK) * 1.e3 + 30.;
      }
      else
	throw std::invalid_argument("Material is not valid");
      z_range.at(0) += muon_position_smeared.Z();
      z_range.at(1) += muon_position_smeared.Z();
      std::array<Double_t, 2 > extrapolate_z = {};
      for ( Int_t ipartner = 0; ipartner < number_of_partners; ipartner++) {
	if ( !partner_detect_flag.at(ipartner) ) continue;
	minimum_distance.at(ipartner) = GetMinimumDistance(muon_position_smeared,
							   partner_position_smeared.at(ipartner),
							   muon_tangent_smeared,
							   partner_tangent_smeared.at(ipartner),
							   z_range, extrapolate_z,
							   partner_recon_vertex.at(ipartner));
	extrapolate_distance.at(ipartner) = extrapolate_z.at(0);
	
	if (!partner_one_seg_flag.at(ipartner)) {
	  second_minimum_distance.at(ipartner) = GetMinimumDistance(muon_position_smeared,
								    partner_second_position_smeared.at(ipartner),
								    muon_tangent_smeared,
								    partner_second_tangent_smeared.at(ipartner),
								    z_range, extrapolate_z,
								    partner_recon_second_vertex.at(ipartner));
	  second_extrapolate_distance.at(ipartner) = extrapolate_z.at(0);
	  }
      }

      true_vertex_position.at(0) = vertex_position.X();
      true_vertex_position.at(1) = vertex_position.Y();
      true_vertex_position.at(2) = vertex_position.Z();
      vertex_ecc = vertex_ecc + 1;
      vertex_plate = vertex_plate + 1;

      otree->Fill();

      entry_in_original_file = -1;
      weight = 0.;
      muon_detection_flag = 0;
      true_vertex_position.clear(); true_vertex_position.resize(3);
      vertex_ecc = -1;
      vertex_plate = -1;
      muon_true_angle = -1.;
      muon_recon_angle = -1.;
      muon_true_momentum = -1.;
      muon_recon_momentum = -1.;
      muon_true_position.clear(); muon_true_position.resize(3);
      muon_recon_position.clear(); muon_recon_position.resize(3);
      number_of_partners = 0;
      partner_pdg_code.clear();
      partner_detect_flag.clear();
      partner_one_seg_flag.clear();
      partner_direction_sign.clear();
      partner_plate.clear();
      partner_second_plate.clear();
      partner_true_angle.clear();
      partner_recon_angle.clear();
      partner_true_position.clear();
      partner_recon_position.clear();
      partner_true_second_angle.clear();
      partner_recon_second_angle.clear();
      partner_true_second_position.clear();
      partner_recon_second_position.clear();
      partner_true_momentum.clear();
      partner_recon_momentum.clear();
      extrapolate_distance.clear();
      second_extrapolate_distance.clear();
      minimum_distance.clear();
      second_minimum_distance.clear();
      partner_recon_vertex.clear();
      partner_recon_second_vertex.clear();
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
  } catch (const std::out_of_range &error) {
    BOOST_LOG_TRIVIAL(fatal) << "Out of range error : " << error.what();
    std::exit(1);
  }

  BOOST_LOG_TRIVIAL(info) << "==========Finish==========";
  std::exit(0);

}
