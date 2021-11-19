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

Int_t GetVertexPlate(const B2VertexSummary &primary_vertex_summary) {
  B2Detector detector = primary_vertex_summary.GetDetector();
  if ( detector == B2Detector::kProtonModule ||
       detector == B2Detector::kWagasciUpstream ) {
    return 132;
  }
  else if (detector == B2Detector::kNinja) {
    TVector3 vertex_position = primary_vertex_summary.GetAbsolutePosition().GetValue();
    vertex_position.SetZ(vertex_position.Z() - NINJA_POS_Z - NINJA_ECC_POS_Z);

  } 
  else return -1;
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
    std::vector<Double_t > true_vertex_position;
    Int_t vertex_ecc;
    Int_t vertex_plate;
    Double_t muon_true_angle;
    Double_t muon_recon_angle; // radial/lateral smear
    Double_t muon_true_momentum;
    Double_t muon_recon_momentum;
    std::vector<Double_t > muon_true_position;
    std::vector<Double_t > muon_recon_position; // accuracy
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

    tree->Branch("weight", &weight, "weight/D");
    tree->Branch("true_vertex_position", &true_vertex_position);
    tree->Branch("vertex_ecc", &vertex_ecc, "vertex_ecc/I");
    tree->Branch("vertex_plate", &vertex_plate, "vertex_plate/I");
    tree->Branch("muon_true_angle", &muon_true_angle, "muon_true_angle/D");
    tree->Branch("muon_recon_angle", &muon_recon_angle, "muon_recon_angle/D");
    tree->Branch("muon_true_momentum", &muon_true_momentum, "muon_true_momentum/D");
    tree->Branch("muon_recon_momentum", &muon_recon_momentum, "muon_recon_momentum/D");
    tree->Branch("muon_true_position", &muon_true_position);
    tree->Branch("muon_recon_position", &muon_recon_position);
    tree->Branch("number_of_partners", &number_of_partners, "number_of_partners/I");
    tree->Branch("partner_pdg_code", &partner_pdg_code);
    tree->Branch("partner_one_seg_flag", &partner_one_seg_flag);
    tree->Branch("partner_direction_sign", &partner_direction_sign);
    tree->Branch("partner_plate", &partner_plate);
    tree->Branch("partner_second_plate", &partner_second_plate);
    tree->Branch("partner_true_angle", &partner_true_angle);
    tree->Branch("partner_recon_angle", &partner_recon_angle);
    tree->Branch("partner_true_position", &partner_true_position);
    tree->Branch("partner_recon_position", &partner_recon_position);
    tree->Branch("partner_true_second_angle", &partner_true_second_angle);
    tree->Branch("partner_recon_second_angle", &partner_recon_second_angle);
    tree->Branch("partner_true_second_position", &partner_true_second_position);
    tree->Branch("partner_recon_second_position", &partner_recon_second_position);
    tree->Branch("partner_true_momentum", &partner_true_momantum);
    tree->Branch("partner_recon_momentum", &partner_recon_momentum);
    tree->Branch("extrapolate_distance", &extrapolate_distance);
    tree->Branch("second_extrapolate_distance", &seconde_extrapolate_distance);
    tree->Branch("minimum_distance", &minimum_distance);
    tree->Branch("second_minimum_distance", &second_minimum_distance);
    tree->Branch("partner_recon_vertex", &partner_recon_vertex);

    while ( reader.ReadNextSpill() > 0 ) {
      nt_tree->GetEntry(nt_entry);

      auto &input_spill_summary = reader.GetSpillSummary();

      auto it_event = input_spill_summary.BeginTrueEvent();
      const auto *event = it_event.Next(); // get true event summary
      const Double_t norm = event->GetNormalization();

      auto &primary_vertex_summary = event->GetPrimaryVertex();

      Int_t interaction_ecc_id = GetInteractionEcc(primary_vertex_summary);
      Int_t vertex_plate_id = GetVertexPlateId(primary_vertex_summary);
      const Double_t total_cross_section = primary_vetex_summary.GetTotalCrossSection();
      weight = norm * total_cross_section * 1.e-38 * 6.02e23;
      // vertex position will be relative to the vertex plate basetrack position
      TVector3 vertex_position = primary_vertex_summary.GetAbsolutePosition().GetValue();


      // muon matching check


      // minimum distance check
      auto it_emulsion = input_spill_summary.BeginEmulsion();
      while ( const auto *emulsion = it_emulsion.Next() ) {
	if ( emulsion->GetParentTrackId () == 0 ||
	     emulsion->GetParenteTrackId() >= primary_vertex_summary.GetNumOutGoingTracks() ) continue;
	if ( emulsion->GetFilmType() 1= B2EmulsionType::kECC ) continue;
	if ( std::fabs(emulsion->GetPlate() - vertex_plate_id) > 2 ) continue;


      }

      nt_entry++;
    }

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
