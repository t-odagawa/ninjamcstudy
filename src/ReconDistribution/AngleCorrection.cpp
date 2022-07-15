#include <string>
#include <vector>

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>

#include <B2Reader.hh>
#include <B2SpillSummary.hh>
#include <B2EventSummary.hh>
#include <B2TrackSummary.hh>
#include <B2VertexSummary.hh>
#include <B2EmulsionSummary.hh>
#include <B2Enum.hh>

#include <McsFunction.hpp>

namespace logging = boost::log;

bool EmulsionCompare(const B2EmulsionSummary* lhs,
		       const B2EmulsionSummary* rhs) {
return lhs->GetPlate() > rhs->GetPlate();
}

int main (int argc, char *argv[]) {

  logging::core::get()->set_filter
    (
     // logging::trivial::severity >= logging::trivial::info
     logging::trivial::severity >= logging::trivial::debug
     );

  if ( argc != 3 ) {
    BOOST_LOG_TRIVIAL(error) << "Usage : " << argv[0]
			     << " <b2filename> <ofilename>";
    std::exit(1);
  }

  std::string b2filename = argv[1];
  B2Reader reader(b2filename);

  std::string ofilename = argv[2];
  TFile *ofile = new TFile(ofilename.c_str(), "recreate");
  TH1D *hist_ax = new TH1D("hist_ax", "ECC rotation x;tan#theta_x;Entries",
			   40, -1, 1);
  TH1D *hist_ay = new TH1D("hist_ay", "ECC rotation y;tan#theta_y;Entries",
			   40, -1, 1);
  TH2D *hist_pos = new TH2D("hist_pos", "Rotation determinating vertices;x [mm];y [mm]",
			    1000, -3000, 11000, 1000, -2000, 12000);
  TH2D *hist_yz = new TH2D("hist_yz", "Rotation determinating vertices;z [mm];y [mm]",
			    1000, -15000, -5000, 1000, -2000, 12000);
  TH2D *hist_xz = new TH2D("hist_xz", "Rotation determinating vertices;z [mm];x [mm]",
			    1000, -15000, -5000, 1000, -3000, 11000);


  // TSpectrum *s = new TSpectrum();

  // ECC rotation to beam direction correction factor
  
  while ( reader.ReadNextSpill() > 0 ) {
    auto &spill_summary = reader.GetSpillSummary();
    auto it_event = spill_summary.BeginTrueEvent();
    const auto *event = it_event.Next();

    auto &vertex = event->GetPrimaryVertex();
    auto vertex_position = vertex.GetAbsolutePosition().GetValue();
    double weight = vertex.GetTotalCrossSection() * event->GetNormalization();

    // true muon track を見つける
    int true_muon_track_id = -1;
    auto it_outgoing_track = vertex.BeginTrack();
    while ( const auto *track = it_outgoing_track.Next() ) {
      if ( track->GetParticlePdg() == 13 ) {
	true_muon_track_id = track->GetTrackId();
	break;
      }
    }

    if ( true_muon_track_id < 0 ) continue;

    // true muon track に対応する ECC5 内の basetrack をすべて集める
    std::vector<const B2EmulsionSummary* > emulsions;
    auto it_emulsion = spill_summary.BeginEmulsion();
    while ( const auto *emulsion = it_emulsion.Next() ) {
      if ( emulsion->GetFilmType() != B2EmulsionType::kECC ) continue;
      if ( emulsion->GetEcc() != 4 ) continue;
      if ( emulsion->GetParentTrackId() != true_muon_track_id ) continue;
      emulsions.push_back(emulsion);
    }

    if ( emulsions.empty() ) continue;
    BOOST_LOG_TRIVIAL(debug) << "Entry :" << reader.GetEntryNumber();

    // std::sort(emulsions.begin(), emulsions.end(), EmulsionCompare);

    // muon chain が ECC penetrate であることを確認
    bool pene_flag = true;
    BOOST_LOG_TRIVIAL(debug) << "Emulsion penetrates from " << emulsions.front()->GetPlate()
			     << " to " << emulsions.back()->GetPlate();
    if ( emulsions.front()->GetPlate() < 130 ) pene_flag = false;
    if ( emulsions.back()->GetPlate() > 3 ) pene_flag = false;

    if ( !pene_flag ) continue;

    // Baby MIND で muon が reconstruct されている
    // shifter/tracker のマッチングや検出は sand muon に関しては角度依存性が
    // 小さいと仮定する
    if ( spill_summary.GetNumReconTracks() < 1 ) continue;
    bool bm_detect_flag = false;

    auto it_recon_vertex = spill_summary.BeginReconVertex();
    while ( const auto *vertex = it_recon_vertex.Next() ) {
      auto it_recon_track = vertex->BeginTrack();
      while ( const auto *track = it_recon_track.Next() ) {
	if ( track->HasDetector(B2Detector::kBabyMind) ) {
	  bm_detect_flag = true;
	  break;
	}
      }
    }

    if ( !bm_detect_flag ) continue;

    // plate3 での角度分布を作成
    // plate に垂直な方向を基準にスメア
    for ( auto emulsion : emulsions ) {
      if ( emulsion->GetPlate() == 3 ) {
	auto tangent = emulsion->GetTangent().GetValue();
	SmearTangentVector(tangent);
	std::cout << tangent.X() << ", " << tangent.Y() << std::endl;
	hist_ax->Fill(tangent.X(), weight);
	hist_ay->Fill(tangent.Y(), weight);
	hist_pos->Fill(vertex_position.X(), vertex_position.Y(), weight);
	hist_yz->Fill(vertex_position.Z(), vertex_position.Y(), weight);
	hist_xz->Fill(vertex_position.Z(), vertex_position.X(), weight);
      }
    }

  }

  ofile->cd();
  hist_ax->Write();
  hist_ay->Write();
  hist_pos->Write();
  hist_yz->Write();
  hist_xz->Write();
  ofile->Close();

  std::exit(0);

}
