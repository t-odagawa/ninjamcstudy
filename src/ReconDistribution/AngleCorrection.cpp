#include <string>
#include <vector>

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

#include <TFile.h>
#include <TH1D.h>

#include <B2Reader.hh>
#include <B2SpillSummary.hh>
#include <B2EventSummary.hh>
#include <B2TrackSummary.hh>
#include <B2VertexSummary.hh>
#include <B2EmulsionSummary.hh>
#include <B2Enum.hh>

namespace logging = boost::log;

bool EmulsionCompare(const B2EmulsionSummary* lhs,
		       const B2EmulsionSummary* rhs) {
return lhs->GetPlate() > rhs->GetPlate();
}

int main (int argc, char *argv[]) {

  logging::core::get()->set_filter
    (
     logging::trivial::severity >= logging::trivial::info
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
			   100, -1, 1);
  TH1D *hist_ay = new TH1D("hist_ay", "ECC rotation y;tan#theta_y;Entries",
			   100, -1, 1);

  // TSpectrum *s = new TSpectrum();

  // ECC rotation to beam direction correction factor
  
  while ( reader.ReadNextSpill() > 0 ) {
    auto &spill_summary = reader.GetSpillSummary();
    auto it_event = spill_summary.BeginTrueEvent();
    const auto *event = it_event.Next();

    auto &vertex = event->GetPrimaryVertex();
    double weight = vertex.GetMcWeight();

    // true muon track を見つける
    int true_muon_track_id = -1;
    auto it_true_track = spill_summary.BeginTrueTrack();    
    while ( const auto *track = it_true_track.Next() ) {
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

    // std::sort(emulsions.begin(), emulsions.end(), EmulsionCompare);

    // muon chain が ECC penetrate であることを確認
    bool pene_flag = true;
    if ( emulsions.front()->GetPlate() < 130 ) pene_flag = false;
    if ( emulsions.back()->GetPlate() > 3 ) pene_flag = false;

    if ( !pene_flag ) continue;

    // plate3 での角度分布を作成
    // plate に垂直な方向を基準にスメア
    for ( auto emulsion : emulsions ) {
      if ( emulsion->GetPlate() == 3 ) {
	hist_ax->Fill(emulsion->GetTangent().GetValue().X(), weight);
	hist_ay->Fill(emulsion->GetTangent().GetValue().Y(), weight);
      }
    }

  }

  ofile->cd();
  hist_ax->Write();
  hist_ay->Write();
  ofile->Close();

  std::exit(0);

}
