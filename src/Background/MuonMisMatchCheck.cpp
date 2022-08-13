#include <iostream>
#include <vector>
#include <cmath>
#include <string>

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/filesystem.hpp>

#include <NTBMSummary.hh>
#include <NTBMConst.hh>

#include <B2Const.hh>
#include <B2Enum.hh>

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>

#include <McsClass.hpp>

namespace logging = boost::log;
namespace fs = boost::filesystem;

int main(int argc, char* argv[]) {

  if ( argc != 4 ) {
    BOOST_LOG_TRIVIAL(error) << "Usage : " << argv[0]
			     << " <Input NTBM file> <Input momch file> <Output root file>";
    std::exit(1);
  }

  TFile *ifile = new TFile(argv[1], "read");
  TTree *tree = (TTree*)ifile->Get("tree");
  NTBMSummary *ntbm = nullptr;
  tree->SetBranchAddress("NTBMSummary", &ntbm);

  const double bm_nt_distance = (BABYMIND_POS_Z + BM_SECOND_LAYER_POS)
    - (NINJA_POS_Z + NINJA_ECC_POS_Z + NINJA_DESIC_DEPTH / 2. - NINJA_DESIC_THICK
       - NINJA_ENV_THICK - 3. * NINJA_FILM_THICK - NINJA_SS_AC_THICK);

  auto ev_vec = Momentum_recon::ReadEventInformationBin(argv[2]);

  TFile *ofile = new TFile(argv[3], "recreate");
  TH1D *hist_signal = new TH1D("hist_signal","", 100, 0, 1500);
  TH1D *hist_mis = new TH1D("hist_mis", "", 100, 0, 1500);

  for ( auto ev : ev_vec ) {

    if ( ev.chains.empty() ) continue;

    tree->GetEntry(ev.groupid);

    double track_length = -1.;
    std::vector<double > muon_position;
    std::vector<double > muon_tangent;

    for ( int itrack = 0; itrack < ntbm->GetNumberOfTracks(); itrack++ ) {
      if ( ntbm->GetTrackLengthTotal(itrack) > track_length ) {
	track_length = ntbm->GetTrackLengthTotal(itrack);
	muon_position = ntbm->GetBabyMindPosition(itrack);
	muon_position.at(B2View::kSideView) += BABYMIND_POS_Y;
	muon_position.at(B2View::kTopView) += BABYMIND_POS_X;
	muon_tangent = ntbm->GetBabyMindTangent(itrack);
      }
    }

    for ( auto chain : ev.chains ) {

      int recon_particle_id = chain.particle_flag % 10000;
      int true_particle_id = chain.particle_flag / 10000;

      if ( recon_particle_id != 13 ) continue;

      double x = chain.base.front().x / 1000.;
      double y = chain.base.front().y / 1000.;
      x = x - 125. + NINJA_POS_X + NINJA_ECC_POS_X;
      y = y - 125. + NINJA_POS_Y + NINJA_FV_IRON_POS_Y;
      double ax = chain.base.front().ax;
      double ay = chain.base.front().ay;

      double dx = x + ax * bm_nt_distance - muon_position.at(B2View::kTopView);
      double dy = y + ay * bm_nt_distance - muon_position.at(B2View::kSideView);

      double position_difference = std::hypot(dx, dy);

      if ( recon_particle_id == true_particle_id ) {
	hist_signal->Fill(position_difference, ev.weight);
      }
      else {
	hist_mis->Fill(position_difference, ev.weight);
      }

    }

  }

  ofile->cd();
  hist_signal->Write();
  hist_mis->Write();
  ofile->Close();
  
  std::exit(0);

}
