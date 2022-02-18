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
#include <B2Pdg.hh>

// root includes
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>

#include "HistogramStyle.hpp"

namespace logging = boost::log;
namespace fs = boost::filesystem;

int main (int argc, char *argv[]) {

  logging::core::get()->set_filter
    (
     //logging::trivial::severity >= logging::trivial::debug
     logging::trivial::severity >= logging::trivial::info
     );

  if (argc != 3) {
    BOOST_LOG_TRIVIAL(error) << "Usage : " << argv[0]
			     << " <input file name> <output file name>";
    std::exit(1);
  }

  BOOST_LOG_TRIVIAL(info) << "==========Start==========";

  const fs::path input_file(argv[1]);

  try {

    B2Reader reader(input_file);
    
    TFile *ofile = new TFile(argv[2], "recreate");
    TH1D *hist_pion_mom = new TH1D("hist_pion_mom", "Pion momentum;p_{p} [GeV/c];Entries/0.1 GeV/c", 40, 0, 4);
    TH1D *hist_pion_cos = new TH1D("hist_pion_cos", "Pion angle;cos#theta_{p};Entries/0.01", 200, -1, 1); // cos theta distribution
    TH1D *hist_pion_deg = new TH1D("hist_pion_deg", "Pion angle;#theta_{p} [degree];Entries/2 deg", 90, 0, 180);// angle [degree] distribution
    TH2D *hist_pion_mom_cos = new TH2D("hist_pion_mom_cos", "Pion 2d kinematics;p_{p} [GeV/c];cos#theta_{p}", 40, 0, 4, 200, -1, 1);
    TH2D *hist_pion_mom_deg = new TH2D("hist_pion_mom_deg", "Pion 2d kinematics;p_{p} [GeV/c];#theta_{p} [deg]",40, 0, 4, 90, 0, 180);

    while (reader.ReadNextSpill() > 0) {
      auto &input_spill_summary = reader.GetSpillSummary();

      auto it_event = input_spill_summary.BeginTrueEvent();
      const auto *event = it_event.Next(); // get true event summary
      const Double_t norm = event->GetNormalization();

      auto &primary_vertex_summary = event->GetPrimaryVertex();

      const Int_t num_outgoing_tracks = primary_vertex_summary.GetNumOutgoingTracks();
      BOOST_LOG_TRIVIAL(debug) << "Number of outgoing tracks : " << num_outgoing_tracks;
      if (num_outgoing_tracks < 1) continue;

      const Double_t total_cross_section = primary_vertex_summary.GetTotalCrossSection();
      const Double_t weight = norm * total_cross_section * 1.e-38 * 6.02e23;
      const Int_t mode = GetNinjaModeId(primary_vertex_summary.GetInteractionType());
      if ( mode < 0 ) continue; // interaction type not set
      if ( mode == 5 ) continue; // NC

      auto it_track = primary_vertex_summary.BeginTrack();
      while (const auto *track = it_track.Next()) {

	const Int_t particle_id = track->GetParticlePdg();
	if (!B2Pdg::IsChargedPion(particle_id)) continue;
	BOOST_LOG_TRIVIAL(debug) << "Found pion!";
	Double_t pion_mom = track->GetInitialAbsoluteMomentum().GetValue() / 1.e3; // GeV/c
	Double_t pion_deg = track->GetAngle().GetValue(); // degree
	Double_t pion_cos = TMath::Cos(pion_deg * TMath::DegToRad());
	BOOST_LOG_TRIVIAL(debug) << "Pion momentum : " << pion_mom << " GeV/c"
				 << " Pion angle : " << pion_deg << " degree";

	hist_pion_mom->Fill(pion_mom, weight);
	hist_pion_cos->Fill(pion_cos, weight);
	hist_pion_deg->Fill(pion_deg, weight);
	hist_pion_mom_cos->Fill(pion_mom, pion_cos, weight);
	hist_pion_mom_deg->Fill(pion_mom, pion_deg, weight);

      }

    }

    ofile->cd();

    hist_pion_mom->Write();
    hist_pion_cos->Write();
    hist_pion_deg->Write();
    hist_pion_mom_cos->Write();
    hist_pion_mom_deg->Write();

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
