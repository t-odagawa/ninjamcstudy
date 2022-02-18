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
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>

#include "HistogramStyle.hpp"

namespace logging = boost::log;
namespace fs = boost::filesystem;

int main (int argc, char *argv[]) {

  logging::core::get()->set_filter
    (
     logging::trivial::severity >= logging::trivial::info
     );

  if ( argc != 3 ) {
    BOOST_LOG_TRIVIAL(error) << "Usage : " << argv[0]
			     << " <input file neme> <outputfile name>";
    std::exit(1);
  }

  BOOST_LOG_TRIVIAL(info) << "==========Start==========";

  const fs::path input_file(argv[1]);

  try {

    B2Reader reader(input_file);
    
    TFile *ofile = new TFile(argv[2], "recreate");
    TH1D *hist_proton_mom_mode[num_ninja_mode];
    TH1D *hist_proton_cos_mode[num_ninja_mode];
    TH1D *hist_proton_deg_mode[num_ninja_mode];
    TH2D *hist_proton_mom_cos_mode[num_ninja_mode];
    TH2D *hist_proton_mom_deg_mode[num_ninja_mode];

    for ( Int_t imode = 0; imode < num_ninja_mode; imode++ ) {
      hist_proton_mom_mode[imode] = new TH1D(Form("hist_mom_%d", imode),
				      "Proton momentum ("+ mode_name[imode] + ");p_{#mu} [GeV/c];Entries/0.1 GeV/c", 40, 0, 4);
      hist_proton_cos_mode[imode] = new TH1D(Form("hist_cos_%d", imode),
				      "Proton angle ("+ mode_name[imode] + ");cos#theta_{#mu};Entries/0.01", 200, -1, 1);
      hist_proton_deg_mode[imode] = new TH1D(Form("hist_deg_%d", imode),
				      "Proton angle ("+ mode_name[imode] + ");#theta_{#mu} [degree];Entries/2 deg", 90, 0, 180);
      hist_proton_mom_cos_mode[imode] = new TH2D(Form("hist_mom_cos_%d", imode),
					       "Proton 2d kinematics (" + mode_name[imode] + ");p_{#mu} [GeV/c];cos#theta_{#mu}", 40, 0, 4, 200, -1, 1);
      hist_proton_mom_deg_mode[imode] = new TH2D(Form("hist_mom_deg_%d", imode),
					       "Proton 2d kinematics (" + mode_name[imode] + ");p_{#mu} [GeV/c];#theta_{#mu} [deg]", 40, 0, 4, 90, 0, 180);
    }

    while ( reader.ReadNextSpill() > 0 ) {
      auto &input_spill_summary = reader.GetSpillSummary();

      auto it_event = input_spill_summary.BeginTrueEvent();
      const auto *event = it_event.Next(); // get true event summary
      const Double_t norm = event->GetNormalization();

      auto &primary_vertex_summary = event->GetPrimaryVertex();
      
      const Int_t num_outgoing_tracks = primary_vertex_summary.GetNumOutgoingTracks();
      if ( num_outgoing_tracks < 1 ) continue;

      const Double_t total_cross_section = primary_vertex_summary.GetTotalCrossSection();
      const Double_t weight = norm * total_cross_section * 1.e-38 * 6.02e23;
      const Int_t mode = GetNinjaModeId(primary_vertex_summary.GetInteractionType());
      if ( mode < 0 ) continue; // interaction type not set
      if ( mode == 5 ) continue; // NC

      auto it_track = primary_vertex_summary.BeginTrack();
      while ( const auto *track = it_track.Next() ) {
	if ( track->GetParticlePdg() != PDG_t::kProton ) continue;
	Double_t proton_mom = track->GetInitialAbsoluteMomentum().GetValue() / 1.e3; // GeV/c
	Double_t proton_deg = track->GetAngle().GetValue(); // degree
	Double_t proton_cos = TMath::Cos(proton_deg * TMath::DegToRad());
	BOOST_LOG_TRIVIAL(debug) << "Proton momentum : " << proton_mom << " GeV/c, "
				 << "Proton angle : " << proton_deg << "degree";

	hist_proton_mom_mode[mode]->Fill(proton_mom, weight);
	hist_proton_cos_mode[mode]->Fill(proton_cos, weight);
	hist_proton_deg_mode[mode]->Fill(proton_deg, weight);
	hist_proton_mom_cos_mode[mode]->Fill(proton_mom, proton_cos, weight);
	hist_proton_mom_deg_mode[mode]->Fill(proton_mom, proton_deg, weight);
      }

    }

    ofile->cd();
    for ( Int_t imode = 0; imode < num_ninja_mode; imode++ ) {
      hist_proton_mom_mode[imode]->Write();
      hist_proton_cos_mode[imode]->Write();
      hist_proton_deg_mode[imode]->Write();
      hist_proton_mom_cos_mode[imode]->Write();
      hist_proton_mom_deg_mode[imode]->Write();
    }
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
