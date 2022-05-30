#include <iostream>
#include <vector>
#include <string>

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/filesystem.hpp>

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>

#include <B2Enum.hh>
#include <B2Reader.hh>
#include <B2SpillSummary.hh>
#include <B2EventSummary.hh>
#include <B2VertexSummary.hh>

#include <McsClass.hpp>

#include "HistogramStyle.hpp"
#include "AnalyzeMomentum.hpp"

namespace logging = boost::log;
namespace fs = boost::filesystem;

void AnalyzeMomentum(std::string b2filename,
		     std::string momchfilename,
		     std::string outputfilename) {
  
  BOOST_LOG_TRIVIAL(info) << "==========Momentum mode==========";

  // input B2 file
  B2Reader reader(b2filename);
  
  // input momch file
  if ( !fs::exists(momchfilename) ) {
    throw std::runtime_error("File not found : " + momchfilename);
  }
  std::vector<Momentum_recon::Event_information > ev_vec = Momentum_recon::ReadEventInformationBin(momchfilename);

  // output file
  TFile *outputfile = new TFile((TString)outputfilename, "recreate");
  BOOST_LOG_TRIVIAL(info) << "Output filename : " << outputfilename;

  TH1D *hist_muon_mom = new TH1D("hist_muon_mom", "Muon reconstructed momentum;p_{#mu};Entries", 50, 0., 1500.);
  TH2D *hist_muon_mom_recon_true = new TH2D("hist_muon_mom_recon_true", "Muon momentum;p_{#mu, true};p_{#mu, recon}", 15, 0., 1500., 15, 0., 1500.);

  TH1D *hist_mode_muon_mom[num_ninja_mode];
  for ( int i = 0; i < num_ninja_mode; i++ ) {
    hist_mode_muon_mom[i] = new TH1D(Form("hist_muon_mom_%d", i), "", 15, 0., 1500.);
    hist_mode_muon_mom[i]->SetFillColor(mode_color[i]);
    hist_mode_muon_mom[i]->SetFillStyle(mode_style[i]);
  }

  for ( auto ev : ev_vec ) {

    reader.ReadSpill(ev.groupid);
    auto &spill_summary = reader.GetSpillSummary();
    auto it_event = spill_summary.BeginTrueEvent();
    const auto *event = it_event.Next();

    auto &vertex = event->GetPrimaryVertex();
    int mode_id = GetNinjaModeId(vertex.GetInteractionType());

    if ( ev.vertex_material != B2Material::kWater ) continue;

    if ( !ev.chains.empty() ) {
      for ( auto chain : ev.chains ) {
	
	// search corresponding true chain
	double true_momentum = -1;
	for ( auto true_chain : ev.true_chains ) {
	  if ( chain.chainid == true_chain.chainid )
	    true_momentum = true_chain.bm_range_mom;
	}

	int particle_id = chain.particle_flag % 10000;
	double recon_momentum = -1;
	if ( particle_id == 13 ) {
	  recon_momentum = chain.ecc_mcs_mom[0];
	  hist_muon_mom->Fill(recon_momentum, ev.weight);
	  hist_muon_mom_recon_true->Fill(true_momentum, recon_momentum, ev.weight);
	  hist_mode_muon_mom[mode_id]->Fill(recon_momentum, ev.weight);
	}
	else if ( particle_id == 211 ) {
	  recon_momentum = chain.bm_curvature_mom;
	}
	else if ( particle_id == 2212 ) {
	  recon_momentum = chain.ecc_mcs_mom[1];
	}

	
      }
    }

  }

  outputfile->cd();
  hist_muon_mom->Write();
  hist_muon_mom_recon_true->Write();
  for ( int i = 0; i < num_ninja_mode; i++ ) {
    hist_mode_muon_mom[i]->Write();
  }
  outputfile->Close();

}
