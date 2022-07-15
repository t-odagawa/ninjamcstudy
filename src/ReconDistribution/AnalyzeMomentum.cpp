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
#include "DrawConst.hpp"
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

  TH1D *hist_muon_mom = new TH1D("hist_muon_mom",
				 "Muon reconstructed momentum;p_{#mu} [MeV/c];Entries",
				 20, 0., 2000.);
  TH1D *hist_muon_mom_mcs = new TH1D("hist_muon_mom_mcs",
				     "Muon MCS reconstructed momentum;p_{#mu, MCS} [MeV/c];Entries",
				     20, 0., 2000.);
  TH1D *hist_muon_mom_range = new TH1D("hist_muon_mom_range",
				       "Muon range reconstructed momentum;p_{#mu, range};Entries",
				       20, 0., 2000.);
  TH2D *hist_muon_mom_recon_true = new TH2D("hist_muon_mom_recon_true",
					    "Muon momentum;p_{#mu, true} [MeV/c];p_{#mu, recon} [MeV/c]",
					    50, 0., 1500., 50, 0., 1500.);
  TH2D *hist_muon_mom_recon_true_mcs = new TH2D("hist_muon_mom_recon_true_mcs",
						"Muon MCS momentum;p_{#mu, true} [MeV/c];p_{#mu, recon} [MeV/c]",
						50, 0., 1500., 50, 0., 1500.);
  TH2D *hist_muon_mom_recon_true_range = new TH2D("hist_muon_mom_recon_true_range",
						  "Muon range momentum;p_{#mu, true} [MeV/c];p_{#mu, recon} [MeV/c]",
						  50, 0., 1500., 50, 0., 1500.);
  TH1D *hist_pion_mom = new TH1D("hist_pion_mom",
				 "Pion reconstructed momentum;p_{#pi} [MeV/c];Entries",
				 50, 0., 1500.);
  TH1D *hist_pion_mom_mcs = new TH1D("hist_pion_mom_mcs",
				     "Pion MCS reconstructed momentum;p_{#pi, MCS} [MeV/c];Entries",
				     50, 0., 1500.);
  TH1D *hist_pion_mom_range = new TH1D("hist_pion_mom_range",
				       "Pion range reconstructed momentum;p_{#pi, range} [MeV/c];Entries",
				       50, 0., 1500.);
  TH1D *hist_proton_mom = new TH1D("hist_proton_mom",
				   "Proton reconstructed momentum;p_{p} [MeV/c];Entries",
				   50, 0., 1500.);
  TH1D *hist_proton_mom_mcs = new TH1D("hist_proton_mom_mcs",
				       "Proton MCS reconstructed momentum;p_{p, MCS} [MeV/c];Entries",
				       50, 0., 1500.);
  TH1D *hist_proton_mom_range = new TH1D("hist_proton_mom_range",
					 "Proton range reconstructed momentum;p_{p, range} [MeV/c];Entries",
					 50, 0., 1500.);
  TH2D *hist_pion_mom_recon_true = new TH2D("hist_pion_mom_recon_true",
					    "Pion momentum;p_{#pi, true} [MeV/c];p_{#pi, recon} [MeV/c]",
					    50, 0., 1500., 50, 0., 1500.);
  TH2D *hist_pion_mom_recon_true_mcs = new TH2D("hist_pion_mom_recon_true_mcs",
						"Pion MCS momentum;p_{#pi, true} [MeV/c];p_{#pi, recon} [MeV/c]",
						50, 0., 1500., 50, 0., 1500.);
  TH2D *hist_pion_mom_recon_true_range = new TH2D("hist_pion_mom_recon_true_range",
						  "Pion range momentum;p_{#pi, true} [MeV/c];p_{#pi, recon} [MeV/c]",
						  50, 0., 1500., 50, 0., 1500.);
  TH2D *hist_proton_mom_recon_true = new TH2D("hist_proton_mom_recon_true",
					      "Proton momentum;p_{p, true} [MeV/c];p_{p, recon} [MeV/c]",
					      50, 0., 1500., 50, 0., 1500.);
  TH2D *hist_proton_mom_recon_true_mcs = new TH2D("hist_proton_mom_recon_true_mcs",
						  "Proton MCS momentum;p_{p, true} [MeV/c];p_{p, recon} [MeV/c]",
						  50, 0., 1500., 50, 0., 1500.);
  TH2D *hist_proton_mom_recon_true_range = new TH2D("hist_proton_mom_recon_true_range",
						    "Proton range momentum;p_{p, true} [MeV/c];p_{p, recon} [MeV/c]",
						    50, 0., 1500., 50, 0., 1500.);
  TH2D *hist_pbeta_muon_proton = new TH2D("hist_pbeta_muon_proton",
					  "Reconstructed pbeta difference;p#beta_{#mu} [MeV/c];p#beta_{p} [MeV/c]",
					  50, 0., 1500., 50, 0., 1500);
  TH1D *hist_muon_misid_mom = new TH1D("hist_muon_misid_mom",
				       "\"Muon\" reconstructed momentum;p_{#mu} [MeV/c];Entries",
				       20, 0., 2000.);
  TH1D *hist_muon_misid_mom_mcs = new TH1D("hist_muon_misid_mom_mcs",
					   "\"Muon\" MCS reconstructed momentum;p_{#mu} [MeV/c];Entries",
					   20, 0., 2000.);
  TH1D *hist_muon_misid_mom_range = new TH1D("hist_muon_misid_mom_range",
					     "\"Muon\" range reconstructed momentum;p_{#mu} [MeV/c];Entries",
					     20, 0., 2000.);
  TH1D *hist_proton_misid_mom = new TH1D("hist_proton_misid_mom",
					 "\"Proton\" reconstructed momentum;p_{p} [MeV/c];Entries",
					 50, 0., 1500.);
  TH1D *hist_proton_misid_mom_mcs = new TH1D("hist_proton_misid_mom_mcs",
					     "\"Proton\" MCS reconstructed momentum;p_{p} [MeV/c];Entries",
					     50, 0., 1500.);
  TH1D *hist_proton_misid_mom_range = new TH1D("hist_proton_misid_mom_range",
					       "\"Proton\" range reconstructed momentum;p_{p} [MeV/c];Entries",
					       50, 0., 1500.);
  TH1D *hist_pion_misid_mom = new TH1D("hist_pion_misid_mom",
				       "\"Pion\" reconstructed momentum;p_{#pi} [MeV/c];Entries",
				       50, 0., 1500.);
  TH1D *hist_pion_misid_mom_mcs = new TH1D("hist_pion_misid_mom_mcs",
					   "\"Pion\" MCS reconstructed momentum;p_{#pi} [MeV/c];Entries",
					   50, 0., 1500.);
  TH1D *hist_pion_misid_mom_range = new TH1D("hist_pion_misid_mom_range",
					      "\"Pion\" range reconstructed momentum;p_{#pi} [MeV/c];Entries",
					      50, 0., 1500.);

  TH1D *hist_muon_mom_single = new TH1D("hist_muon_mom_single", "Muon reconstructed momentum;p_{#mu} [MeV/c];Entries",
					20, 0., 2000.);
  TH1D *hist_muon_mom_single_mcs = new TH1D("hist_muon_mom_single_mcs", "Muon reconstructed momentum;p_{#mu, MCS} [MeV/c];Entries",
					    20, 0., 2000.);
  TH1D *hist_muon_mom_single_range = new TH1D("hist_muon_mom_single_range", "Muon reconstructed momentum;p_{#mu, range} [MeV/c];Entries",
					      20, 0., 2000.);

  TH1D *hist_mode_muon_mom[num_ninja_mode];
  TH1D *hist_mode_pion_mom[num_ninja_mode];
  TH1D *hist_mode_proton_mom[num_ninja_mode];
  for ( int i = 0; i < num_ninja_mode; i++ ) {
    hist_mode_muon_mom[i] = new TH1D(Form("hist_muon_mom_%d", i), "", 20, 0., 2000.);
    hist_mode_muon_mom[i]->SetFillColor(mode_color[i]);
    hist_mode_muon_mom[i]->SetFillStyle(mode_style[i]);
    hist_mode_pion_mom[i] = new TH1D(Form("hist_pion_mom_%d", i), "", 50, 0., 1500.);
    hist_mode_pion_mom[i]->SetFillColor(mode_color[i]);
    hist_mode_pion_mom[i]->SetFillStyle(mode_style[i]);
    hist_mode_proton_mom[i] = new TH1D(Form("hist_proton_mom_%d", i), "", 50, 0., 1500.);
    hist_mode_proton_mom[i]->SetFillColor(mode_color[i]);
    hist_mode_proton_mom[i]->SetFillStyle(mode_style[i]);
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
	  if ( chain.chainid == true_chain.chainid ) {
	    true_momentum = true_chain.bm_range_mom;
	    break;
	  }
	}

	int particle_id = chain.particle_flag % 10000;
	int true_particle_id = chain.particle_flag / 10000;

	double recon_momentum = -1;
	if ( particle_id == true_particle_id ) {
	  if ( particle_id == 13 ) {
	    if ( chain.stop_flag == 1 ) {
	      recon_momentum = chain.bm_range_mom;
	      hist_muon_mom_range->Fill(recon_momentum, ev.weight);
	      hist_muon_mom_recon_true_range->Fill(true_momentum, recon_momentum, ev.weight);
	      if ( ev.chains.size() == 1 ) {
		hist_muon_mom_single_range->Fill(recon_momentum, ev.weight);
	      }
	    }
	    else if (chain.stop_flag == 0 ) {
	      recon_momentum = chain.ecc_mcs_mom[0];
	      hist_muon_mom_mcs->Fill(recon_momentum, ev.weight);
	      hist_muon_mom_recon_true_mcs->Fill(true_momentum, recon_momentum, ev.weight);
	      if ( ev.chains.size() == 1 ) {
		hist_muon_mom_single_mcs->Fill(recon_momentum, ev.weight);
	      }
	    }
	    hist_muon_mom->Fill(recon_momentum, ev.weight);
	    hist_muon_mom_recon_true->Fill(true_momentum, recon_momentum, ev.weight);
	    hist_mode_muon_mom[mode_id]->Fill(recon_momentum, ev.weight);
	    if ( ev.chains.size() == 1 ) {
	      hist_muon_mom_single->Fill(recon_momentum, ev.weight);
	    }
	  }
	  else if ( particle_id == 211 ) {
	    if ( chain.stop_flag == 0 ) {
	      recon_momentum = chain.bm_curvature_mom;
	      hist_pion_mom_mcs->Fill(recon_momentum, ev.weight);
	      hist_pion_mom_recon_true_mcs->Fill(true_momentum, recon_momentum, ev.weight);
	    }
	    else if ( chain.stop_flag == 2 ) {
	      recon_momentum = chain.ecc_range_mom[0];
	      hist_pion_mom_range->Fill(recon_momentum, ev.weight);
	      hist_pion_mom_recon_true_range->Fill(true_momentum, recon_momentum, ev.weight);
	    }
	    hist_pion_mom->Fill(recon_momentum, ev.weight);
	    hist_pion_mom_recon_true->Fill(true_momentum, recon_momentum, ev.weight);
	    hist_mode_pion_mom[mode_id]->Fill(recon_momentum, ev.weight);
	  }
	  else if ( particle_id == 2212 ) {
	    double pbeta_mu = chain.ecc_mcs_mom[0];
	    pbeta_mu = pbeta_mu * pbeta_mu / std::sqrt(pbeta_mu * pbeta_mu + muon_mass * muon_mass);
	    double pbeta_p = chain.ecc_mcs_mom[1];
	    pbeta_p = pbeta_p * pbeta_p / std::sqrt(pbeta_p * pbeta_p + proton_mass * proton_mass);
	    hist_pbeta_muon_proton->Fill(pbeta_mu, pbeta_p, ev.weight);
	    if ( chain.stop_flag == 0 ) {
	      recon_momentum = chain.ecc_mcs_mom[1];
	      hist_proton_mom_mcs->Fill(recon_momentum, ev.weight);
	      hist_proton_mom_recon_true_mcs->Fill(true_momentum, recon_momentum, ev.weight);
	    }
	    else if ( chain.stop_flag == 2 ) {
	      recon_momentum = chain.ecc_range_mom[1];
	      hist_proton_mom_range->Fill(recon_momentum, ev.weight);
	      hist_proton_mom_recon_true_range->Fill(true_momentum, recon_momentum, ev.weight);
	    }
	    hist_proton_mom->Fill(recon_momentum, ev.weight);
	    hist_proton_mom_recon_true->Fill(true_momentum, recon_momentum, ev.weight);
	    hist_mode_proton_mom[mode_id]->Fill(recon_momentum, ev.weight);
	  }
	}
	else if ( particle_id != true_particle_id ) {
	  if ( particle_id == 13 ) {
	    if ( chain.stop_flag == 0 ) {
	      recon_momentum = chain.ecc_mcs_mom[0];
	      hist_muon_misid_mom_mcs->Fill(recon_momentum, ev.weight);
	    }
	    else if ( chain.stop_flag == 1 ) {
	      recon_momentum = chain.bm_range_mom;
	      hist_muon_misid_mom_range->Fill(recon_momentum, ev.weight);
	    }
	    hist_muon_misid_mom->Fill(recon_momentum, ev.weight);
	  }
	  else if ( particle_id == 211 ) {
	    if ( chain.stop_flag == 0 ) {
	      recon_momentum = chain.bm_curvature_mom;
	      hist_pion_misid_mom_mcs->Fill(recon_momentum, ev.weight);
	    }
	    else if ( chain.stop_flag == 2 ) {
	      recon_momentum = chain.ecc_range_mom[0];
	      hist_pion_misid_mom_range->Fill(recon_momentum, ev.weight);
	    }
	    hist_pion_misid_mom->Fill(recon_momentum, ev.weight);
	  }
	  else if ( particle_id == 2212 ) {
	    if ( chain.stop_flag == 0 ) {
	      recon_momentum = chain.ecc_mcs_mom[1];
	      hist_proton_misid_mom_mcs->Fill(recon_momentum, ev.weight);	      
	    }
	    else if ( chain.stop_flag == 2 ) {
	      recon_momentum = chain.ecc_range_mom[1];
	      hist_proton_misid_mom_range->Fill(recon_momentum, ev.weight);
	    }
	    hist_proton_misid_mom->Fill(recon_momentum, ev.weight);
	  }
	}
      }
    }

  }

  outputfile->cd();
  hist_muon_mom->Write();
  hist_muon_mom_mcs->Write();
  hist_muon_mom_range->Write();
  hist_muon_mom_recon_true->Write();
  hist_muon_mom_recon_true_mcs->Write();
  hist_muon_mom_recon_true_range->Write();
  hist_pion_mom->Write();
  hist_pion_mom_mcs->Write();
  hist_pion_mom_range->Write();
  hist_pion_mom_recon_true->Write();
  hist_pion_mom_recon_true_mcs->Write();
  hist_pion_mom_recon_true_range->Write();
  hist_proton_mom->Write();
  hist_proton_mom_mcs->Write();
  hist_proton_mom_range->Write();
  hist_proton_mom_recon_true->Write();
  hist_proton_mom_recon_true_mcs->Write();
  hist_proton_mom_recon_true_range->Write();
  hist_pbeta_muon_proton->Write();

  hist_muon_misid_mom->Write();
  hist_muon_misid_mom_mcs->Write();
  hist_muon_misid_mom_range->Write();
  hist_pion_misid_mom->Write();
  hist_pion_misid_mom_mcs->Write();
  hist_pion_misid_mom_range->Write();
  hist_proton_misid_mom->Write();
  hist_proton_misid_mom_mcs->Write();
  hist_proton_misid_mom_range->Write();

  hist_muon_mom_single->Write();
  hist_muon_mom_single_mcs->Write();
  hist_muon_mom_single_range->Write();

  for ( int i = 0; i < num_ninja_mode; i++ ) {
    hist_mode_muon_mom[i]->Write();
    hist_mode_pion_mom[i]->Write();
    hist_mode_proton_mom[i]->Write();
  }
  outputfile->Close();

}
