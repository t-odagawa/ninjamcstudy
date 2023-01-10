#include <iostream>
#include <vector>
#include <string>
#include <cmath>

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/filesystem.hpp>

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>

#include <B2Enum.hh>
#include <B2Reader.hh>
#include <B2SpillSummary.hh>
#include <B2EventSummary.hh>
#include <B2VertexSummary.hh>

#include <McsClass.hpp>

#include "HistogramStyle.hpp"
#include "DrawConst.hpp"

namespace logging = boost::log;
namespace fs = boost::filesystem;


void AnalyzeAngle(std::string b2filename,
		  std::string momchfilename,
		  std::string outputfilename) {

  BOOST_LOG_TRIVIAL(info) << "==========Angle mode===========";

  // input B2 file
  B2Reader reader(b2filename);

  // input momch file
  if ( !fs::exists(momchfilename) ) {
    throw std::runtime_error("File not found : " + momchfilename);   
  }
  auto ev_vec = Momentum_recon::ReadEventInformationBin(momchfilename);

  // output file
  TFile *outputfile = new TFile((TString)outputfilename, "recreate");
  BOOST_LOG_TRIVIAL(info) << "Output filename : " << outputfilename;

  TH1D *hist_muon_ang_cos = new TH1D("hist_muon_ang_cos", 
				     "Muon reconstructed angle;cos#theta_{#mu};Entries",
				     muon_cos_bin_size - 1, muon_cos_bins);
  TH1D *hist_muon_ang_deg = new TH1D("hist_muon_ang_deg",
				     "Muon reconstructed angle;#theta_{#mu} [deg];Entries",
				     muon_deg_bin_size - 1, muon_deg_bins);
  TH1D *hist_pion_ang_cos = new TH1D("hist_pion_ang_cos",
				     "Pion reconstructed angle;cos#theta_{#pi};Entries",
				     hadron_cos_bin_size - 1, hadron_cos_bins);
  TH1D *hist_pion_ang_deg = new TH1D("hist_pion_ang_deg",
				     "Pion reconstructed angle;#theta_{#pi} [deg];Entries",
				     hadron_deg_bin_size - 1, hadron_deg_bins);
  TH1D *hist_proton_ang_cos = new TH1D("hist_proton_ang_cos",
				       "Proton reconstructed angle;cos#theta_{p};Entries",
				       hadron_cos_bin_size - 1, hadron_cos_bins);
  TH1D *hist_proton_ang_deg = new TH1D("hist_proton_ang_deg",
				       "Proton reconstructed angle;#theta_{p} [deg];Entries",
				       hadron_deg_bin_size - 1, hadron_deg_bins);
  TH1D *hist_other_ang_cos = new TH1D("hist_other_ang_cos",
				      "",
				      hadron_cos_bin_size - 1, hadron_cos_bins);
  TH1D *hist_other_ang_deg = new TH1D("hist_other_ang_deg",
				      "",
				      hadron_deg_bin_size - 1, hadron_deg_bins);
  
  TH1D *hist_muon_ang_cos_single = new TH1D("hist_muon_ang_cos_single",
					    "Muon reconstructed angle;cos#theta_{#mu};Entries",
					    muon_cos_bin_size - 1, muon_cos_bins);
  TH1D *hist_muon_ang_deg_single = new TH1D("hist_muon_ang_deg_single",
					    "Muon reconstructed angle;#theta_{#mu} [deg];Entries",
					    muon_deg_bin_size - 1, muon_deg_bins);

  // Muon is correctly id-ed
  TH1D *hist_mode_muon_ang_cos[num_ninja_mode];
  TH1D *hist_mode_muon_ang_deg[num_ninja_mode];
  TH1D *hist_mode_pion_ang_cos[num_ninja_mode];
  TH1D *hist_mode_pion_ang_deg[num_ninja_mode];
  TH1D *hist_mode_proton_ang_cos[num_ninja_mode];
  TH1D *hist_mode_proton_ang_deg[num_ninja_mode];
  for ( int i = 0; i < num_ninja_mode; i++ ) {
    hist_mode_muon_ang_cos[i] = new TH1D(Form("hist_muon_ang_cos_%d", i), "", muon_cos_bin_size - 1, muon_cos_bins);
    hist_mode_muon_ang_cos[i]->SetFillColor(mode_color[i]);
    hist_mode_muon_ang_cos[i]->SetFillStyle(mode_style[i]);
    hist_mode_muon_ang_deg[i] = new TH1D(Form("hist_muon_ang_deg_%d", i), "", muon_deg_bin_size - 1, muon_deg_bins);
    hist_mode_muon_ang_deg[i]->SetFillColor(mode_color[i]);
    hist_mode_muon_ang_deg[i]->SetFillStyle(mode_style[i]);
    hist_mode_pion_ang_cos[i] = new TH1D(Form("hist_pion_ang_cos_%d", i), "", hadron_cos_bin_size - 1, hadron_cos_bins);
    hist_mode_pion_ang_cos[i]->SetFillColor(mode_color[i]);
    hist_mode_pion_ang_cos[i]->SetFillStyle(mode_style[i]);
    hist_mode_pion_ang_deg[i] = new TH1D(Form("hist_pion_ang_deg_%d", i), "", hadron_deg_bin_size - 1, hadron_deg_bins);
    hist_mode_pion_ang_deg[i]->SetFillColor(mode_color[i]);
    hist_mode_pion_ang_deg[i]->SetFillStyle(mode_style[i]);
    hist_mode_proton_ang_cos[i] = new TH1D(Form("hist_proton_ang_cos_%d", i), "", hadron_cos_bin_size - 1, hadron_cos_bins);
    hist_mode_proton_ang_cos[i]->SetFillColor(mode_color[i]);
    hist_mode_proton_ang_cos[i]->SetFillStyle(mode_style[i]);
    hist_mode_proton_ang_deg[i] = new TH1D(Form("hist_proton_ang_deg_%d", i), "", hadron_deg_bin_size - 1, hadron_deg_bins);
    hist_mode_proton_ang_deg[i]->SetFillColor(mode_color[i]);
    hist_mode_proton_ang_deg[i]->SetFillStyle(mode_style[i]);
  }

  // Muon is not id-ed
  TH1D *hist_muon_misid_ang_cos = new TH1D("hist_muon_misid_ang_cos",
					   "\"Muon\" reconstructed angle;cos#theta_{#mu};Entries",
					   muon_cos_bin_size - 1, muon_cos_bins);
  TH1D *hist_muon_misid_ang_deg = new  TH1D("hist_muon_misid_ang_deg",
					    "\"Muon\" reconstructed angle;#theta_{#mu} [deg];Entries",
					    muon_deg_bin_size - 1, muon_deg_bins);
  TH1D *hist_pion_ang_cos_muon_misid = new TH1D("hist_pion_ang_cos_muon_misid",
						"", hadron_cos_bin_size - 1, hadron_cos_bins);
  TH1D *hist_pion_ang_deg_muon_misid = new TH1D("hist_pion_ang_deg_muon_misid",
						"", hadron_deg_bin_size - 1, hadron_deg_bins);
  TH1D *hist_proton_ang_cos_muon_misid = new TH1D("hist_proton_ang_cos_muon_misid",
						  "", hadron_cos_bin_size - 1, hadron_cos_bins);
  TH1D *hist_proton_ang_deg_muon_misid = new TH1D("hist_proton_ang_deg_muon_misid",
						  "", hadron_deg_bin_size - 1, hadron_deg_bins);

  // Muon is correctly id-ed but partner is not
  TH1D *hist_proton_misid_ang_cos = new TH1D("hist_proton_misid_ang_cos",
					     "\"Proton\" reconstructed angle;cos#theta_{p};Entries",
					     hadron_cos_bin_size - 1, hadron_cos_bins);
  TH1D *hist_proton_misid_ang_deg = new TH1D("hist_proton_misid_ang_deg",
					     "\"Proton\" reconstructed angle;#theta_{p} [deg];Entries",
					     hadron_deg_bin_size - 1, hadron_deg_bins);
  TH1D *hist_pion_misid_ang_cos = new TH1D("hist_pion_misid_ang_cos",
					   "\"Pion\" reconstructed angle;cos#theta_{#pi};Entries",
					   hadron_cos_bin_size - 1, hadron_cos_bins);
  TH1D *hist_pion_misid_ang_deg = new TH1D("hist_pion_misid_ang_deg",
					   "\"Pion\" reconstructed angle;#theta_{#pi} [deg];Entries",
					   hadron_deg_bin_size - 1, hadron_deg_bins);

  // For flux systematic uncertainty study
  TH2D *hist_flux_mu_deg = new TH2D("hist_flux_mu_deg", ";#theta_{#mu} [deg];E_{#nu} [GeV]",
				    muon_deg_bin_size - 1, muon_deg_bins, nu_ene_bin_size - 1, nu_ene_bins);
  TH2D *hist_flux_mu_cos = new TH2D("hist_flux_mu_cos", ";cos#theta_{#mu};E_{#nu} [GeV]",
				    muon_cos_bin_size - 1, muon_cos_bins, nu_ene_bin_size - 1, nu_ene_bins);
  TH2D *hist_flux_p_deg = new TH2D("hist_flux_p_deg", ";#theta_{p} [deg];E_{#nu} [GeV]",
				   hadron_deg_bin_size - 1, hadron_deg_bins, nu_ene_bin_size - 1, nu_ene_bins);
  TH2D *hist_flux_p_cos = new TH2D("hist_flux_p_cos", ";cos#theta_{p};E_{#nu} [GeV]",
				   hadron_cos_bin_size - 1, hadron_cos_bins, nu_ene_bin_size - 1, nu_ene_bins);
  TH2D *hist_flux_pi_deg = new TH2D("hist_flux_pi_deg", ";#theta_{#pi} [deg];E_{#nu} [GeV]",
				    hadron_deg_bin_size - 1, hadron_deg_bins, nu_ene_bin_size - 1, nu_ene_bins);
  TH2D *hist_flux_pi_cos = new TH2D("hist_flux_pi_cos", ";cos#theta_{#pi};E_{#nu} [GeV]",
				    hadron_cos_bin_size - 1, hadron_cos_bins, nu_ene_bin_size - 1, nu_ene_bins);

  for ( auto ev : ev_vec ) {

    reader.ReadSpill(ev.groupid);
    auto &spill_summary = reader.GetSpillSummary();
    auto it_event = spill_summary.BeginTrueEvent();
    const auto *event = it_event.Next();

    auto &vertex = event->GetPrimaryVertex();
    int mode_id = GetNinjaModeId(vertex.GetInteractionType());

    if ( ev.vertex_material != B2Material::kWater ) continue;

    if ( !ev.chains.empty() ) {      

      bool muon_correct_id_flag = false;
      int num_proton_water = 0;
      int num_true_proton_water = 0;
      int num_pion_water = 0;
      int num_true_pion_water = 0;
      for ( auto chain : ev.chains ) {
	if ( chain.particle_flag % 10000 == 13 ) {
	  if ( chain.particle_flag / 10000 == 13 ) muon_correct_id_flag = true;
	}
	else if ( chain.particle_flag % 10000 == 2212 ) {
	  if ( chain.particle_flag / 10000 == 2212 ) num_true_proton_water++;
	  num_proton_water++;
	}
	else if ( chain.particle_flag % 10000 == 211 ) {
	  if ( chain.particle_flag / 10000 == 211 ) num_true_pion_water++;
	  num_pion_water++;
	}
      }

      if ( ev.chains.size() != num_proton_water + 1 ) continue;

      for ( auto chain : ev.chains ) {

	int particle_id = chain.particle_flag % 10000;
	int true_particle_id = chain.particle_flag / 10000;

	// vertex に一番近い basetrack の角度を使う
	double ax = 0.;
	double ay = 0.;
	if ( chain.direction == 1 ) {
	  ax = chain.base.back().ax;
	  ay = chain.base.back().ay;
	}
	else if ( chain.direction == -1 ) {
	  ax = chain.base.front().ax;
	  ay = chain.base.front().ay;
	}
	
	// Neutrino beam の方向の補正を行う
	double thetax = std::atan(ax);
	double thetay = std::atan(ay);
	
	thetax -= neutrino_beam_thetax;
	thetay -= neutrino_beam_thetay;
	
	double tangent = chain.direction * std::hypot(std::tan(thetax), std::tan(thetay));
	double theta = std::atan(tangent);
	double theta_deg = theta * TMath::RadToDeg();
	if ( chain.direction == -1 ) {
	  theta_deg += 180.;
	}
	double cosine = std::cos(theta_deg * TMath::DegToRad());

	// For packing evaluation
	if ( ev.chains.size() == 1 && particle_id == 13 ) {
	  hist_muon_ang_deg_single->Fill(theta_deg, ev.weight);
	  hist_muon_ang_cos_single->Fill(cosine, ev.weight);
	}
	
	if ( muon_correct_id_flag ) {
	  if ( particle_id == 13 ) {
	    hist_muon_ang_deg->Fill(theta_deg, ev.weight);
	    hist_muon_ang_cos->Fill(cosine, ev.weight);
	    hist_mode_muon_ang_deg[mode_id]->Fill(theta_deg, ev.weight);
	    hist_mode_muon_ang_cos[mode_id]->Fill(cosine, ev.weight);
	    hist_flux_mu_deg->Fill(theta_deg, ev.nu_energy / 1000., ev.weight);
	    hist_flux_mu_cos->Fill(cosine, ev.nu_energy / 1000., ev.weight);
	  }
	  else if ( particle_id == 2212 ) {
	    if ( particle_id == true_particle_id ) {
	      hist_proton_ang_deg->Fill(theta_deg, ev.weight);
	      hist_proton_ang_cos->Fill(cosine, ev.weight);
	      hist_mode_proton_ang_deg[mode_id]->Fill(theta_deg, ev.weight);
	      hist_mode_proton_ang_cos[mode_id]->Fill(cosine, ev.weight);
	      hist_flux_p_deg->Fill(theta_deg, ev.nu_energy / 1000., ev.weight);
	      hist_flux_p_cos->Fill(cosine, ev.nu_energy / 1000., ev.weight);
	    }
	    else {
	      hist_proton_misid_ang_deg->Fill(theta_deg, ev.weight);
	      hist_proton_misid_ang_cos->Fill(cosine, ev.weight);
	    }	    
	  }
	  else if ( particle_id == 211 ) {
	    if ( particle_id == true_particle_id ) {
	      hist_pion_ang_deg->Fill(theta_deg, ev.weight);
	      hist_pion_ang_cos->Fill(cosine, ev.weight);
	      hist_mode_pion_ang_deg[mode_id]->Fill(theta_deg, ev.weight);
	      hist_mode_pion_ang_cos[mode_id]->Fill(cosine, ev.weight);
	      hist_flux_pi_deg->Fill(theta_deg, ev.nu_energy / 1000., ev.weight);
	      hist_flux_pi_cos->Fill(cosine, ev.nu_energy / 1000., ev.weight);
	    }
	    else {
	      hist_pion_misid_ang_deg->Fill(theta_deg, ev.weight);
	      hist_pion_misid_ang_cos->Fill(cosine, ev.weight);
	    }
	  }
	  else {
	    hist_other_ang_deg->Fill(theta_deg, ev.weight);
	    hist_other_ang_cos->Fill(cosine, ev.weight);
	  }
	}
	else {
	  if ( particle_id == 13 ) {
	    hist_muon_misid_ang_deg->Fill(theta_deg, ev.weight);
	    hist_muon_misid_ang_cos->Fill(cosine, ev.weight);
	  }
	  else if ( particle_id == 2212 ) {
	    hist_proton_ang_deg_muon_misid->Fill(theta_deg, ev.weight);
	    hist_proton_ang_cos_muon_misid->Fill(cosine, ev.weight);
	  }
	  else if ( particle_id == 211 ) {
	    hist_pion_ang_deg_muon_misid->Fill(theta_deg, ev.weight);
	    hist_pion_ang_cos_muon_misid->Fill(cosine, ev.weight);
	  }
	  else {
	    hist_other_ang_deg->Fill(theta_deg, ev.weight);
	    hist_other_ang_cos->Fill(cosine, ev.weight);
	  }
	}
      }
    }    
  }

  outputfile->cd();
  hist_muon_ang_cos->Write();
  hist_muon_ang_deg->Write();
  hist_pion_ang_cos->Write();
  hist_pion_ang_deg->Write();
  hist_proton_ang_cos->Write();
  hist_proton_ang_deg->Write();
  hist_other_ang_cos->Write();
  hist_other_ang_deg->Write();

  hist_muon_misid_ang_cos->Write();
  hist_muon_misid_ang_deg->Write();
  hist_proton_ang_cos_muon_misid->Write();
  hist_proton_ang_deg_muon_misid->Write();
  hist_pion_ang_cos_muon_misid->Write();
  hist_pion_ang_deg_muon_misid->Write();
  hist_proton_misid_ang_cos->Write();
  hist_proton_misid_ang_deg->Write();
  hist_pion_misid_ang_cos->Write();
  hist_pion_misid_ang_deg->Write();

  hist_muon_ang_cos_single->Write();
  hist_muon_ang_deg_single->Write();
  for ( int i = 0; i < num_ninja_mode; i++ ) {
    hist_mode_muon_ang_cos[i]->Write();
    hist_mode_muon_ang_deg[i]->Write();
    hist_mode_pion_ang_cos[i]->Write();
    hist_mode_pion_ang_deg[i]->Write();
    hist_mode_proton_ang_cos[i]->Write();
    hist_mode_proton_ang_deg[i]->Write();
  }

  hist_flux_mu_deg->Write();
  hist_flux_mu_cos->Write();
  hist_flux_p_deg->Write();
  hist_flux_p_cos->Write();
  hist_flux_pi_deg->Write();
  hist_flux_pi_cos->Write();

  outputfile->Close();

}
