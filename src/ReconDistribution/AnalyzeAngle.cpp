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
				     20, 0., 1.);
  TH1D *hist_muon_ang_deg = new TH1D("hist_muon_ang_deg",
				     "Muon reconstructed angle;#theta_{#mu} [deg];Entries",
				     18, 0., 90.);
  TH1D *hist_pion_ang_cos = new TH1D("hist_pion_ang_cos",
				     "Pion reconstructed angle;cos#theta_{#pi};Entries",
				     40, -1., 1.);
  TH1D *hist_pion_ang_deg = new TH1D("hist_pion_ang_deg",
				     "Pion reconstructed angle;#theta_{#pi} [deg];Entries",
				     36, 0., 180.);
  TH1D *hist_proton_ang_cos = new TH1D("hist_proton_ang_cos",
				       "Proton reconstructed angle;cos#theta_{p};Entries",
				       40, -1., 1.);
  TH1D *hist_proton_ang_deg = new TH1D("hist_proton_ang_deg",
				       "Proton reconstructed angle;#theta_{p} [deg];Entries",
				       36, 0., 180.);

  TH1D *hist_mode_muon_ang_cos[num_ninja_mode];
  TH1D *hist_mode_muon_ang_deg[num_ninja_mode];
  TH1D *hist_mode_pion_ang_cos[num_ninja_mode];
  TH1D *hist_mode_pion_ang_deg[num_ninja_mode];
  TH1D *hist_mode_proton_ang_cos[num_ninja_mode];
  TH1D *hist_mode_proton_ang_deg[num_ninja_mode];
  for ( int i = 0; i < num_ninja_mode; i++ ) {
    hist_mode_muon_ang_cos[i] = new TH1D(Form("hist_muon_ang_cos_%d", i), "", 20, 0., 1.);
    hist_mode_muon_ang_cos[i]->SetFillColor(mode_color[i]);
    hist_mode_muon_ang_cos[i]->SetFillStyle(mode_style[i]);
    hist_mode_muon_ang_deg[i] = new TH1D(Form("hist_muon_ang_deg_%d", i), "", 18, 0., 90.);
    hist_mode_muon_ang_deg[i]->SetFillColor(mode_color[i]);
    hist_mode_muon_ang_deg[i]->SetFillStyle(mode_style[i]);
    hist_mode_pion_ang_cos[i] = new TH1D(Form("hist_pion_ang_cos_%d", i), "", 40, -1., 1.);
    hist_mode_pion_ang_cos[i]->SetFillColor(mode_color[i]);
    hist_mode_pion_ang_cos[i]->SetFillStyle(mode_style[i]);
    hist_mode_pion_ang_deg[i] = new TH1D(Form("hist_pion_ang_deg_%d", i), "", 36, 0., 180.);
    hist_mode_pion_ang_deg[i]->SetFillColor(mode_color[i]);
    hist_mode_pion_ang_deg[i]->SetFillStyle(mode_style[i]);
    hist_mode_proton_ang_cos[i] = new TH1D(Form("hist_proton_ang_cos_%d", i), "", 40, -1., 1.);
    hist_mode_proton_ang_cos[i]->SetFillColor(mode_color[i]);
    hist_mode_proton_ang_cos[i]->SetFillStyle(mode_style[i]);
    hist_mode_proton_ang_deg[i] = new TH1D(Form("hist_proton_ang_deg_%d", i), "", 36, 0., 180.);
    hist_mode_proton_ang_deg[i]->SetFillColor(mode_color[i]);
    hist_mode_proton_ang_deg[i]->SetFillStyle(mode_style[i]);
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
	int particle_id = chain.particle_flag % 10000;
	int true_particle_id = chain.particle_flag / 10000;
	if ( particle_id != true_particle_id ) continue;
	
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
	
	if ( particle_id == 13 ) {
	  hist_muon_ang_cos->Fill(cosine, ev.weight);
	  hist_muon_ang_deg->Fill(theta_deg, ev.weight);
	  hist_mode_muon_ang_cos[mode_id]->Fill(cosine, ev.weight);
	  hist_mode_muon_ang_deg[mode_id]->Fill(theta_deg, ev.weight);	
	}
	else if ( particle_id == 211 ) {
	  hist_pion_ang_cos->Fill(cosine, ev.weight);
	  hist_pion_ang_deg->Fill(theta_deg, ev.weight);
	  hist_mode_pion_ang_cos[mode_id]->Fill(cosine, ev.weight);
	  hist_mode_pion_ang_deg[mode_id]->Fill(theta_deg, ev.weight);
	}
	else if ( particle_id == 2212 ) {
	  hist_proton_ang_cos->Fill(cosine, ev.weight);
	  hist_proton_ang_deg->Fill(theta_deg, ev.weight);
	  hist_mode_proton_ang_cos[mode_id]->Fill(cosine, ev.weight);
	  hist_mode_proton_ang_deg[mode_id]->Fill(theta_deg, ev.weight);
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
  for ( int i = 0; i < num_ninja_mode; i++ ) {
    hist_mode_muon_ang_cos[i]->Write();
    hist_mode_muon_ang_deg[i]->Write();
    hist_mode_pion_ang_cos[i]->Write();
    hist_mode_pion_ang_deg[i]->Write();
    hist_mode_proton_ang_cos[i]->Write();
    hist_mode_proton_ang_deg[i]->Write();
  }
  outputfile->Close();

}
