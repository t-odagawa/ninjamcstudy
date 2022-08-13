#include <iostream>
#include <vector>
#include <string>

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/filesystem.hpp>

#include <TRandom.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>

#include <B2Reader.hh>
#include <B2SpillSummary.hh>
#include <B2EventSummary.hh>
#include <B2VertexSummary.hh>

#include <McsClass.hpp>
#include <McsConst.hpp>
#include <McsFunction.hpp>
#include <PidData.hpp>
#include <PidFunction.hpp>

#include "AnalyzePid.hpp"

namespace logging = boost::log;
namespace fs = boost::filesystem;

void AnalyzePid(std::string b2filename,
		std::string momchfilename,
		std::string outputfilename) {

  BOOST_LOG_TRIVIAL(info) << "==========Particle ID mode==========";

  // input B2 file
  B2Reader reader(b2filename);
   
  // input momch file
  if ( !fs::exists(momchfilename) ) {
    throw std::runtime_error("File not found : " + momchfilename);
  }
  auto ev_vec = Momentum_recon::ReadEventInformationBin(momchfilename);

  // outputfilename
  TFile *outputfile = new TFile((TString)outputfilename, "recreate");
  BOOST_LOG_TRIVIAL(info) << "Output filename : " << outputfilename;

  TH1D *hist_likelihood_ratio = new TH1D("hist_likelihood_ratio", "", 20, 0, 1);
  TH1D *hist_likelihood_ratio_pi = new TH1D("hist_likelihood_ratio_pi", "", 20, 0, 1);
  TH1D *hist_likelihood_ratio_p = new TH1D("hist_likelihood_ratio_p", "", 20, 0, 1);

  TH2D *hist_vph_mom_ang[18];
  TH1D *hist_vph_ang_lowmom[18];
  TH1D *hist_likelihood_ratio_ang[18];
  TH1D *hist_likelihood_ratio_ang_pi[18];
  TH1D *hist_likelihood_ratio_ang_p[18];

  for ( int i = 0; i < 18; i++ ) {
    double ang_min, ang_max;
    if ( i < 7 ) {
      ang_min = 0.1 * i;
      ang_max = 0.1 * (i + 1);
    }
    else if ( i < 11 ) {
      ang_min = 0.7 + 0.2 * (i - 7);
      ang_max = 0.7 + 0.2 * (i - 7 + 1);
    }
    else if ( i < 14 ) {
      ang_min = 1.5 + 0.4 * (i - 11);
      ang_max = 1.5 + 0.4 * (i - 11 + 1);
    }
    else {
      ang_min = 2.7 + 0.6 * (i - 14);
      ang_max = 2.7 + 0.6 * (i - 14 + 1);
    }
    
    BOOST_LOG_TRIVIAL(debug) << ang_min << ", " << ang_max;
  }

  for ( int i = 0; i < 18; i++ ) {
    hist_vph_mom_ang[i] = new TH2D(Form("hist_vph_mom_ang_%d", i), "",
				   100, 0, 1500, 350, 0, 350);
    hist_vph_ang_lowmom[i] = new TH1D(Form("hist_vph_ang_lowmom_%d", i), "",
				      50, 0, 350);
    hist_likelihood_ratio_ang[i] = new TH1D(Form("hist_likelihood_ratio_ang_%d", i), "",
					 20, 0, 1);
    hist_likelihood_ratio_ang_pi[i] = new TH1D(Form("hist_likelihood_ratio_ang_pi_%d", i), "",
					 20, 0, 1);
    hist_likelihood_ratio_ang_p[i] = new TH1D(Form("hist_likelihood_ratio_ang_p_%d", i), "",
					 20, 0, 1);
  }

  const std::string data_path = "/home/t2k/odagawa/NinjaMomentumRecon/data";
  PidData pid_data_(data_path);
  PidFunction pid_function_(pid_data_);

  for ( auto ev : ev_vec ) {
    
    reader.ReadSpill(ev.groupid);
    auto &spill_summary = reader.GetSpillSummary();
    auto it_event = spill_summary.BeginTrueEvent();
    const auto *event = it_event.Next();

    auto &vertex = event->GetPrimaryVertex();
    // int mode_id = GetNinjaModeId(vertex.GetInteractionType());

    if ( ev.vertex_material != B2Material::kWater ) continue;

    if ( !ev.chains.empty() ) {
      for ( auto chain : ev.chains ) {

	// particle id (only reconstructed pion or proton)
	int true_particle_id = chain.particle_flag / 10000;
	int recon_particle_id = chain.particle_flag % 10000;
	if ( recon_particle_id != 2212 &&
	     TMath::Abs(recon_particle_id) != 211 ) continue;

	double vph = chain.base.front().m[0].ph;
	double ax, ay;
	if ( chain.direction == 1 ) {
	  ax = chain.base.back().ax;
	  ay = chain.base.back().ay;
	}
	else if ( chain.direction == -1 ) {
	  ax = chain.base.front().ax;
	  ay = chain.base.front().ay;
	}
	else continue; // direction should be +/-1
	double pbeta = chain.ecc_mcs_mom[0];
	pbeta = CalculatePBetaFromMomentum(pbeta, MCS_MUON_MASS);
	double tangent = TMath::Hypot(ax, ay);

	// likelihood ratio
	double vph_mean_pi = pid_function_.CalcVphMuon(pbeta, tangent);
	double vph_mean_p = pid_function_.CalcVphProton(pbeta, tangent);

	double l_ratio;
	if ( vph < vph_mean_pi ) l_ratio = 0.999;
	else if ( vph > vph_mean_p ) l_ratio = 0.;
	else {
	  l_ratio = chain.muon_likelihood / (chain.muon_likelihood + chain.proton_likelihood);
	}

	hist_likelihood_ratio->Fill(l_ratio, ev.weight);
	if ( true_particle_id == 2212 ) 
	  hist_likelihood_ratio_p->Fill(l_ratio, ev.weight);
	else if ( TMath::Abs(true_particle_id) == 211 )
	  hist_likelihood_ratio_pi->Fill(l_ratio, ev.weight);

	for ( int i = 0; i < 18; i++ ) {
	  double ang_min, ang_max;
	  if ( i < 7 ) {
	    ang_min = 0.1 * i;
	    ang_max = 0.1 * (i + 1);
	  }
	  else if ( i < 11 ) {
	    ang_min = 0.7 + 0.2 * (i - 7);
	    ang_max = 0.7 + 0.2 * (i - 7 + 1);
	  }
	  else if ( i < 14 ) {
	    ang_min = 1.5 + 0.4 * (i - 11);
	    ang_max = 1.5 + 0.4 * (i - 11 + 1);
	  }
	  else {
	    ang_min = 2.7 + 0.6 * (i - 14);
	    ang_max = 2.7 + 0.6 * (i - 14 + 1);
	  }

	  if ( ang_min <= tangent &&
	       tangent < ang_max ) {
	    hist_vph_mom_ang[i]->Fill(chain.ecc_mcs_mom[0], vph, ev.weight);
	    if ( chain.ecc_mcs_mom[0] >= 200. &&
		 chain.ecc_mcs_mom[0] < 300. ) {
	      hist_vph_ang_lowmom[i]->Fill(vph, ev.weight);
	    }
	    hist_likelihood_ratio_ang[i]->Fill(l_ratio, ev.weight);
	    if ( true_particle_id == 2212 ) 
	      hist_likelihood_ratio_ang_p[i]->Fill(l_ratio, ev.weight);
	    else if ( TMath::Abs(true_particle_id) == 211 )
	      hist_likelihood_ratio_ang_pi[i]->Fill(l_ratio, ev.weight);
	    break;
	  }

	}

	
      }
    }

  }

  outputfile->cd();
  hist_likelihood_ratio->Write();
  hist_likelihood_ratio_pi->Write();
  hist_likelihood_ratio_p->Write();
  for ( int i = 0; i < 18; i++ ) {
    hist_vph_mom_ang[i]->Write();
    hist_vph_ang_lowmom[i]->Write();
    hist_likelihood_ratio_ang[i]->Write();
    hist_likelihood_ratio_ang_pi[i]->Write();
    hist_likelihood_ratio_ang_p[i]->Write();
  }
  outputfile->Close();

}
