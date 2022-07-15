#include <cmath>

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/filesystem.hpp>

#include <B2Reader.hh>
#include <B2Enum.hh>
#include <B2EventSummary.hh>
#include <B2VertexSummary.hh>

#include <TFile.h>
#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TVector2.h>
#include <TVector3.h>

#include <McsClass.hpp>

#include "HistogramStyle.hpp"
#include "DrawConst.hpp"

#include "Analyze0pi1p.hpp"

namespace logging = boost::log;
namespace fs = boost::filesystem;

void Analyze0pi1p(std::string b2filename,
		  std::string momchfilename,
		  std::string outputfilename) {

  BOOST_LOG_TRIVIAL(info) << "==========0pi1p mode==========";

  // input B2 file
  B2Reader reader(b2filename);

  if ( !fs::exists(momchfilename) ) {
    throw std::runtime_error("File not found : " + momchfilename);
  }
  auto ev_vec = Momentum_recon::ReadEventInformationBin(momchfilename);

  // output file
  TFile *outputfile = new TFile((TString)outputfilename, "recreate");
  BOOST_LOG_TRIVIAL(info) << "Output filename : " << outputfilename;

  // total
  TH1D *hist_muon_mom = new TH1D("hist_muon_mom",
				 "Muon reconstructed momentum;p_{#mu} [MeV/c];Entries",
				 50, 0., 2000.);
  TH1D *hist_proton_mom = new TH1D("hist_proton_mom",
				   "Proton reconstructed momentum;p_{p} [MeV/c];Entries",
				   50, 0., 1500.);
  TH1D *hist_muon_ang = new TH1D("hist_muon_ang",
				 "Muon reconstructed angle;#theta_{#mu} [deg];Entries",
				 18, 0., 90.);
  TH1D *hist_muon_cos = new TH1D("hist_muon_cos",
				 "Muon reconstructed angle;cos#theta_{#mu};Entries",
				 20, 0., 1.);
  TH1D *hist_proton_ang = new TH1D("hist_proton_ang",
				   "Proton reconstructed angle;#theta_{#mu} [deg];Entries",
				   36, 0., 180.);
  TH1D *hist_proton_cos = new TH1D("hist_proton_cos",
				   "Proton reconstructed angle;cos#theta_{#mu};Entries",
				   40, -1., 1.);
  TH2D *hist_muon_mom_ang = new TH2D("hist_muon_mom_ang",
				     "Muon momentum vs angle;p_{#mu} [MeV/c];#theta_{#mu} [deg]",
				     50, 0., 2000., 18, 0., 90.);
  TH2D *hist_muon_mom_cos = new TH2D("hist_muon_mom_cos",
				     "Muon momentum vs angle;p_{#mu} [MeV/c];cos#theta_{#mu}",
				     50, 0., 2000., 20, 0., 1.);
  TH2D *hist_proton_mom_ang = new TH2D("hist_proton_mom_ang",
				       "Proton momentum vs angle;p_{p} [MeV/c];#theta_{p} [deg]",
				       50., 0., 1500., 36, 0., 180.);
  TH2D *hist_proton_mom_cos = new TH2D("hist_proton_mom_cos",
				       "Proton momentum vs angle;p_{p} [MeV/c];cos#theta_{p}",
				       50, 0., 1500., 40, -1., 1.);
  TH1D *hist_q2 = new TH1D("hist_q2",
			   "Q^{2};Q^{2} [MeV^{2}/c^{2}];Entries",
			   100, 0, 1e6);
  TH1D *hist_nu_ene_bias = new TH1D("hist_nu_ene_bias",
				    "Energy bias assume CCQE;(E_{rec} - E_{true})/E_{true};Entries",
				    100, -1., 1.);
  TH1D *hist_nu_ene_recon = new TH1D("hist_nu_ene_recon",
				     "Reconstructed neutrino energy;E_{rec} [MeV];Entries",
				     100, 0, 2000);
  TH2D *hist_nu_ene_recon_true = new TH2D("hist_nu_ene_recon_true",
					  "Neutrino energy;E_{rec} [MeV];E_{true} [MeV]",
					  100, 0, 2000, 100, 0, 2000);

  TH1D *hist_muon_mom_mcs = new TH1D("hist_muon_mom_mcs",
				     "Muon MCS momentum;p_{#mu, MCS} [MeV/c];Entries",
				     50, 0., 2000.);
  TH1D *hist_muon_mom_range = new TH1D("hist_muon_mom_range",
				       "Muon rangemomentum;p_{#mu, range} [MeV/c];Entries",
				       50, 0., 2000.);
  TH1D *hist_proton_mom_mcs = new TH1D("hist_proton_mom_mcs",
				       "Proton MCS momentum;p_{p, MCS} [MeV/c];Entries",
				       50, 0., 1500.);
  TH1D *hist_proton_mom_range = new TH1D("hist_proton_mom_range",
					 "Proton range momentum;p_{p, range} [MeV/c];Entries",
					 50, 0., 1500.);
  TH2D *hist_muon_mom_ang_mcs = new TH2D("hist_muon_mom_ang_mcs",
					 "Muon momentum vs angle;p_{#mu, MCS} [MeV/c];#theta_{#mu} [deg]",
					 50, 0., 2000., 18, 0., 90.);
  TH2D *hist_muon_mom_cos_mcs = new TH2D("hist_muon_mom_cos_mcs",
					 "Muon momentum vs angle;p_{#mu, MCS} [MeV/c];cos#theta_{#mu}",
					 50, 0., 2000., 18, 0., 90.);
  TH2D *hist_muon_mom_ang_range = new TH2D("hist_muon_mom_ang_range",
					 "Muon momentum vs angle;p_{#mu, range} [MeV/c];#theta_{#mu} [deg]",
					 50, 0., 2000., 18, 0., 90.);
  TH2D *hist_muon_mom_cos_range = new TH2D("hist_muon_mom_cos_range",
					 "Muon momentum vs angle;p_{#mu, range} [MeV/c];cos#theta_{#mu}",
					 50, 0., 2000., 18, 0., 90.);
  TH2D *hist_proton_mom_ang_mcs = new TH2D("hist_proton_mom_ang_mcs",
					   "Proton momentum vs angle;p_{p, MCS} [MeV/c];#theta_{p} [deg]",
					   50, 0., 1500., 36, 0., 180.);
  TH2D *hist_proton_mom_cos_mcs = new TH2D("hist_proton_mom_cos_mcs",
					   "Proton momentum vs angle;p_{p, MCS} [MeV/c];cos#theta_{p}",
					   50, 0., 1500., 36, 0., 180.);
  TH2D *hist_proton_mom_ang_range = new TH2D("hist_proton_mom_ang_range",
					     "Proton momentum vs angle;p_{p, range} [MeV/c];#theta_{p} [deg]",
					     50, 0., 1500., 36, 0., 180.);
  TH2D *hist_proton_mom_cos_range = new TH2D("hist_proton_mom_cos_range",
					     "Proton momentum vs angle;p_{p, range} [MeV/c];cos#theta_{p}",
					     50, 0., 1500., 36, 0., 180.);
  
  TH2D *hist_muon_mom_recon_true = new TH2D("hist_muon_mom_recon_true",
					    "Muon momentum;p_{#mu, true} [MeV/c];p_{#mu, recon} [MeV/c]",
					    50, 0., 2000., 50, 0., 2000.);
  TH2D *hist_muon_mom_recon_true_mcs = new TH2D("hist_muon_mom_recon_true_mcs",
						"Muon MCS momentum;p_{#mu, true} [MeV/c];p_{#mu, recon} [MeV/c]",
						50, 0., 2000., 50, 0., 2000.); 
  TH2D *hist_muon_mom_recon_true_range = new TH2D("hist_muon_mom_recon_true_range",
						  "Muon rangemomentum;p_{#mu, true} [MeV/c];p_{#mu, recon} [MeV/c]",
						  50, 0., 2000., 50, 0., 2000.);
  TH2D *hist_proton_mom_recon_true = new TH2D("hist_proton_mom_recon_true",
					      "Proton momentum;p_{p, true} [MeV/c];p_{p, recon} [MeV/c]",
					      50, 0., 1500., 50, 0., 1500.);
  TH2D *hist_proton_mom_recon_true_mcs = new TH2D("hist_proton_mom_recon_true_mcs",
						  "Proton MCS momentum;p_{p, true} [MeV/c];p_{p, recon} [MeV/c]",
						  50, 0., 1500., 50, 0., 1500.); 
  TH2D *hist_proton_mom_recon_true_range = new TH2D("hist_proton_mom_recon_true_range",
						    "Proton rangemomentum;p_{p, true} [MeV/c];p_{p, recon} [MeV/c]",
						    50, 0., 1500., 50, 0., 1500.);

  // TKI
  TH1D *hist_dpt = new TH1D("hist_dpt",
			    "#deltap_{T};#deltap_{T} [MeV/c];Entries",
			    200, 0., 2000.);
  TH1D *hist_dalphat = new TH1D("hist_dalphat",
				"#delta#alpha_{T};#delta#alpha_{T};Entries",
				90, 0., 180.);
  TH1D *hist_cosdat = new TH1D("hist_cosdat",
			       "cos#delta#alpha_{T};cos#delta#alpha_{T};Entries",
			       100, -1., 1.);
  TH1D *hist_dphit = new TH1D("hist_dphit",
			      "#delta#phi_{T};#delta#phi_{T};Entries",
			      90, 0., 180.);
  TH1D *hist_cosdphit = new TH1D("hist_cosdphit",
				 "cos#delta#phi_{T};cos#delta#phi_{T};Entries",
				 100, -1., 1.);
  TH1D *hist_dptx = new TH1D("hist_dptx",
			     "#deltap_{Tx};#deltap_{Tx} [MeV/c];Entries",
			     200, -2., 2.);
  TH1D *hist_dpty = new TH1D("hist_dpty",
			     "#deltap_{Ty};#deltap_{Ty}[MeV/c];Entries",
			     200, -2., 2.);

  // mode
  TH1D *hist_mode_muon_mom[num_ninja_mode];
  TH1D *hist_mode_proton_mom[num_ninja_mode];
  TH1D *hist_mode_muon_ang[num_ninja_mode];
  TH1D *hist_mode_muon_cos[num_ninja_mode];
  TH1D *hist_mode_proton_ang[num_ninja_mode];
  TH1D *hist_mode_proton_cos[num_ninja_mode];
  TH1D *hist_mode_q2[num_ninja_mode];
  TH1D *hist_mode_nu_ene_bias[num_ninja_mode];
  TH1D *hist_mode_nu_ene_recon[num_ninja_mode];
  TH1D *hist_mode_dpt[num_ninja_mode];
  TH1D *hist_mode_dalphat[num_ninja_mode];
  TH1D *hist_mode_cosdat[num_ninja_mode];
  TH1D *hist_mode_dphit[num_ninja_mode];
  TH1D *hist_mode_cosdphit[num_ninja_mode];
  TH1D *hist_mode_dptx[num_ninja_mode];
  TH1D *hist_mode_dpty[num_ninja_mode];
  for ( int i = 0; i < num_ninja_mode; i++ ) {
    hist_mode_muon_mom[i] = new TH1D(Form("hist_muon_mom_%d", i), "", 15, 0., 1500.);
    hist_mode_muon_mom[i]->SetFillColor(mode_color[i]);
    hist_mode_muon_mom[i]->SetFillStyle(mode_style[i]);
    hist_mode_proton_mom[i] = new TH1D(Form("hist_proton_mom_%d", i), "", 15, 0., 1500.);
    hist_mode_proton_mom[i]->SetFillColor(mode_color[i]);
    hist_mode_proton_mom[i]->SetFillStyle(mode_style[i]);
    hist_mode_muon_ang[i] = new TH1D(Form("hist_muon_ang_%d", i), "", 18, 0., 90.);
    hist_mode_muon_ang[i]->SetFillColor(mode_color[i]);
    hist_mode_muon_ang[i]->SetFillStyle(mode_style[i]);
    hist_mode_muon_cos[i] = new TH1D(Form("hist_muon_cos_%d", i), "", 20, 0., 1.);
    hist_mode_muon_cos[i]->SetFillColor(mode_color[i]);
    hist_mode_muon_cos[i]->SetFillStyle(mode_style[i]);
    hist_mode_proton_ang[i] = new TH1D(Form("hist_proton_ang_%d", i), "", 36, 0., 180.);
    hist_mode_proton_ang[i]->SetFillColor(mode_color[i]);
    hist_mode_proton_ang[i]->SetFillStyle(mode_style[i]);
    hist_mode_proton_cos[i] = new TH1D(Form("hist_proton_cos_%d", i), "", 40, -1., 1.);
    hist_mode_proton_cos[i]->SetFillColor(mode_color[i]);
    hist_mode_proton_cos[i]->SetFillStyle(mode_style[i]);
    hist_mode_q2[i] = new TH1D(Form("hist_q2_%d", i), "", 100, 0, 1e6);
    hist_mode_q2[i]->SetFillColor(mode_color[i]);
    hist_mode_q2[i]->SetFillStyle(mode_style[i]);
    hist_mode_nu_ene_bias[i] = new TH1D(Form("hist_nu_ene_bias_%d", i), "", 100, -1., 1.);
    hist_mode_nu_ene_bias[i]->SetFillColor(mode_color[i]);
    hist_mode_nu_ene_bias[i]->SetFillStyle(mode_style[i]);
    hist_mode_nu_ene_recon[i] = new TH1D(Form("hist_nu_ene_recon_%d", i), "", 100, 0, 2000);
    hist_mode_nu_ene_recon[i]->SetFillColor(mode_color[i]);
    hist_mode_nu_ene_recon[i]->SetFillStyle(mode_style[i]);
    hist_mode_dpt[i] = new TH1D(Form("hist_dpt_%d", i), "", 200, 0., 2000.);
    hist_mode_dpt[i]->SetFillColor(mode_color[i]);
    hist_mode_dpt[i]->SetFillStyle(mode_style[i]);
    hist_mode_dalphat[i] = new TH1D(Form("hist_dalphat_%d", i), "", 90, 0., 180.);
    hist_mode_dalphat[i]->SetFillColor(mode_color[i]);
    hist_mode_dalphat[i]->SetFillStyle(mode_style[i]);
    hist_mode_cosdat[i] = new TH1D(Form("hist_cosdat_%d", i), "", 100, -1., 1.);
    hist_mode_cosdat[i]->SetFillColor(mode_color[i]);
    hist_mode_cosdat[i]->SetFillStyle(mode_style[i]);
    hist_mode_dphit[i] = new TH1D(Form("hist_dphit_%d", i), "", 90, 0., 180.);
    hist_mode_dphit[i]->SetFillColor(mode_color[i]);
    hist_mode_dphit[i]->SetFillStyle(mode_style[i]);
    hist_mode_cosdphit[i] = new TH1D(Form("hist_cosdphit_%d", i), "", 100, -1., 1.);
    hist_mode_cosdphit[i]->SetFillColor(mode_color[i]);
    hist_mode_cosdphit[i]->SetFillStyle(mode_style[i]);
    hist_mode_dptx[i] = new TH1D(Form("hist_dptx_%d", i), "", 200, -2000., 2000.);
    hist_mode_dptx[i]->SetFillColor(mode_color[i]);
    hist_mode_dptx[i]->SetFillStyle(mode_style[i]);
    hist_mode_dpty[i] = new TH1D(Form("hist_dpty_%d", i), "", 200, -2000., 2000.);
    hist_mode_dpty[i]->SetFillColor(mode_color[i]);
    hist_mode_dpty[i]->SetFillStyle(mode_style[i]);
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
      int num_muon = 0;
      int num_proton = 0;
      int num_pion = 0;
      
      bool muon_stop_flag = false;
      bool proton_stop_flag = false;
      TVector3 muon_tangent;
      TVector3 proton_tangent;
      int proton_direction;
      double muon_momentum;
      double proton_momentum;
      TVector3 muon_momentum_vec;
      TVector3 proton_momentum_vec;
      double true_muon_momentum;
      double true_proton_momentum;

      for ( auto chain : ev.chains ) {

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

	double thetax = std::atan(ax);
	double thetay = std::atan(ay);

	thetax -= neutrino_beam_thetax;
	thetay -= neutrino_beam_thetay;
	
	ax = std::tan(thetax);
	ay = std::tan(thetay);

	int particle_id = chain.particle_flag % 10000;
	int true_particle_id = chain.particle_flag / 10000;
	if ( particle_id != true_particle_id ) continue;

	for ( auto true_chain : ev.true_chains ) {
	  if ( chain.chainid == true_chain.chainid ) {
	    if ( particle_id == 13 )
	      true_muon_momentum = true_chain.bm_range_mom;
	    else if ( particle_id == 2212 )
	      true_proton_momentum = true_chain.bm_range_mom;
	    break;
	  }
	}

	if ( particle_id == 13 ) {
	  num_muon++;
	  muon_tangent.SetXYZ(ax, ay, 1.);
	  if ( chain.stop_flag == 1 ) {
	    muon_momentum = chain.bm_range_mom;
	    muon_stop_flag = true;
	  }
	  else if ( chain.stop_flag == 0 ) {
	    muon_momentum = chain.ecc_mcs_mom[0];
	  }
	  muon_momentum_vec = (muon_momentum / muon_tangent.Mag()) * muon_tangent;
	}
	else if ( particle_id == 2212 ) {
	  num_proton++;
	  if ( chain.direction == 1 ) {
	    proton_tangent.SetXYZ(ax, ay, 1.);
	    proton_direction = 1;
	  }
	  else if ( chain.direction == -1 ) {
	    proton_tangent.SetXYZ(ax, ay, 1.);
	    proton_direction = -1;
	  }
	  if ( chain.stop_flag == 2 ) {
	    proton_momentum = chain.ecc_range_mom[1];
	    proton_stop_flag = true;
	  }
	  else {
	    proton_momentum = chain.ecc_mcs_mom[1];
	  }
	  proton_momentum_vec = (proton_momentum * chain.direction / proton_tangent.Mag()) * proton_tangent;
	}
	else if ( particle_id == 211 ) {
	  num_pion++;
	}
      }

      if ( ev.chains.size() != 2 ) continue;
      if ( num_muon != 1 || num_proton != 1 || num_pion != 0 ) continue;

      hist_muon_mom->Fill(muon_momentum, ev.weight);
      hist_mode_muon_mom[mode_id]->Fill(muon_momentum, ev.weight);
      hist_proton_mom->Fill(proton_momentum, ev.weight);
      hist_mode_proton_mom[mode_id]->Fill(proton_momentum, ev.weight);

      double muon_rad = std::atan(std::hypot(muon_tangent.X(), muon_tangent.Y()));
      double muon_ang = muon_rad * TMath::RadToDeg();
      hist_muon_ang->Fill(muon_ang, ev.weight);
      hist_mode_muon_ang[mode_id]->Fill(muon_ang, ev.weight);
      hist_muon_cos->Fill(std::cos(muon_rad), ev.weight);
      hist_mode_muon_cos[mode_id]->Fill(std::cos(muon_rad), ev.weight);
      
      double proton_rad = std::atan(proton_direction * std::hypot(proton_tangent.X(), proton_tangent.Y()));
      double proton_ang = proton_rad * TMath::RadToDeg();
      if ( proton_direction == -1 ) {
	proton_ang += 180.;
      }
      hist_proton_ang->Fill(proton_ang, ev.weight);
      hist_mode_proton_ang[mode_id]->Fill(proton_ang, ev.weight);
      hist_proton_cos->Fill(std::cos(proton_ang * TMath::DegToRad()));
      hist_mode_proton_cos[mode_id]->Fill(std::cos(proton_ang * TMath::DegToRad()));

      hist_muon_mom_ang->Fill(muon_momentum, muon_ang, ev.weight);
      hist_muon_mom_cos->Fill(muon_momentum, std::cos(muon_rad), ev.weight);
      hist_proton_mom_ang->Fill(proton_momentum, proton_ang, ev.weight);
      hist_proton_mom_cos->Fill(proton_momentum, std::cos(proton_ang * TMath::DegToRad()), ev.weight);

      hist_muon_mom_recon_true->Fill(true_muon_momentum, muon_momentum, ev.weight);
      hist_proton_mom_recon_true->Fill(true_proton_momentum, proton_momentum, ev.weight);

      if ( muon_stop_flag ) {
	hist_muon_mom_range->Fill(muon_momentum, ev.weight);
	hist_muon_mom_ang_range->Fill(muon_momentum, muon_ang, ev.weight);
	hist_muon_mom_cos_range->Fill(muon_momentum, std::cos(muon_rad), ev.weight);
	hist_muon_mom_recon_true_range->Fill(true_muon_momentum, muon_momentum, ev.weight);
      }
      else {
	hist_muon_mom_mcs->Fill(muon_momentum, ev.weight);
	hist_muon_mom_ang_mcs->Fill(muon_momentum, muon_ang, ev.weight);
	hist_muon_mom_cos_mcs->Fill(muon_momentum, std::cos(muon_rad), ev.weight);
	hist_muon_mom_recon_true_mcs->Fill(true_muon_momentum, muon_momentum, ev.weight);
      }

      if ( proton_stop_flag ) {
	hist_proton_mom_range->Fill(proton_momentum, ev.weight);
	hist_proton_mom_ang_range->Fill(proton_momentum, proton_ang, ev.weight);
	hist_proton_mom_cos_range->Fill(proton_momentum, std::cos(proton_ang * TMath::RadToDeg()), ev.weight);
	hist_proton_mom_recon_true_range->Fill(true_proton_momentum, proton_momentum, ev.weight);
      }
      else {
	hist_proton_mom_mcs->Fill(proton_momentum, ev.weight);
	hist_proton_mom_ang_mcs->Fill(proton_momentum, proton_ang, ev.weight);
	hist_proton_mom_cos_mcs->Fill(proton_momentum, std::cos(proton_ang * TMath::RadToDeg()), ev.weight);
	hist_proton_mom_recon_true_mcs->Fill(true_proton_momentum, proton_momentum, ev.weight);
      }             
      
      double muon_energy = std::sqrt(muon_momentum * muon_momentum + muon_mass * muon_mass);      
      double recon_nu_ene = (neutron_mass * muon_energy - muon_mass * muon_mass / 2.
			     + (proton_mass * proton_mass - neutron_mass * neutron_mass) / 2.)
	/ (neutron_mass - muon_energy + muon_momentum * std::cos(muon_rad));
      double recon_q2 = - muon_mass * muon_mass + 2 * recon_nu_ene * (muon_energy - muon_momentum * std::cos(muon_rad));
      hist_q2->Fill(recon_q2, ev.weight);
      hist_mode_q2[mode_id]->Fill(recon_q2, ev.weight);
      hist_nu_ene_bias->Fill((recon_nu_ene - ev.nu_energy) / ev.nu_energy, ev.weight);
      hist_mode_nu_ene_bias[mode_id]->Fill((recon_nu_ene - ev.nu_energy) / ev.nu_energy, ev.weight);
      hist_nu_ene_recon->Fill(recon_nu_ene, ev.weight);
      hist_mode_nu_ene_recon[mode_id]->Fill(recon_nu_ene, ev.weight);
      hist_nu_ene_recon_true->Fill(recon_nu_ene, ev.nu_energy, ev.weight);

      TVector2 mumom_vec_2d(muon_momentum_vec.X(), muon_momentum_vec.Y());
      TVector2 promom_vec_2d(proton_momentum_vec.X(), proton_momentum_vec.Y());
      TVector2 dpt_vec = mumom_vec_2d + promom_vec_2d;
      Double_t dalphat = std::acos((-1. * mumom_vec_2d * dpt_vec) / mumom_vec_2d.Mod() / dpt_vec.Mod());
      hist_dpt->Fill(dpt_vec.Mod(), ev.weight);
      hist_mode_dpt[mode_id]->Fill(dpt_vec.Mod(), ev.weight);
      hist_dalphat->Fill(dalphat * TMath::RadToDeg(),
			 ev.weight);
      hist_mode_dalphat[mode_id]->Fill(dalphat * TMath::RadToDeg(),
				  ev.weight);
      hist_cosdat->Fill((-1. * mumom_vec_2d * dpt_vec) / mumom_vec_2d.Mod() / dpt_vec.Mod(), ev.weight);
      hist_mode_cosdat[mode_id]->Fill((-1. * mumom_vec_2d * dpt_vec) / mumom_vec_2d.Mod() / dpt_vec.Mod(), ev.weight);
      hist_dphit->Fill(std::acos((-1. * mumom_vec_2d * promom_vec_2d) / mumom_vec_2d.Mod() / promom_vec_2d.Mod()) * TMath::RadToDeg(),
		       ev.weight);
      hist_mode_dphit[mode_id]->Fill(std::acos((-1. * mumom_vec_2d * promom_vec_2d) / mumom_vec_2d.Mod() / promom_vec_2d.Mod()) * TMath::RadToDeg(),
				     ev.weight);
      hist_cosdphit->Fill((-1. * mumom_vec_2d * promom_vec_2d) / mumom_vec_2d.Mod() / promom_vec_2d.Mod(),
			  ev.weight);
      hist_mode_cosdphit[mode_id]->Fill((-1. * mumom_vec_2d * promom_vec_2d) / mumom_vec_2d.Mod() / promom_vec_2d.Mod(),
					ev.weight);
      hist_dptx->Fill(TMath::Sign(1., mumom_vec_2d.X()) * dpt_vec.Mod() * std::sin(dalphat),
		      ev.weight);
      hist_mode_dptx[mode_id]->Fill(TMath::Sign(1., mumom_vec_2d.X()) * dpt_vec.Mod() * std::sin(dalphat),
				    ev.weight);
      hist_dpty->Fill(dpt_vec.Mod() * std::cos(dalphat), ev.weight);
      hist_mode_dpty[mode_id]->Fill(dpt_vec.Mod() * std::cos(dalphat), ev.weight);
    }
 
  }

  outputfile->cd();
  hist_muon_mom->Write();
  hist_proton_mom->Write();
  hist_muon_ang->Write();
  hist_muon_cos->Write();
  hist_proton_ang->Write();
  hist_proton_cos->Write();
  hist_muon_mom_ang->Write();
  hist_muon_mom_cos->Write();
  hist_proton_mom_ang->Write();
  hist_proton_mom_cos->Write();
  hist_q2->Write();
  hist_nu_ene_bias->Write();
  hist_nu_ene_recon->Write();
  hist_nu_ene_recon_true->Write();
  hist_muon_mom_mcs->Write();
  hist_muon_mom_range->Write();
  hist_proton_mom_range->Write();
  hist_proton_mom_mcs->Write();
  hist_muon_mom_ang_mcs->Write();
  hist_muon_mom_cos_mcs->Write();
  hist_muon_mom_ang_range->Write();
  hist_muon_mom_cos_range->Write();
  hist_proton_mom_ang_mcs->Write();
  hist_proton_mom_cos_mcs->Write();
  hist_proton_mom_ang_range->Write();
  hist_proton_mom_cos_range->Write();
  hist_muon_mom_recon_true->Write();
  hist_muon_mom_recon_true_mcs->Write();
  hist_muon_mom_recon_true_range->Write();
  hist_proton_mom_recon_true->Write();
  hist_proton_mom_recon_true_mcs->Write();
  hist_proton_mom_recon_true_range->Write();
  hist_dpt->Write();
  hist_dalphat->Write();
  hist_cosdat->Write();
  hist_dphit->Write();
  hist_cosdphit->Write();
  hist_dptx->Write();
  hist_dpty->Write();
  for ( int i = 0; i < num_ninja_mode; i++ ) {
    hist_mode_muon_mom[i]->Write();
    hist_mode_proton_mom[i]->Write();
    hist_mode_muon_ang[i]->Write();
    hist_mode_muon_cos[i]->Write();
    hist_mode_proton_ang[i]->Write();
    hist_mode_proton_cos[i]->Write();
    hist_mode_q2[i]->Write();
    hist_mode_nu_ene_bias[i]->Write();
    hist_mode_nu_ene_recon[i]->Write();
    hist_mode_dpt[i]->Write();
    hist_mode_dalphat[i]->Write();
    hist_mode_cosdat[i]->Write();
    hist_mode_dphit[i]->Write();
    hist_mode_cosdphit[i]->Write();
    hist_mode_dptx[i]->Write();
    hist_mode_dpty[i]->Write();
  }
  outputfile->Close();

}
