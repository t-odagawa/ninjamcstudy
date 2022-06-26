#include <iostream>
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

#include "Analyze0pi2p.hpp"

namespace logging = boost::log;
namespace fs = boost::filesystem;

void Analyze0pi2p(std::string b2filename,
		  std::string momchfilename,
		  std::string outputfilename) {

  BOOST_LOG_TRIVIAL(info) << "==========0pi2p mode==========";

  // input B2 file
  B2Reader reader(b2filename);

  if ( !fs::exists(momchfilename) ) {
    throw std::runtime_error("File not found : " + momchfilename);
  }
  auto ev_vec = Momentum_recon::ReadEventInformationBin(momchfilename);

  // outpu file
  TFile *outputfile = new TFile((TString)outputfilename, "recreate");
  BOOST_LOG_TRIVIAL(info) << "Output filename : " << outputfilename;

  // total
  TH1D *hist_muon_mom = new TH1D("hist_muon_mom",
				 "Muon reconstructed momentum;p_{#mu} [MeV/c];Entries",
				 50, 0., 1500.);
  TH1D *hist_proton_mom = new TH1D("hist_proton_mom",
				   "Proton reconstructed momentum;p_{p} [MeV/c];Entries",
				   50, 0., 1500.);
  TH1D *hist_proton_mom_high = new TH1D("hist_proton_mom_high",
					"Higher proton momentum;p_{p1} [MeV/c];Entries",
					50, 0., 1500.);
  TH1D *hist_proton_mom_low = new TH1D("hist_proton_mom_low",
				       "Lower proton momentum;p_{p2} [MeV/c];Entries",
				       50, 0., 1500.);
  TH1D *hist_muon_ang = new TH1D("hist_muon_ang",
				 "Muon reconstructed angle;#theta_{#mu} [deg];Entries",
				 18, 0., 90.);
  TH1D *hist_muon_cos = new TH1D("hist_muon_cos",
				 "Muon reconstructed angle;cos#theta_{#mu};Entries",
				 20, 0., 1.);
  TH1D *hist_proton_ang = new TH1D("hist_proton_ang",
				   "Proton reconstructed angle;#theta_{p} [deg];Entries",
				   36, 0., 180.);
  TH1D *hist_proton_ang_high = new TH1D("hist_proton_ang_high",
					"Higher proton angle;#theta_{p1} [deg];Entries",
					36, 0., 180.);
  TH1D *hist_proton_ang_low = new TH1D("hist_proton_ang_low",
				       "Lower proton angle;#theta_{p2} [deg];Entries",
				       36, 0., 180.);
  TH1D *hist_proton_cos = new TH1D("hist_proton_cos",
				   "Proton reconstructed angle;cos#theta_{p};Entries",
				   40, -1., 1.);
  TH1D *hist_proton_cos_high = new TH1D("hist_proton_cos_high",
					"Higher proton angle;cos#theta_{p1};Entries",
					40, -1., 1.);
  TH1D *hist_proton_cos_low = new TH1D("hist_proton_cos_low",
				       "Lower proton angle;cos#theta_{p2};Entries",
				       40, -1., 1.);
  TH2D *hist_muon_mom_ang = new TH2D("hist_muon_mom_ang",
				     "Muon momentum vs angle;p_{#mu} [MeV/c];#theta_{#mu} [deg]",
				     15, 0., 1500., 18, 0., 90.);
  TH2D *hist_muon_mom_cos = new TH2D("hist_muon_mom_cos",
				     "Muon momentum vs angle;p_{#mu} [MeV/c];cos#theta_{#mu}",
				     15, 0., 1500., 20, 0., 1.);
  TH2D *hist_proton_mom_ang = new TH2D("hist_proton_mom_ang",
				       "Proton momentum vs angle;p_{p} [MeV/c];#theta_{p} [deg]",
				       15, 0., 1500., 36, 0., 180.);
  TH2D *hist_proton_mom_ang_high = new TH2D("hist_proton_mom_ang_high",
					    "Higher proton momentum vs angle;p_{p1} [MeV/c];#theta_{p1} [deg]",
					    15, 0., 1500., 36, 0., 180.);
  TH2D *hist_proton_mom_ang_low = new TH2D("hist_proton_mom_ang_low",
					    "Lower proton momentum vs angle;p_{p2} [MeV/c];#theta_{p2} [deg]",
					    15, 0., 1500., 36, 0., 180.);
  TH2D *hist_proton_mom_cos = new TH2D("hist_proton_mom_cos",
				       "Proton momentum vs angle;p_{p} [MeV/c];cos#theta_{p}",
				       15, 0., 1500., 40, -1., 1.);
  TH2D *hist_proton_mom_cos_high = new TH2D("hist_proton_mom_cos_high",
					    "Higher proton momentum vs angle;p_{p1} [MeV/c];cos#theta_{p1}",
					    15, 0., 1500., 40, -1., 1.);
  TH2D *hist_proton_mom_cos_low = new TH2D("hist_proton_mom_cos_low",
					    "Lower proton momentum vs angle;p_{p2} [MeV/c];cos#theta_{p2}",
					   15, 0., 1500., 40, -1., 1.);

  TH1D *hist_open_ang = new TH1D("hist_open_ang",
				 "Proton-proton opening angle;#theta_{pp} [deg];Entries",
				 36, 0., 180.);
  TH1D *hist_open_cos = new TH1D("hist_open_cos",
				 "Proton-proton opening angle;cos#theta_{pp};Entries",
				 40, -1., 1.);
  TH1D *hist_mom_ratio = new TH1D("hist_mom_ratio",
				  "Proton-proton momentum ratio;p_{p2}/p_{p1};Entries",
				  10, 0., 0.5);
  TH1D *hist_mom_vecsum = new TH1D("hist_mom_vecsum",
				   "Proton momentum vector sum;|#vec{p}_{p1} + #vec{p}_{p2}| [MeV/c];Entries",
				   100, 0., 2000.);
  TH1D *hist_mom_scasum = new TH1D("hist_mom_scasum",
				   "Proton momentum scalar sum;p_{p1} + p_{p2} [MeV/c];Entries",
				   100, 0., 2000.);
  TH2D *hist_mom_high_low = new TH2D("hist_mom_high_low",
				     "Proton momentum correlation;p_{p1} [MeV/c];p_{p2} [MeV/c]",
				     100, 0., 2000., 100, 0., 2000.);
  TH2D *hist_ang_high_low = new TH2D("hist_ang_high_low",
				     "Proton angle correlation;#theta_{p1} [deg];#theta_{p2} [deg]",
				     36, 0., 180., 36, 0., 180.);
  TH2D *hist_cos_high_low = new TH2D("hist_cos_high_low",
				     "Proton angle correlation;cos#theta_{p1};cos#theat_{p2}",
				     40, -1., 1., 40, -1., 1.);

  // TKI
  TH1D *hist_dptt = new TH1D("hist_dptt",
			     "#deltap_{TT};#deltap_{TT} [MeV/c];Entries",
			     100, -500., 500.);
  TH1D *hist_dpt = new TH1D("hist_dpt",
			    "#deltap_{T};#deltap_{T} [MeV/c];Entries",
			    200, 0., 2000.);
  TH1D *hist_pn = new TH1D("hist_pn",
			   "p_{N};p_{N} [MeV/c];Entries",
			   200, 0., 2000.);
  TH1D *hist_dalphat = new TH1D("hist_dalphat",
				"#delat#alpha_{T};#delta#alpha_{T} [deg];Entries",
				36, 0., 180.);
  TH1D *hist_cosdat = new TH1D("hist_cosdat",
			       "cos#delta#alpha_{T};cos#delta#alpha_{T};Entries",
			       40, -1., 1.);


  // mode
  TH1D *hist_mode_muon_mom[num_ninja_mode];
  TH1D *hist_mode_proton_mom[num_ninja_mode];
  TH1D *hist_mode_proton_mom_high[num_ninja_mode];
  TH1D *hist_mode_proton_mom_low[num_ninja_mode];
  TH1D *hist_mode_muon_ang[num_ninja_mode];
  TH1D *hist_mode_muon_cos[num_ninja_mode];
  TH1D *hist_mode_proton_ang[num_ninja_mode];
  TH1D *hist_mode_proton_ang_high[num_ninja_mode];
  TH1D *hist_mode_proton_ang_low[num_ninja_mode];
  TH1D *hist_mode_proton_cos[num_ninja_mode];
  TH1D *hist_mode_proton_cos_high[num_ninja_mode];
  TH1D *hist_mode_proton_cos_low[num_ninja_mode];
  TH1D *hist_mode_open_ang[num_ninja_mode];
  TH1D *hist_mode_open_cos[num_ninja_mode];
  TH1D *hist_mode_mom_ratio[num_ninja_mode];
  TH1D *hist_mode_mom_vecsum[num_ninja_mode];
  TH1D *hist_mode_mom_scasum[num_ninja_mode];
  TH1D *hist_mode_dptt[num_ninja_mode];
  TH1D *hist_mode_dpt[num_ninja_mode];
  TH1D *hist_mode_pn[num_ninja_mode];
  TH1D *hist_mode_dalphat[num_ninja_mode];
  TH1D *hist_mode_cosdat[num_ninja_mode];
  for ( int i = 0; i < num_ninja_mode; i++ ) {
    hist_mode_muon_mom[i] = new TH1D(Form("hist_muon_mom_%d", i), "", 15, 0., 1500.);
    hist_mode_muon_mom[i]->SetFillColor(mode_color[i]);
    hist_mode_muon_mom[i]->SetFillStyle(mode_style[i]);
    hist_mode_proton_mom[i] = new TH1D(Form("hist_proton_mom_%d", i), "", 15, 0., 1500.);
    hist_mode_proton_mom[i]->SetFillColor(mode_color[i]);
    hist_mode_proton_mom[i]->SetFillStyle(mode_style[i]);
    hist_mode_proton_mom_high[i] = new TH1D(Form("hist_proton_mom_high_%d", i), "", 15, 0., 1500.);
    hist_mode_proton_mom_high[i]->SetFillColor(mode_color[i]);
    hist_mode_proton_mom_high[i]->SetFillStyle(mode_style[i]);
    hist_mode_proton_mom_low[i] = new TH1D(Form("hist_proton_mom_low_%d", i), "", 15, 0., 1500.);
    hist_mode_proton_mom_low[i]->SetFillColor(mode_color[i]);
    hist_mode_proton_mom_low[i]->SetFillStyle(mode_style[i]);
    hist_mode_muon_ang[i] = new TH1D(Form("hist_muon_ang_%d", i), "", 18, 0., 90.);
    hist_mode_muon_ang[i]->SetFillColor(mode_color[i]);
    hist_mode_muon_ang[i]->SetFillStyle(mode_style[i]);
    hist_mode_proton_ang[i] = new TH1D(Form("hist_proton_ang_%d", i), "", 36, 0., 180.);
    hist_mode_proton_ang[i]->SetFillColor(mode_color[i]);
    hist_mode_proton_ang[i]->SetFillStyle(mode_style[i]);
    hist_mode_proton_ang_high[i] = new TH1D(Form("hist_proton_ang_high_%d", i), "", 36, 0., 180.);
    hist_mode_proton_ang_high[i]->SetFillColor(mode_color[i]);
    hist_mode_proton_ang_high[i]->SetFillStyle(mode_style[i]);
    hist_mode_proton_ang_low[i] = new TH1D(Form("hist_proton_ang_low_%d", i), "", 36, 0., 180.);
    hist_mode_proton_ang_low[i]->SetFillColor(mode_color[i]);
    hist_mode_proton_ang_low[i]->SetFillStyle(mode_style[i]);
    hist_mode_muon_cos[i] = new TH1D(Form("hist_muon_cos_%d", i), "", 20, 0., 1.);
    hist_mode_muon_cos[i]->SetFillColor(mode_color[i]);
    hist_mode_muon_cos[i]->SetFillStyle(mode_style[i]);
    hist_mode_proton_cos[i] = new TH1D(Form("hist_proton_cos_%d", i), "", 40, -1., 1.);
    hist_mode_proton_cos[i]->SetFillColor(mode_color[i]);
    hist_mode_proton_cos[i]->SetFillStyle(mode_style[i]);
    hist_mode_proton_cos_high[i] = new TH1D(Form("hist_proton_cos_high_%d", i), "", 40, -1., 1.);
    hist_mode_proton_cos_high[i]->SetFillColor(mode_color[i]);
    hist_mode_proton_cos_high[i]->SetFillStyle(mode_style[i]);
    hist_mode_proton_cos_low[i] = new TH1D(Form("hist_proton_cos_low_%d", i), "", 40, -1., 1.);
    hist_mode_proton_cos_low[i]->SetFillColor(mode_color[i]);
    hist_mode_proton_cos_low[i]->SetFillStyle(mode_style[i]);
    hist_mode_open_ang[i] = new TH1D(Form("hist_open_ang_%d", i), "", 36, 0., 180.);
    hist_mode_open_ang[i]->SetFillColor(mode_color[i]);
    hist_mode_open_ang[i]->SetFillStyle(mode_style[i]);
    hist_mode_open_cos[i] = new TH1D(Form("hist_open_cos_%d", i), "", 40, -1., 1.);
    hist_mode_open_cos[i]->SetFillColor(mode_color[i]);
    hist_mode_open_cos[i]->SetFillStyle(mode_style[i]);
    hist_mode_mom_ratio[i] = new TH1D(Form("hist_mom_ratio_%d", i), "", 10, 0., 0.5);
    hist_mode_mom_ratio[i]->SetFillColor(mode_color[i]);
    hist_mode_mom_ratio[i]->SetFillStyle(mode_style[i]);
    hist_mode_mom_vecsum[i] = new TH1D(Form("hist_mom_vecsum_%d", i), "", 100, 0., 2000.);
    hist_mode_mom_vecsum[i]->SetFillColor(mode_color[i]);
    hist_mode_mom_vecsum[i]->SetFillStyle(mode_style[i]);
    hist_mode_mom_scasum[i] = new TH1D(Form("hist_mom_scasum_%d", i), "", 100, 0., 2000.);
    hist_mode_mom_scasum[i]->SetFillColor(mode_color[i]);
    hist_mode_mom_scasum[i]->SetFillStyle(mode_style[i]);
    hist_mode_dptt[i] = new TH1D(Form("hist_dptt_%d", i), "", 100, -500., 500.);
    hist_mode_dptt[i]->SetFillColor(mode_color[i]);
    hist_mode_dptt[i]->SetFillStyle(mode_style[i]);
    hist_mode_dpt[i] = new TH1D(Form("hist_dpt_%d", i), "", 200, 0., 2000.);
    hist_mode_dpt[i]->SetFillColor(mode_color[i]);
    hist_mode_dpt[i]->SetFillStyle(mode_style[i]);
    hist_mode_pn[i] = new TH1D(Form("hist_pn_%d", i), "", 100, 0., 2000.);
    hist_mode_pn[i]->SetFillColor(mode_color[i]);
    hist_mode_pn[i]->SetFillStyle(mode_style[i]);
    hist_mode_dalphat[i] = new TH1D(Form("hist_dalphat_%d", i), "", 36, 0., 180.);
    hist_mode_dalphat[i]->SetFillColor(mode_color[i]);
    hist_mode_dalphat[i]->SetFillStyle(mode_style[i]);
    hist_mode_cosdat[i] = new TH1D(Form("hist_cosdat_%d", i), "", 40, -1., 1.);
    hist_mode_cosdat[i]->SetFillColor(mode_color[i]);
    hist_mode_cosdat[i]->SetFillStyle(mode_style[i]);
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
      
      TVector3 muon_tangent;
      TVector3 proton_tangent_tmp;
      TVector3 proton_tangent_high;
      TVector3 proton_tangent_low;
      int proton_direction_tmp;
      int proton_direction_high;
      int proton_direction_low;
      double muon_momentum;
      double proton_momentum_tmp;
      double proton_momentum_high;
      double proton_momentum_low;
      TVector3 muon_momentum_vec;
      TVector3 proton_momentum_vec_tmp;
      TVector3 proton_momentum_vec_high;
      TVector3 proton_momentum_vec_low;

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

	if ( particle_id == 13 ) {
	  num_muon++;
	  muon_tangent.SetXYZ(ax, ay, 1.);
	  muon_momentum = chain.ecc_mcs_mom[0]; // Baby MIND range?
	  muon_momentum_vec = (muon_momentum / muon_tangent.Mag()) * muon_tangent;
	}
	else if ( particle_id == 2212 ) {
	  num_proton++;
	  if ( chain.direction == 1 ) {
	    proton_tangent_tmp.SetXYZ(ax, ay, 1.);
	    proton_direction_tmp = 1;
	  }
	  else if ( chain.direction == -1 ) {
	    proton_tangent_tmp.SetXYZ(ax, ay, 1);
	    proton_direction_tmp = -1;
	  }

	  if ( chain.stop_flag == 2 ) {
	    proton_momentum_tmp = chain.ecc_range_mom[1];
	  }
	  else {
	    proton_momentum_tmp = chain.ecc_mcs_mom[1];
	  }
	  proton_momentum_vec_tmp = (proton_momentum_tmp * chain.direction / proton_tangent_tmp.Mag()) * proton_tangent_tmp;

	  if ( num_proton == 1 ) {
	    proton_tangent_high = proton_tangent_tmp;
	    proton_direction_high = proton_direction_tmp;
	    proton_momentum_high = proton_momentum_tmp;
	    proton_momentum_vec_high = proton_momentum_vec_tmp;
	  }
	  else if ( num_proton == 2 ) {
	    proton_tangent_low = proton_tangent_tmp;
	    proton_direction_low = proton_direction_tmp;
	    proton_momentum_low = proton_momentum_tmp;
	    proton_momentum_vec_low = proton_momentum_vec_tmp;
	  }
	}
	else if ( particle_id == 211 ) {
	  num_pion++;
	}

      }

      if ( num_muon != 1 || num_proton != 2 || num_pion != 0 ) continue;
      if ( proton_momentum_low < 100 || proton_momentum_high < 100 )
	std::cout << "Momentum error :"  << ev.groupid << std::endl;
      if ( proton_momentum_low > proton_momentum_high ) {
	std::swap(proton_tangent_high, proton_tangent_low);
	std::swap(proton_direction_high, proton_direction_low);
	std::swap(proton_momentum_high, proton_momentum_low);
	std::swap(proton_momentum_vec_high, proton_momentum_vec_low);
      }

      hist_muon_mom->Fill(muon_momentum, ev.weight);
      hist_mode_muon_mom[mode_id]->Fill(muon_momentum, ev.weight);
      hist_proton_mom->Fill(proton_momentum_high, ev.weight);
      hist_proton_mom->Fill(proton_momentum_low, ev.weight);
      hist_mode_proton_mom[mode_id]->Fill(proton_momentum_high, ev.weight);
      hist_mode_proton_mom[mode_id]->Fill(proton_momentum_low, ev.weight);
      hist_proton_mom_high->Fill(proton_momentum_high, ev.weight);
      hist_mode_proton_mom_high[mode_id]->Fill(proton_momentum_high, ev.weight);
      hist_proton_mom_low->Fill(proton_momentum_low, ev.weight);
      hist_mode_proton_mom_low[mode_id]->Fill(proton_momentum_low, ev.weight);

      double muon_rad = std::atan(std::hypot(muon_tangent.X(), muon_tangent.Y()));
      double muon_ang = muon_rad * TMath::RadToDeg();
      hist_muon_ang->Fill(muon_ang, ev.weight);
      hist_mode_muon_ang[mode_id]->Fill(muon_ang, ev.weight);
      hist_muon_cos->Fill(std::cos(muon_ang), ev.weight);
      hist_mode_muon_cos[mode_id]->Fill(std::cos(muon_ang), ev.weight);

      double proton_rad_high = std::atan(proton_direction_high * std::hypot(proton_tangent_high.X(),
									    proton_tangent_high.Y()));
      double proton_ang_high = proton_rad_high * TMath::RadToDeg();
      if ( proton_direction_high == -1 ) {
	proton_ang_high += 180.;
      }
      double proton_rad_low = std::atan(proton_direction_low * std::hypot(proton_tangent_low.X(),
									  proton_tangent_low.Y()));
      double proton_ang_low = proton_rad_low * TMath::RadToDeg();
      if ( proton_direction_low == -1 ) {
	proton_ang_low += 180.;
      }

      hist_proton_ang->Fill(proton_ang_high, ev.weight);
      hist_proton_ang->Fill(proton_ang_low, ev.weight);
      hist_mode_proton_ang[mode_id]->Fill(proton_ang_high, ev.weight);
      hist_mode_proton_ang[mode_id]->Fill(proton_ang_low, ev.weight);
      hist_proton_ang_high->Fill(proton_ang_high, ev.weight);
      hist_mode_proton_ang_high[mode_id]->Fill(proton_ang_high, ev.weight);
      hist_proton_ang_low->Fill(proton_ang_low, ev.weight);
      hist_mode_proton_ang_low[mode_id]->Fill(proton_ang_low, ev.weight);

      hist_proton_cos->Fill(std::cos(proton_ang_high * TMath::DegToRad()), ev.weight);
      hist_proton_cos->Fill(std::cos(proton_ang_low * TMath::DegToRad()), ev.weight);
      hist_mode_proton_cos[mode_id]->Fill(std::cos(proton_ang_high * TMath::DegToRad()), ev.weight);
      hist_mode_proton_cos[mode_id]->Fill(std::cos(proton_ang_low * TMath::DegToRad()), ev.weight);
      hist_proton_cos_high->Fill(std::cos(proton_ang_high * TMath::DegToRad()), ev.weight);
      hist_mode_proton_cos_high[mode_id]->Fill(std::cos(proton_ang_high * TMath::DegToRad()), ev.weight);
      hist_proton_cos_low->Fill(std::cos(proton_ang_low * TMath::DegToRad()), ev.weight);
      hist_mode_proton_cos_low[mode_id]->Fill(std::cos(proton_ang_low * TMath::DegToRad()), ev.weight);

      hist_muon_mom_ang->Fill(muon_momentum, muon_ang, ev.weight);
      hist_muon_mom_cos->Fill(muon_momentum, std::cos(muon_rad), ev.weight);
      hist_proton_mom_ang->Fill(proton_momentum_high, proton_ang_high, ev.weight);
      hist_proton_mom_ang->Fill(proton_momentum_low, proton_ang_low, ev.weight);
      hist_proton_mom_ang_high->Fill(proton_momentum_high, proton_ang_high, ev.weight);
      hist_proton_mom_ang_low->Fill(proton_momentum_low, proton_ang_low, ev.weight);
      hist_proton_mom_cos->Fill(proton_momentum_high, std::cos(proton_ang_high * TMath::DegToRad()), ev.weight);
      hist_proton_mom_cos->Fill(proton_momentum_low, std::cos(proton_ang_low * TMath::DegToRad()), ev.weight);
      hist_proton_mom_cos_high->Fill(proton_momentum_high, std::cos(proton_ang_high * TMath::DegToRad()), ev.weight);
      hist_proton_mom_cos_low->Fill(proton_momentum_low, std::cos(proton_ang_low * TMath::DegToRad()), ev.weight);

      double open_cos = (proton_momentum_vec_high * proton_momentum_vec_low)
	/ proton_momentum_vec_high.Mag() / proton_momentum_vec_low.Mag();
      double open_ang = std::acos(open_cos) * TMath::RadToDeg();
      hist_open_ang->Fill(open_ang, ev.weight);
      hist_mode_open_ang[mode_id]->Fill(open_ang, ev.weight);
      hist_open_cos->Fill(open_cos, ev.weight);
      hist_mode_open_cos[mode_id]->Fill(open_cos, ev.weight);
      hist_mom_ratio->Fill(proton_momentum_low / proton_momentum_high, ev.weight);
      hist_mode_mom_ratio[mode_id]->Fill(proton_momentum_low / proton_momentum_high, ev.weight);
      hist_mom_vecsum->Fill((proton_momentum_vec_high + proton_momentum_vec_low).Mag(), ev.weight);
      hist_mode_mom_vecsum[mode_id]->Fill((proton_momentum_vec_high + proton_momentum_vec_low).Mag(), ev.weight);
      hist_mom_scasum->Fill(proton_momentum_vec_high.Mag() + proton_momentum_vec_low.Mag(), ev.weight);
      hist_mode_mom_scasum[mode_id]->Fill(proton_momentum_vec_high.Mag() + proton_momentum_vec_low.Mag(), ev.weight);

      hist_mom_high_low->Fill(proton_momentum_high, proton_momentum_low, ev.weight);
      hist_ang_high_low->Fill(proton_ang_high, proton_ang_low, ev.weight);
      hist_cos_high_low->Fill(std::cos(proton_ang_high * TMath::DegToRad()),
			      std::cos(proton_ang_low * TMath::DegToRad()),
			      ev.weight);

      TVector3 ztt(-1. * muon_momentum_vec.Y(), muon_momentum_vec.X(), 0.);
      ztt = (1. / std::hypot(muon_momentum_vec.X(), muon_momentum_vec.Y())) * ztt;
      // ztt = pnu x pmu / |pnu x pmu|, where pnu = (0,0,Enu), pmu = muon_momentum_vec.
      double dptt = ztt * proton_momentum_vec_high
	+ ztt * proton_momentum_vec_low;
      hist_dptt->Fill(dptt, ev.weight);
      hist_mode_dptt[mode_id]->Fill(dptt, ev.weight);
      
      TVector2 mumom_vec_2d(muon_momentum_vec.X(), muon_momentum_vec.Y());
      TVector2 promom_vec_2d_high(proton_momentum_vec_high.X(), proton_momentum_vec_high.Y());
      TVector2 promom_vec_2d_low(proton_momentum_vec_low.X(), proton_momentum_vec_low.Y());
      TVector2 dpt_vec = mumom_vec_2d + promom_vec_2d_high + promom_vec_2d_low;
      hist_dpt->Fill(dpt_vec.Mod(), ev.weight);
      hist_mode_dpt[mode_id]->Fill(dpt_vec.Mod(), ev.weight);

      double muon_energy = std::sqrt(muon_momentum * muon_momentum + muon_mass * muon_mass);
      double proton_energy_high = std::sqrt(proton_momentum_high * proton_momentum_high
					    + proton_mass * proton_mass);
      double proton_energy_low = std::sqrt(proton_momentum_low * proton_momentum_low
					   + proton_mass * proton_mass);
      double residual_mass = 23.; // MeV?
      double pl = 0.5 * (proton_mass + muon_momentum_vec.Z() + proton_momentum_vec_high.Z() + proton_momentum_vec_low.Z()
			 - muon_energy - proton_energy_high - proton_energy_low)
	- 0.5 * (dpt_vec.Mod2() + residual_mass * residual_mass)
	/ (proton_mass + muon_momentum_vec.Z() + proton_momentum_vec_high.Z() + proton_momentum_vec_low.Z()
	   - muon_energy - proton_energy_high - proton_energy_low);
      double pn = std::sqrt(dpt_vec.Mod2() + pl * pl);
      hist_pn->Fill(pn, ev.weight);
      hist_mode_pn[mode_id]->Fill(pn, ev.weight);

      double dalphat = std::acos((-1. * mumom_vec_2d * dpt_vec) / mumom_vec_2d.Mod() / dpt_vec.Mod());
      hist_dalphat->Fill(dalphat * TMath::RadToDeg(), ev.weight);
      hist_mode_dalphat[mode_id]->Fill(dalphat * TMath::RadToDeg(), ev.weight);
      hist_cosdat->Fill((-1. * mumom_vec_2d * dpt_vec) / mumom_vec_2d.Mod() / dpt_vec.Mod(), ev.weight);
      hist_mode_cosdat[mode_id]->Fill((-1. * mumom_vec_2d * dpt_vec) / mumom_vec_2d.Mod() / dpt_vec.Mod(), ev.weight);
    }

  }

  outputfile->cd();
  hist_muon_mom->Write();
  hist_proton_mom->Write();
  hist_proton_mom_high->Write();
  hist_proton_mom_low->Write();
  hist_muon_ang->Write();
  hist_muon_cos->Write();
  hist_proton_ang->Write();
  hist_proton_ang_high->Write();
  hist_proton_ang_low->Write();
  hist_proton_cos->Write();
  hist_proton_cos_high->Write();
  hist_proton_cos_low->Write();
  hist_muon_mom_ang->Write();
  hist_muon_mom_cos->Write();
  hist_proton_mom_ang->Write();
  hist_proton_mom_ang_high->Write();
  hist_proton_mom_ang_low->Write();
  hist_proton_mom_cos->Write();
  hist_proton_mom_cos_high->Write();
  hist_proton_mom_cos_low->Write();
  hist_open_ang->Write();
  hist_open_cos->Write();
  hist_mom_ratio->Write();
  hist_mom_vecsum->Write();
  hist_mom_scasum->Write();
  hist_mom_high_low->Write();
  hist_ang_high_low->Write();
  hist_cos_high_low->Write();
  hist_dptt->Write();
  hist_dpt->Write();
  hist_pn->Write();
  hist_dalphat->Write();
  hist_cosdat->Write();
  for ( int i = 0; i < num_ninja_mode; i++ ) {
    hist_mode_muon_mom[i]->Write();
    hist_mode_proton_mom[i]->Write();
    hist_mode_proton_mom_high[i]->Write();
    hist_mode_proton_mom_low[i]->Write();
    hist_mode_muon_ang[i]->Write();
    hist_mode_proton_ang[i]->Write();
    hist_mode_proton_ang_high[i]->Write();
    hist_mode_proton_ang_low[i]->Write();
    hist_mode_muon_cos[i]->Write();
    hist_mode_proton_cos[i]->Write();
    hist_mode_proton_cos_high[i]->Write();
    hist_mode_proton_cos_low[i]->Write();
    hist_mode_open_ang[i]->Write();
    hist_mode_open_cos[i]->Write();
    hist_mode_mom_ratio[i]->Write();
    hist_mode_mom_vecsum[i]->Write();
    hist_mode_mom_scasum[i]->Write();
    hist_mode_dptt[i]->Write();
    hist_mode_dpt[i]->Write();
    hist_mode_pn[i]->Write();
    hist_mode_dalphat[i]->Write();
    hist_mode_cosdat[i]->Write();
  }
  outputfile->Close();

}
