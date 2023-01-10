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

  // output file
  TFile *outputfile = new TFile((TString)outputfilename, "recreate");
  BOOST_LOG_TRIVIAL(info) << "Output filename : " << outputfilename;

  // total
  TH1D *hist_muon_mom = new TH1D("hist_muon_mom",
				 "Muon reconstructed momentum;p_{#mu} [MeV/c];Entries",
				 muon_mom_bin_size-1, muon_mom_bins);
  TH1D *hist_proton_mom = new TH1D("hist_proton_mom",
				   "Proton reconstructed momentum;p_{p} [MeV/c];Entries",
				   hadron_mom_bin_size-1, hadron_mom_bins);
  TH1D *hist_proton_mom_high = new TH1D("hist_proton_mom_high",
					"Higher proton momentum;p_{p1} [MeV/c];Entries",
					hadron_mom_bin_size-1, hadron_mom_bins);
  TH1D *hist_proton_mom_low = new TH1D("hist_proton_mom_low",
				       "Lower proton momentum;p_{p2} [MeV/c];Entries",
				       hadron_mom_bin_size-1, hadron_mom_bins);
  TH1D *hist_muon_ang = new TH1D("hist_muon_ang",
				 "Muon reconstructed angle;#theta_{#mu} [deg];Entries",
				 muon_deg_bin_size-1, muon_deg_bins);
  TH1D *hist_muon_cos = new TH1D("hist_muon_cos",
				 "Muon reconstructed angle;cos#theta_{#mu};Entries",
				 muon_cos_bin_size-1, muon_cos_bins);
  TH1D *hist_proton_ang = new TH1D("hist_proton_ang",
				   "Proton reconstructed angle;#theta_{p} [deg];Entries",
				   hadron_deg_bin_size-1, hadron_deg_bins);
  TH1D *hist_proton_ang_high = new TH1D("hist_proton_ang_high",
					"Higher proton angle;#theta_{p1} [deg];Entries",
					hadron_deg_bin_size-1, hadron_deg_bins);
  TH1D *hist_proton_ang_low = new TH1D("hist_proton_ang_low",
				       "Lower proton angle;#theta_{p2} [deg];Entries",
				       hadron_deg_bin_size-1, hadron_deg_bins);
  TH1D *hist_proton_cos = new TH1D("hist_proton_cos",
				   "Proton reconstructed angle;cos#theta_{p};Entries",
				   hadron_cos_bin_size-1, hadron_cos_bins);
  TH1D *hist_proton_cos_high = new TH1D("hist_proton_cos_high",
					"Higher proton angle;cos#theta_{p1};Entries",
					hadron_cos_bin_size-1, hadron_cos_bins);
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

  TH1D *hist_muon_mom_mcs = new TH1D("hist_muon_mom_mcs", "", muon_mom_bin_size-1, muon_mom_bins);
  TH1D *hist_muon_mom_range = new TH1D("hist_muon_mom_range", "", muon_mom_bin_size-1, muon_mom_bins);
  TH1D *hist_proton_mom_mcs = new TH1D("hist_proton_mom_mcs", "", hadron_mom_bin_size-1, hadron_mom_bins);
  TH1D *hist_proton_mom_mcs_high = new TH1D("hist_proton_mom_mcs_high", "", hadron_mom_bin_size-1, hadron_mom_bins);
  TH1D *hist_proton_mom_mcs_low = new TH1D("hist_proton_mom_mcs_low", "", hadron_mom_bin_size-1, hadron_mom_bins);
  TH1D *hist_proton_mom_range = new TH1D("hist_proton_mom_range", "", hadron_mom_bin_size-1, hadron_mom_bins);
  TH1D *hist_proton_mom_range_high = new TH1D("hist_proton_mom_range_high", "", hadron_mom_bin_size-1, hadron_mom_bins);
  TH1D *hist_proton_mom_range_low = new TH1D("hist_proton_mom_range_low", "", hadron_mom_bin_size-1, hadron_mom_bins);
  TH2D *hist_muon_mom_ang_mcs = new TH2D("hist_muon_mom_ang_mcs", "", 200, 0., 2000., 90, 0., 90.);
  TH2D *hist_muon_mom_cos_mcs = new TH2D("hist_muon_mom_cos_mcs", "", 200, 0., 2000., 100, 0., 1.);
  TH2D *hist_muon_mom_ang_range = new TH2D("hist_muon_mom_ang_range", "", 200, 0., 2000., 90, 0., 90.);
  TH2D *hist_muon_mom_cos_range = new TH2D("hist_muon_cos_range", "", 200, 0., 2000., 100, 0., 1.);
  TH2D *hist_proton_mom_ang_mcs = new TH2D("hist_proton_mom_ang_mcs", "", 150, 0., 1500., 180, 0., 180.);
  TH2D *hist_proton_mom_ang_mcs_high = new TH2D("hist_proton_mom_ang_mcs_high", "", 150, 0., 1500., 180, 0., 180.);
  TH2D *hist_proton_mom_ang_mcs_low = new TH2D("hist_proton_mom_ang_mcs_low", "", 150, 0., 1500., 180, 0., 180.);
  TH2D *hist_proton_mom_cos_mcs = new TH2D("hist_proton_mom_cos_mcs", "", 150, 0., 1500., 200, -1., 1.);
  TH2D *hist_proton_mom_cos_mcs_high = new TH2D("hist_proton_mom_cos_mcs_high", "", 150, 0., 1500., 200, -1., 1.);
  TH2D *hist_proton_mom_cos_mcs_low = new TH2D("hist_proton_mom_cos_mcs_low", "", 150, 0., 1500., 200, -1., 1.);
  TH2D *hist_proton_mom_ang_range = new TH2D("hist_proton_mom_ang_range", "", 150, 0., 1500., 180, 0., 180.);
  TH2D *hist_proton_mom_ang_range_high = new TH2D("hist_proton_mom_ang_range_high", "", 150, 0., 1500., 180, 0., 180.);
  TH2D *hist_proton_mom_ang_range_low = new TH2D("hist_proton_mom_ang_range_low", "", 150, 0., 1500., 180, 0., 180.);
  TH2D *hist_proton_mom_cos_range = new TH2D("hist_proton_mom_cos_range", "", 150, 0., 1500., 200, -1., 1.);
  TH2D *hist_proton_mom_cos_range_high = new TH2D("hist_proton_mom_cos_range_high", "", 150, 0., 1500., 200, -1., 1.);
  TH2D *hist_proton_mom_cos_range_low = new TH2D("hist_proton_mom_cos_range_low", "", 150, 0., 1500., 200, -1., 1.);

  TH1D *hist_open_ang = new TH1D("hist_open_ang",
				 "Proton-proton opening angle;#theta_{pp} [deg];Entries",
				 open_deg_bin_size-1, open_deg_bins);
  TH1D *hist_open_cos = new TH1D("hist_open_cos",
				 "Proton-proton opening angle;cos#theta_{pp};Entries",
				 open_cos_bin_size-1, open_cos_bins);
  TH1D *hist_mom_ratio = new TH1D("hist_mom_ratio",
				  "Proton-proton momentum ratio;p_{p2}/p_{p1};Entries",
				  mom_ratio_bin_size-1, mom_ratio_bins);
  TH1D *hist_mom_vecsum = new TH1D("hist_mom_vecsum",
				   "Proton momentum vector sum;|#vec{p}_{p1} + #vec{p}_{p2}| [MeV/c];Entries",
				   mom_vecsum_bin_size-1, mom_vecsum_bins);
  TH1D *hist_mom_scasum = new TH1D("hist_mom_scasum",
				   "Proton momentum scalar sum;p_{p1} + p_{p2} [MeV/c];Entries",
				   mom_scasum_bin_size-1, mom_scasum_bins);
  TH2D *hist_mom_high_low = new TH2D("hist_mom_high_low",
				     "Proton momentum correlation;p_{p1} [MeV/c];p_{p2} [MeV/c]",
				     150, 0., 1500., 150, 0., 1500.);
  TH2D *hist_ang_high_low = new TH2D("hist_ang_high_low",
				     "Proton angle correlation;#theta_{p1} [deg];#theta_{p2} [deg]",
				     180, 0., 180., 180, 0., 180.);
  TH2D *hist_cos_high_low = new TH2D("hist_cos_high_low",
				     "Proton angle correlation;cos#theta_{p1};cos#theat_{p2}",
				     200, -1., 1., 200, -1., 1.);

  // TKI
  TH1D *hist_dptt = new TH1D("hist_dptt",
			     "#deltap_{TT};#deltap_{TT} [MeV/c];Entries",
			     dptt_bin_size-1, dptt_bins);
  TH1D *hist_dpt = new TH1D("hist_dpt",
			    "#deltap_{T};#deltap_{T} [MeV/c];Entries",
			    dpt_2p_bin_size-1, dpt_2p_bins);
  TH1D *hist_pn = new TH1D("hist_pn",
			   "p_{N};p_{N} [MeV/c];Entries",
			   pn_bin_size-1, pn_bins);
  TH1D *hist_dalphat = new TH1D("hist_dalphat",
				"#delat#alpha_{T};#delta#alpha_{T} [deg];Entries",
				dalphat_2p_bin_size-1, dalphat_2p_bins);
  TH1D *hist_cosdat = new TH1D("hist_cosdat",
			       "cos#delta#alpha_{T};cos#delta#alpha_{T};Entries",
			       cosdat_2p_bin_size-1, cosdat_2p_bins);

  TH1D *hist_dptt_range = new TH1D("hist_dptt_range", "", dptt_bin_size-1, dptt_bins);
  TH1D *hist_dpt_range = new TH1D("hist_dpt_range", "", dpt_2p_bin_size-1, dpt_2p_bins);
  TH1D *hist_pn_range = new TH1D("hist_pn_range", "", pn_bin_size-1, pn_bins);
  TH1D *hist_dalphat_range = new TH1D("hist_dalphat_range", "", dalphat_2p_bin_size-1, dalphat_2p_bins);
  TH1D *hist_cosdat_range = new TH1D("hist_cosdat_range", "", cosdat_2p_bin_size-1, cosdat_2p_bins);

  // mode
  TH1D *hist_mode_muon_mom[num_ninja_mode];
  TH1D *hist_mode_proton_mom[num_ninja_mode];
  TH1D *hist_mode_proton_mom_high[num_ninja_mode];
  TH1D *hist_mode_proton_mom_low[num_ninja_mode];
  TH1D *hist_mode_muon_mom_mcs[num_ninja_mode];
  TH1D *hist_mode_proton_mom_mcs[num_ninja_mode];
  TH1D *hist_mode_proton_mom_mcs_high[num_ninja_mode];
  TH1D *hist_mode_proton_mom_mcs_low[num_ninja_mode];
  TH1D *hist_mode_muon_mom_range[num_ninja_mode];
  TH1D *hist_mode_proton_mom_range[num_ninja_mode];
  TH1D *hist_mode_proton_mom_range_high[num_ninja_mode];
  TH1D *hist_mode_proton_mom_range_low[num_ninja_mode];
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
  TH1D *hist_mode_dptt_range[num_ninja_mode];
  TH1D *hist_mode_dpt_range[num_ninja_mode];
  TH1D *hist_mode_pn_range[num_ninja_mode];
  TH1D *hist_mode_dalphat_range[num_ninja_mode];
  TH1D *hist_mode_cosdat_range[num_ninja_mode];
  for ( int i = 0; i < num_ninja_mode; i++ ) {
    hist_mode_muon_mom[i] = new TH1D(Form("hist_muon_mom_%d", i), "", muon_mom_bin_size-1, muon_mom_bins);
    hist_mode_muon_mom[i]->SetFillColor(mode_color[i]);
    hist_mode_muon_mom[i]->SetFillStyle(mode_style[i]);
    hist_mode_proton_mom[i] = new TH1D(Form("hist_proton_mom_%d", i), "", hadron_mom_bin_size-1, hadron_mom_bins);
    hist_mode_proton_mom[i]->SetFillColor(mode_color[i]);
    hist_mode_proton_mom[i]->SetFillStyle(mode_style[i]);
    hist_mode_proton_mom_high[i] = new TH1D(Form("hist_proton_mom_high_%d", i), "", hadron_mom_bin_size-1, hadron_mom_bins);
    hist_mode_proton_mom_high[i]->SetFillColor(mode_color[i]);
    hist_mode_proton_mom_high[i]->SetFillStyle(mode_style[i]);
    hist_mode_proton_mom_low[i] = new TH1D(Form("hist_proton_mom_low_%d", i), "", hadron_mom_bin_size-1, hadron_mom_bins);
    hist_mode_proton_mom_low[i]->SetFillColor(mode_color[i]);
    hist_mode_proton_mom_low[i]->SetFillStyle(mode_style[i]);
    hist_mode_muon_mom_mcs[i] = new TH1D(Form("hist_muon_mom_mcs_%d", i), "", muon_mom_bin_size-1, muon_mom_bins);
    hist_mode_muon_mom_mcs[i]->SetFillColor(mode_color[i]);
    hist_mode_muon_mom_mcs[i]->SetFillStyle(mode_style[i]);
    hist_mode_proton_mom_mcs[i] = new TH1D(Form("hist_proton_mom_mcs_%d", i), "", hadron_mom_bin_size-1, hadron_mom_bins);
    hist_mode_proton_mom_mcs[i]->SetFillColor(mode_color[i]);
    hist_mode_proton_mom_mcs[i]->SetFillStyle(mode_style[i]);
    hist_mode_proton_mom_mcs_high[i] = new TH1D(Form("hist_proton_mom_mcs_high_%d", i), "", hadron_mom_bin_size-1, hadron_mom_bins);
    hist_mode_proton_mom_mcs_high[i]->SetFillColor(mode_color[i]);
    hist_mode_proton_mom_mcs_high[i]->SetFillStyle(mode_style[i]);
    hist_mode_proton_mom_mcs_low[i] = new TH1D(Form("hist_proton_mom_mcs_low_%d", i), "", hadron_mom_bin_size-1, hadron_mom_bins);
    hist_mode_proton_mom_mcs_low[i]->SetFillColor(mode_color[i]);
    hist_mode_proton_mom_mcs_low[i]->SetFillStyle(mode_style[i]);
    hist_mode_muon_mom_range[i] = new TH1D(Form("hist_muon_mom_range_%d", i), "", muon_mom_bin_size-1, muon_mom_bins);
    hist_mode_muon_mom_range[i]->SetFillColor(mode_color[i]);
    hist_mode_muon_mom_range[i]->SetFillStyle(mode_style[i]);
    hist_mode_proton_mom_range[i] = new TH1D(Form("hist_proton_mom_range_%d", i), "", hadron_mom_bin_size-1, hadron_mom_bins);
    hist_mode_proton_mom_range[i]->SetFillColor(mode_color[i]);
    hist_mode_proton_mom_range[i]->SetFillStyle(mode_style[i]);
    hist_mode_proton_mom_range_high[i] = new TH1D(Form("hist_proton_mom_range_high_%d", i), "", hadron_mom_bin_size-1, hadron_mom_bins);
    hist_mode_proton_mom_range_high[i]->SetFillColor(mode_color[i]);
    hist_mode_proton_mom_range_high[i]->SetFillStyle(mode_style[i]);
    hist_mode_proton_mom_range_low[i] = new TH1D(Form("hist_proton_mom_range_low_%d", i), "", hadron_mom_bin_size-1, hadron_mom_bins);
    hist_mode_proton_mom_range_low[i]->SetFillColor(mode_color[i]);
    hist_mode_proton_mom_range_low[i]->SetFillStyle(mode_style[i]);
    hist_mode_muon_ang[i] = new TH1D(Form("hist_muon_ang_%d", i), "", muon_deg_bin_size-1, muon_deg_bins);
    hist_mode_muon_ang[i]->SetFillColor(mode_color[i]);
    hist_mode_muon_ang[i]->SetFillStyle(mode_style[i]);
    hist_mode_proton_ang[i] = new TH1D(Form("hist_proton_ang_%d", i), "", hadron_deg_bin_size-1, hadron_deg_bins);
    hist_mode_proton_ang[i]->SetFillColor(mode_color[i]);
    hist_mode_proton_ang[i]->SetFillStyle(mode_style[i]);
    hist_mode_proton_ang_high[i] = new TH1D(Form("hist_proton_ang_high_%d", i), "", hadron_deg_bin_size-1, hadron_deg_bins);
    hist_mode_proton_ang_high[i]->SetFillColor(mode_color[i]);
    hist_mode_proton_ang_high[i]->SetFillStyle(mode_style[i]);
    hist_mode_proton_ang_low[i] = new TH1D(Form("hist_proton_ang_low_%d", i), "", hadron_deg_bin_size-1, hadron_deg_bins);
    hist_mode_proton_ang_low[i]->SetFillColor(mode_color[i]);
    hist_mode_proton_ang_low[i]->SetFillStyle(mode_style[i]);
    hist_mode_muon_cos[i] = new TH1D(Form("hist_muon_cos_%d", i), "", muon_cos_bin_size-1, muon_cos_bins);
    hist_mode_muon_cos[i]->SetFillColor(mode_color[i]);
    hist_mode_muon_cos[i]->SetFillStyle(mode_style[i]);
    hist_mode_proton_cos[i] = new TH1D(Form("hist_proton_cos_%d", i), "", hadron_cos_bin_size-1, hadron_cos_bins);
    hist_mode_proton_cos[i]->SetFillColor(mode_color[i]);
    hist_mode_proton_cos[i]->SetFillStyle(mode_style[i]);
    hist_mode_proton_cos_high[i] = new TH1D(Form("hist_proton_cos_high_%d", i), "", hadron_cos_bin_size-1, hadron_cos_bins);
    hist_mode_proton_cos_high[i]->SetFillColor(mode_color[i]);
    hist_mode_proton_cos_high[i]->SetFillStyle(mode_style[i]);
    hist_mode_proton_cos_low[i] = new TH1D(Form("hist_proton_cos_low_%d", i), "", hadron_cos_bin_size-1, hadron_cos_bins);
    hist_mode_proton_cos_low[i]->SetFillColor(mode_color[i]);
    hist_mode_proton_cos_low[i]->SetFillStyle(mode_style[i]);
    hist_mode_open_ang[i] = new TH1D(Form("hist_open_ang_%d", i), "", open_deg_bin_size-1, open_deg_bins);
    hist_mode_open_ang[i]->SetFillColor(mode_color[i]);
    hist_mode_open_ang[i]->SetFillStyle(mode_style[i]);
    hist_mode_open_cos[i] = new TH1D(Form("hist_open_cos_%d", i), "", open_cos_bin_size-1, open_cos_bins);
    hist_mode_open_cos[i]->SetFillColor(mode_color[i]);
    hist_mode_open_cos[i]->SetFillStyle(mode_style[i]);
    hist_mode_mom_ratio[i] = new TH1D(Form("hist_mom_ratio_%d", i), "", mom_ratio_bin_size-1, mom_ratio_bins);
    hist_mode_mom_ratio[i]->SetFillColor(mode_color[i]);
    hist_mode_mom_ratio[i]->SetFillStyle(mode_style[i]);
    hist_mode_mom_vecsum[i] = new TH1D(Form("hist_mom_vecsum_%d", i), "", mom_vecsum_bin_size-1, mom_vecsum_bins);
    hist_mode_mom_vecsum[i]->SetFillColor(mode_color[i]);
    hist_mode_mom_vecsum[i]->SetFillStyle(mode_style[i]);
    hist_mode_mom_scasum[i] = new TH1D(Form("hist_mom_scasum_%d", i), "", mom_scasum_bin_size-1, mom_scasum_bins);
    hist_mode_mom_scasum[i]->SetFillColor(mode_color[i]);
    hist_mode_mom_scasum[i]->SetFillStyle(mode_style[i]);
    hist_mode_dptt[i] = new TH1D(Form("hist_dptt_%d", i), "", dptt_bin_size-1, dptt_bins);
    hist_mode_dptt[i]->SetFillColor(mode_color[i]);
    hist_mode_dptt[i]->SetFillStyle(mode_style[i]);
    hist_mode_dpt[i] = new TH1D(Form("hist_dpt_%d", i), "", dpt_2p_bin_size-1, dpt_2p_bins);
    hist_mode_dpt[i]->SetFillColor(mode_color[i]);
    hist_mode_dpt[i]->SetFillStyle(mode_style[i]);
    hist_mode_pn[i] = new TH1D(Form("hist_pn_%d", i), "", pn_bin_size-1, pn_bins);
    hist_mode_pn[i]->SetFillColor(mode_color[i]);
    hist_mode_pn[i]->SetFillStyle(mode_style[i]);
    hist_mode_dalphat[i] = new TH1D(Form("hist_dalphat_%d", i), "", dalphat_2p_bin_size-1, dalphat_2p_bins);
    hist_mode_dalphat[i]->SetFillColor(mode_color[i]);
    hist_mode_dalphat[i]->SetFillStyle(mode_style[i]);
    hist_mode_cosdat[i] = new TH1D(Form("hist_cosdat_%d", i), "", cosdat_2p_bin_size-1, cosdat_2p_bins);
    hist_mode_cosdat[i]->SetFillColor(mode_color[i]);
    hist_mode_cosdat[i]->SetFillStyle(mode_style[i]);
    hist_mode_dptt_range[i] = new TH1D(Form("hist_dptt_%d", i), "", dptt_bin_size-1, dptt_bins);
    hist_mode_dptt_range[i]->SetFillColor(mode_color[i]);
    hist_mode_dptt_range[i]->SetFillStyle(mode_style[i]);
    hist_mode_dpt_range[i] = new TH1D(Form("hist_dpt_%d", i), "", dpt_2p_bin_size-1, dpt_2p_bins);
    hist_mode_dpt_range[i]->SetFillColor(mode_color[i]);
    hist_mode_dpt_range[i]->SetFillStyle(mode_style[i]);
    hist_mode_pn_range[i] = new TH1D(Form("hist_pn_%d", i), "", pn_bin_size-1, pn_bins);
    hist_mode_pn_range[i]->SetFillColor(mode_color[i]);
    hist_mode_pn_range[i]->SetFillStyle(mode_style[i]);
    hist_mode_dalphat_range[i] = new TH1D(Form("hist_dalphat_%d", i), "", dalphat_2p_bin_size-1, dalphat_2p_bins);
    hist_mode_dalphat_range[i]->SetFillColor(mode_color[i]);
    hist_mode_dalphat_range[i]->SetFillStyle(mode_style[i]);
    hist_mode_cosdat_range[i] = new TH1D(Form("hist_cosdat_%d", i), "", cosdat_2p_bin_size-1, cosdat_2p_bins);
    hist_mode_cosdat_range[i]->SetFillColor(mode_color[i]);
    hist_mode_cosdat_range[i]->SetFillStyle(mode_style[i]);
  }

  // Muon mis-id
  TH1D *hist_muon_misid_muon_mom = new TH1D("hist_muon_misid_muon_mom", "", muon_mom_bin_size-1, muon_mom_bins);
  TH1D *hist_muon_misid_proton_mom = new TH1D("hist_muon_misid_proton_mom", "", hadron_mom_bin_size-1, hadron_mom_bins);
  TH1D *hist_muon_misid_proton_mom_high = new TH1D("hist_muon_misid_proton_mom_high", "", hadron_mom_bin_size-1, hadron_mom_bins);
  TH1D *hist_muon_misid_proton_mom_low = new TH1D("hist_muon_misid_proton_mom_low", "", hadron_mom_bin_size-1, hadron_mom_bins); 
  TH1D *hist_muon_misid_muon_ang = new TH1D("hist_muon_misid_muon_ang", "", muon_deg_bin_size-1, muon_deg_bins); 
  TH1D *hist_muon_misid_muon_cos = new TH1D("hist_muon_misid_muon_cos", "", muon_cos_bin_size-1, muon_cos_bins);
  TH1D *hist_muon_misid_proton_ang = new TH1D("hist_muon_misid_proton_ang", "", hadron_deg_bin_size-1, hadron_deg_bins);
  TH1D *hist_muon_misid_proton_ang_high = new TH1D("hist_muon_misid_proton_ang_high", "", hadron_deg_bin_size-1, hadron_deg_bins);
  TH1D *hist_muon_misid_proton_ang_low = new TH1D("hist_muon_misid_proton_ang_low", "", hadron_deg_bin_size-1, hadron_deg_bins); 
  TH1D *hist_muon_misid_proton_cos = new TH1D("hist_muon_misid_proton_cos", "", hadron_cos_bin_size-1, hadron_cos_bins);
  TH1D *hist_muon_misid_proton_cos_high = new TH1D("hist_muon_misid_proton_cos_high", "", hadron_cos_bin_size-1, hadron_cos_bins);
  TH1D *hist_muon_misid_proton_cos_low = new TH1D("hist_muon_misid_proton_cos_low", "", hadron_cos_bin_size-1, hadron_cos_bins);
  TH1D *hist_muon_misid_muon_mom_mcs = new TH1D("hist_muon_misid_muon_mom_mcs", "", muon_mom_bin_size-1, muon_mom_bins);
  TH1D *hist_muon_misid_muon_mom_range = new TH1D("hist_muon_misid_muon_mom_range", "", muon_mom_bin_size-1, muon_mom_bins);
  TH1D *hist_muon_misid_proton_mom_mcs = new TH1D("hist_muon_misid_proton_mom_mcs", "", hadron_mom_bin_size-1, hadron_mom_bins);
  TH1D *hist_muon_misid_proton_mom_mcs_high = new TH1D("hist_muon_misid_proton_mom_mcs_high", "", hadron_mom_bin_size-1, hadron_mom_bins);
  TH1D *hist_muon_misid_proton_mom_mcs_low = new TH1D("hist_muon_misid_proton_mom_mcs_low", "", hadron_mom_bin_size-1, hadron_mom_bins); 
  TH1D *hist_muon_misid_proton_mom_range = new TH1D("hist_muon_misid_proton_mom_range", "", hadron_mom_bin_size-1, hadron_mom_bins);
  TH1D *hist_muon_misid_proton_mom_range_high = new TH1D("hist_muon_misid_proton_mom_range_high", "", hadron_mom_bin_size-1, hadron_mom_bins);
  TH1D *hist_muon_misid_proton_mom_range_low = new TH1D("hist_muon_misid_proton_mom_range_low", "", hadron_mom_bin_size-1, hadron_mom_bins);

  TH1D *hist_muon_misid_open_ang = new TH1D("hist_muon_misid_open_ang", "", open_deg_bin_size-1, open_deg_bins);
  TH1D *hist_muon_misid_open_cos = new TH1D("hist_muon_misid_open_cos", "", open_cos_bin_size-1, open_cos_bins);
  TH1D *hist_muon_misid_mom_ratio = new TH1D("hist_muon_misid_mom_ratio", "", mom_ratio_bin_size-1, mom_ratio_bins);
  TH1D *hist_muon_misid_mom_vecsum = new TH1D("hist_muon_misid_mom_vecsum", "", mom_vecsum_bin_size-1, mom_vecsum_bins);
  TH1D *hist_muon_misid_mom_scasum = new TH1D("hist_muon_misid_mom_scasum", "", mom_scasum_bin_size-1, mom_scasum_bins);
  
  TH1D *hist_muon_misid_dptt = new TH1D("hist_muon_misid_dptt", "", dptt_bin_size-1, dptt_bins);
  TH1D *hist_muon_misid_dpt = new TH1D("hist_muon_misid_dpt", "", dpt_2p_bin_size-1, dpt_2p_bins);
  TH1D *hist_muon_misid_pn = new TH1D("hist_muon_misid_pn", "", pn_bin_size-1, pn_bins);
  TH1D *hist_muon_misid_dalphat = new TH1D("hist_muon_misid_dalphat", "", dalphat_2p_bin_size-1, dalphat_2p_bins);
  TH1D *hist_muon_misid_cosdat = new TH1D("hist_muon_misid_cosdat", "", cosdat_2p_bin_size-1, cosdat_2p_bins);
  TH1D *hist_muon_misid_dptt_range = new TH1D("hist_muon_misid_dptt_range", "", dptt_bin_size-1, dptt_bins);
  TH1D *hist_muon_misid_dpt_range = new TH1D("hist_muon_misid_dpt_range", "", dpt_2p_bin_size-1, dpt_2p_bins);
  TH1D *hist_muon_misid_pn_range = new TH1D("hist_muon_misid_pn_range", "", pn_bin_size-1, pn_bins);
  TH1D *hist_muon_misid_dalphat_range = new TH1D("hist_muon_misid_dalphat_range", "", dalphat_2p_bin_size-1, dalphat_2p_bins);
  TH1D *hist_muon_misid_cosdat_range = new TH1D("hist_muon_misid_cosdat_range", "", cosdat_2p_bin_size-1, cosdat_2p_bins);

  // Proton mis-id
  TH1D *hist_proton_misid_muon_mom = new TH1D("hist_proton_misid_muon_mom", "", muon_mom_bin_size-1, muon_mom_bins);
  TH1D *hist_proton_misid_proton_mom = new TH1D("hist_proton_misid_proton_mom", "", hadron_mom_bin_size-1, hadron_mom_bins);
  TH1D *hist_proton_misid_proton_mom_high = new TH1D("hist_proton_misid_proton_mom_high", "", hadron_mom_bin_size-1, hadron_mom_bins);
  TH1D *hist_proton_misid_proton_mom_low = new TH1D("hist_proton_misid_proton_mom_low", "", hadron_mom_bin_size-1, hadron_mom_bins); 
  TH1D *hist_proton_misid_muon_ang = new TH1D("hist_proton_misid_muon_ang", "", muon_deg_bin_size-1, muon_deg_bins); 
  TH1D *hist_proton_misid_muon_cos = new TH1D("hist_proton_misid_muon_cos", "", muon_cos_bin_size-1, muon_cos_bins);
  TH1D *hist_proton_misid_proton_ang = new TH1D("hist_proton_misid_proton_ang", "", hadron_deg_bin_size-1, hadron_deg_bins);
  TH1D *hist_proton_misid_proton_ang_high = new TH1D("hist_proton_misid_proton_ang_high", "", hadron_deg_bin_size-1, hadron_deg_bins);
  TH1D *hist_proton_misid_proton_ang_low = new TH1D("hist_proton_misid_proton_ang_low", "", hadron_deg_bin_size-1, hadron_deg_bins); 
  TH1D *hist_proton_misid_proton_cos = new TH1D("hist_proton_misid_proton_cos", "", hadron_cos_bin_size-1, hadron_cos_bins);
  TH1D *hist_proton_misid_proton_cos_high = new TH1D("hist_proton_misid_proton_cos_high", "", hadron_cos_bin_size-1, hadron_cos_bins);
  TH1D *hist_proton_misid_proton_cos_low = new TH1D("hist_proton_misid_proton_cos_low", "", hadron_cos_bin_size-1, hadron_cos_bins);
  TH1D *hist_proton_misid_muon_mom_mcs = new TH1D("hist_proton_misid_muon_mom_mcs", "", muon_mom_bin_size-1, muon_mom_bins);
  TH1D *hist_proton_misid_muon_mom_range = new TH1D("hist_proton_misid_muon_mom_range", "", muon_mom_bin_size-1, muon_mom_bins);
  TH1D *hist_proton_misid_proton_mom_mcs = new TH1D("hist_proton_misid_proton_mom_mcs", "", hadron_mom_bin_size-1, hadron_mom_bins);
  TH1D *hist_proton_misid_proton_mom_mcs_high = new TH1D("hist_proton_misid_proton_mom_mcs_high", "", hadron_mom_bin_size-1, hadron_mom_bins);
  TH1D *hist_proton_misid_proton_mom_mcs_low = new TH1D("hist_proton_misid_proton_mom_mcs_low", "", hadron_mom_bin_size-1, hadron_mom_bins); 
  TH1D *hist_proton_misid_proton_mom_range = new TH1D("hist_proton_misid_proton_mom_range", "", hadron_mom_bin_size-1, hadron_mom_bins);
  TH1D *hist_proton_misid_proton_mom_range_high = new TH1D("hist_proton_misid_proton_mom_range_high", "", hadron_mom_bin_size-1, hadron_mom_bins);
  TH1D *hist_proton_misid_proton_mom_range_low = new TH1D("hist_proton_misid_proton_mom_range_low", "", hadron_mom_bin_size-1, hadron_mom_bins);

  TH1D *hist_proton_misid_open_ang = new TH1D("hist_proton_misid_open_ang", "", open_deg_bin_size-1, open_deg_bins);
  TH1D *hist_proton_misid_open_cos = new TH1D("hist_proton_misid_open_cos", "", open_cos_bin_size-1, open_cos_bins);
  TH1D *hist_proton_misid_mom_ratio = new TH1D("hist_proton_misid_mom_ratio", "", mom_ratio_bin_size-1, mom_ratio_bins);
  TH1D *hist_proton_misid_mom_vecsum = new TH1D("hist_proton_misid_mom_vecsum", "", mom_vecsum_bin_size-1, mom_vecsum_bins);
  TH1D *hist_proton_misid_mom_scasum = new TH1D("hist_proton_misid_mom_scasum", "", mom_scasum_bin_size-1, mom_scasum_bins);
  
  TH1D *hist_proton_misid_dptt = new TH1D("hist_proton_misid_dptt", "", dptt_bin_size-1, dptt_bins);
  TH1D *hist_proton_misid_dpt = new TH1D("hist_proton_misid_dpt", "", dpt_2p_bin_size-1, dpt_2p_bins);
  TH1D *hist_proton_misid_pn = new TH1D("hist_proton_misid_pn", "", pn_bin_size-1, pn_bins);
  TH1D *hist_proton_misid_dalphat = new TH1D("hist_proton_misid_dalphat", "", dalphat_2p_bin_size-1, dalphat_2p_bins);
  TH1D *hist_proton_misid_cosdat = new TH1D("hist_proton_misid_cosdat", "", cosdat_2p_bin_size-1, cosdat_2p_bins);
  TH1D *hist_proton_misid_dptt_range = new TH1D("hist_proton_misid_dptt_range", "", dptt_bin_size-1, dptt_bins);
  TH1D *hist_proton_misid_dpt_range = new TH1D("hist_proton_misid_dpt_range", "", dpt_2p_bin_size-1, dpt_2p_bins);
  TH1D *hist_proton_misid_pn_range = new TH1D("hist_proton_misid_pn_range", "", pn_bin_size-1, pn_bins);
  TH1D *hist_proton_misid_dalphat_range = new TH1D("hist_proton_misid_dalphat_range", "", dalphat_2p_bin_size-1, dalphat_2p_bins);
  TH1D *hist_proton_misid_cosdat_range = new TH1D("hist_proton_misid_cosdat_range", "", cosdat_2p_bin_size-1, cosdat_2p_bins);

  // Flux systematics
  TH2D *hist_flux_mu_mom_total = new TH2D("hist_flux_mu_mom_total", "",
					  muon_mom_bin_size-1, muon_mom_bins,
					  nu_ene_bin_size-1, nu_ene_bins);
  TH2D *hist_flux_mu_mom_mcs = new TH2D("hist_flux_mu_mom_mcs", "",
					muon_mom_bin_size-1, muon_mom_bins,
					nu_ene_bin_size-1, nu_ene_bins);
  TH2D *hist_flux_mu_mom_range = new TH2D("hist_flux_mu_mom_range", "",
					  muon_mom_bin_size-1, muon_mom_bins,
					  nu_ene_bin_size-1, nu_ene_bins);
  TH2D *hist_flux_p_mom_total = new TH2D("hist_flux_p_mom_total", "",
					 hadron_mom_bin_size-1, hadron_mom_bins,
					 nu_ene_bin_size-1, nu_ene_bins);
  TH2D *hist_flux_p_mom_total_high = new TH2D("hist_flux_p_mom_total_high", "",
					      hadron_mom_bin_size-1, hadron_mom_bins,
					      nu_ene_bin_size-1, nu_ene_bins);
  TH2D *hist_flux_p_mom_total_low = new TH2D("hist_flux_p_mom_total_low", "",
					     hadron_mom_bin_size-1, hadron_mom_bins,
					     nu_ene_bin_size-1, nu_ene_bins);
  TH2D *hist_flux_p_mom_mcs = new TH2D("hist_flux_p_mom_mcs", "",
				       hadron_mom_bin_size-1, hadron_mom_bins,
				       nu_ene_bin_size-1, nu_ene_bins);
  TH2D *hist_flux_p_mom_mcs_high = new TH2D("hist_flux_p_mom_mcs_high", "",
					    hadron_mom_bin_size-1, hadron_mom_bins,
					    nu_ene_bin_size-1, nu_ene_bins);
  TH2D *hist_flux_p_mom_mcs_low = new TH2D("hist_flux_p_mom_mcs_low", "",
					   hadron_mom_bin_size-1, hadron_mom_bins,
					   nu_ene_bin_size-1, nu_ene_bins);
  TH2D *hist_flux_p_mom_range = new TH2D("hist_flux_p_mom_range", "",
					 hadron_mom_bin_size-1, hadron_mom_bins,
					 nu_ene_bin_size-1, nu_ene_bins);
  TH2D *hist_flux_p_mom_range_high = new TH2D("hist_flux_p_mom_range_high", "",
					      hadron_mom_bin_size-1, hadron_mom_bins,
					      nu_ene_bin_size-1, nu_ene_bins);
  TH2D *hist_flux_p_mom_range_low = new TH2D("hist_flux_p_mom_range_low", "",
					     hadron_mom_bin_size-1, hadron_mom_bins,
					     nu_ene_bin_size-1, nu_ene_bins);
  TH2D *hist_flux_mu_deg = new TH2D("hist_flux_mu_deg", "",
				    muon_deg_bin_size-1, muon_deg_bins,
				    nu_ene_bin_size-1, nu_ene_bins);
  TH2D *hist_flux_mu_cos = new TH2D("hist_flux_mu_cos", "",
				    muon_cos_bin_size-1, muon_cos_bins,
				    nu_ene_bin_size-1, nu_ene_bins);
  TH2D *hist_flux_p_deg = new TH2D("hist_flux_p_deg", "",
				   hadron_deg_bin_size-1, hadron_deg_bins,
				   nu_ene_bin_size-1, nu_ene_bins);
  TH2D *hist_flux_p_deg_high = new TH2D("hist_flux_p_deg_high", "",
					hadron_deg_bin_size-1, hadron_deg_bins,
					nu_ene_bin_size-1, nu_ene_bins);
  TH2D *hist_flux_p_deg_low = new TH2D("hist_flux_p_deg_low", "",
				       hadron_deg_bin_size-1, hadron_deg_bins,
				       nu_ene_bin_size-1, nu_ene_bins);
  TH2D *hist_flux_p_cos = new TH2D("hist_flux_p_cos", "",
				   hadron_cos_bin_size-1, hadron_cos_bins,
				   nu_ene_bin_size-1, nu_ene_bins);
  TH2D *hist_flux_p_cos_high = new TH2D("hist_flux_p_cos_high", "",
					hadron_cos_bin_size-1, hadron_cos_bins,
					nu_ene_bin_size-1, nu_ene_bins);
  TH2D *hist_flux_p_cos_low = new TH2D("hist_flux_p_cos_low", "",
				       hadron_cos_bin_size-1, hadron_cos_bins,
				       nu_ene_bin_size-1, nu_ene_bins);
  TH2D *hist_flux_open_ang = new TH2D("hist_flux_open_ang", "",
				      open_deg_bin_size-1, open_deg_bins,
				      nu_ene_bin_size-1, nu_ene_bins);
  TH2D *hist_flux_open_cos = new TH2D("hist_flux_open_cos", "",
				      open_cos_bin_size-1, open_cos_bins,
				      nu_ene_bin_size-1, nu_ene_bins);
  TH2D *hist_flux_mom_ratio = new TH2D("hist_flux_mom_ratio", "",
				       mom_ratio_bin_size-1, mom_ratio_bins,
				       nu_ene_bin_size-1, nu_ene_bins);
  TH2D *hist_flux_mom_vecsum = new TH2D("hist_flux_mom_vecsum", "",
					mom_vecsum_bin_size-1, mom_vecsum_bins,
					nu_ene_bin_size-1, nu_ene_bins);
  TH2D *hist_flux_mom_scasum = new TH2D("hist_flux_mom_scasum", "",
					mom_scasum_bin_size-1, mom_scasum_bins,
					nu_ene_bin_size-1, nu_ene_bins);
  TH2D *hist_flux_dptt = new TH2D("hist_flux_dptt", "",
				  dptt_bin_size-1, dptt_bins,
				  nu_ene_bin_size-1, nu_ene_bins);
  TH2D *hist_flux_dpt = new TH2D("hist_flux_dpt", "",
				 dpt_2p_bin_size-1, dpt_2p_bins,
				 nu_ene_bin_size-1, nu_ene_bins);
  TH2D *hist_flux_pn = new TH2D("hist_flux_pn", "",
				pn_bin_size-1, pn_bins,
				nu_ene_bin_size-1, nu_ene_bins);
  TH2D *hist_flux_dalphat = new TH2D("hist_flux_dalphat", "",
				     dalphat_2p_bin_size-1, dalphat_2p_bins,
				     nu_ene_bin_size-1, nu_ene_bins);
  TH2D *hist_flux_cosdat = new TH2D("hist_flux_cosdat", "",
				    cosdat_2p_bin_size-1, cosdat_2p_bins,
				    nu_ene_bin_size-1, nu_ene_bins);
  TH2D *hist_flux_dptt_range = new TH2D("hist_flux_dptt_range", "",
					dptt_bin_size-1, dptt_bins,
					nu_ene_bin_size-1, nu_ene_bins);
  TH2D *hist_flux_dpt_range = new TH2D("hist_flux_dpt_range", "",
				       dpt_2p_bin_size-1, dpt_2p_bins,
				       nu_ene_bin_size-1, nu_ene_bins);
  TH2D *hist_flux_pn_range = new TH2D("hist_flux_pn_range", "",
				      pn_bin_size-1, pn_bins,
				      nu_ene_bin_size-1, nu_ene_bins);
  TH2D *hist_flux_dalphat_range = new TH2D("hist_flux_dalphat_range", "",
					   dalphat_2p_bin_size-1, dalphat_2p_bins,
					   nu_ene_bin_size-1, nu_ene_bins);
  TH2D *hist_flux_cosdat_range = new TH2D("hist_flux_cosdat_range", "",
					  cosdat_2p_bin_size-1, cosdat_2p_bins,
					  nu_ene_bin_size-1, nu_ene_bins);

  for ( auto ev : ev_vec ) {

    reader.ReadSpill(ev.groupid);
    auto &spill_summary = reader.GetSpillSummary();
    auto it_event = spill_summary.BeginTrueEvent();
    const auto *event = it_event.Next();

    auto &vertex = event->GetPrimaryVertex();
    int mode_id = GetNinjaModeId(vertex.GetInteractionType());

    if ( ev.vertex_material != B2Material::kWater ) continue;

    if ( ev.chains.size() != 3 ) continue;

    bool muon_correct_flag = true;
    bool proton_correct_flag = true;
    int num_muon = 0;
    int num_proton = 0;
    int num_nonproton = 0;
    for ( auto chain : ev.chains ) {
      int recon_particle_flag = chain.particle_flag % 10000;
      int true_particle_flag = chain.particle_flag / 10000;
      if ( recon_particle_flag == 13 ) {
	num_muon++;
	if ( recon_particle_flag != true_particle_flag ) muon_correct_flag = false;
      }
      else if ( recon_particle_flag == 2212 ) {
	num_proton++;
	if ( recon_particle_flag != true_particle_flag ) proton_correct_flag = false;
      }
      else {
	num_nonproton++;
      }
    }

    // CC0pi2p0other
    if ( num_proton != 2 || num_nonproton != 0 || num_muon != 1 ) continue;

    // kinematics calculation
    bool muon_stop_flag = false;
    bool proton_stop_flag_high = false;
    bool proton_stop_flag_low = false;
    TVector3 muon_tangent;
    TVector3 proton_tangent_high;
    TVector3 proton_tangent_low;
    int proton_direction_high;
    int proton_direction_low;
    double muon_momentum;
    double proton_momentum_high;
    double proton_momentum_low;
    TVector3 muon_momentum_vec;
    TVector3 proton_momentum_vec_high;
    TVector3 proton_momentum_vec_low;
    double muon_deg;
    double proton_deg_high;
    double proton_deg_low;
    double muon_cosine;
    double proton_cosine_high;
    double proton_cosine_low;
    num_proton = 0;

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
      
      double tangent = chain.direction * std::hypot(ax, ay);
      double theta = std::atan(tangent);
      double theta_deg = theta * TMath::RadToDeg();
      if ( chain.direction == -1 ) theta_deg += 180.;
      double cosine = std::cos(theta_deg * TMath::DegToRad());

      if ( chain.particle_flag % 10000 == 13 ) {
	muon_tangent.SetXYZ(ax, ay, 1.);
	muon_deg = theta_deg;
	muon_cosine = cosine;
	if ( chain.stop_flag == 0 ) {
	  muon_momentum = chain.ecc_mcs_mom[0];
	  muon_stop_flag = true;
	}
	else if ( chain.stop_flag == 1 )
	  muon_momentum = chain.bm_range_mom;
	muon_momentum_vec = (muon_momentum / muon_tangent.Mag()) * muon_tangent;
      }
      else if ( chain.particle_flag % 10000 == 2212 ) {
	num_proton++;
	if ( num_proton == 1 ) {
	  proton_tangent_high.SetXYZ(ax, ay, 1);
	  proton_direction_high = chain.direction;
	  proton_deg_high = theta_deg;
	  proton_cosine_high = cosine;
	  if ( chain.stop_flag == 2 ) {
	    proton_momentum_high = chain.ecc_range_mom[1];
	    proton_stop_flag_high = true;
	  }
	  else {
	    proton_momentum_high = chain.ecc_mcs_mom[1];
	  }
	  proton_momentum_vec_high = (proton_momentum_high * chain.direction / proton_tangent_high.Mag()) * proton_tangent_high;
	}
	else if ( num_proton == 2 ) {
	  proton_tangent_low.SetXYZ(ax, ay, 1);
	  proton_direction_low = chain.direction;
	  proton_deg_low = theta_deg;
	  proton_cosine_low = cosine;
	  if ( chain.stop_flag == 2 ) {
	    proton_momentum_low = chain.ecc_range_mom[1];
	    proton_stop_flag_low = true;
	  }
	  else {
	    proton_momentum_low = chain.ecc_mcs_mom[1];
	  }
	  proton_momentum_vec_low = (proton_momentum_low * chain.direction / proton_tangent_low.Mag()) * proton_tangent_low;
	}
      }
    }

    if ( proton_momentum_low > proton_momentum_high ) {
      std::swap(proton_direction_low, proton_direction_high);
      std::swap(proton_tangent_low, proton_tangent_high);
      std::swap(proton_deg_low, proton_deg_high);
      std::swap(proton_cosine_low, proton_cosine_high);
      std::swap(proton_stop_flag_low, proton_stop_flag_high);
      std::swap(proton_momentum_low, proton_momentum_high);
      std::swap(proton_momentum_vec_low, proton_momentum_vec_high);
    }

    // Calculate
    double open_cos = (proton_momentum_vec_high * proton_momentum_vec_low)
      / proton_momentum_vec_high.Mag() / proton_momentum_vec_low.Mag();
    double open_ang = std::acos(open_cos) * TMath::RadToDeg();

    TVector3 ztt(-1. * muon_momentum_vec.Y(), muon_momentum_vec.X(), 0.);
    ztt = (1. / std::hypot(muon_momentum_vec.X(), muon_momentum_vec.Y())) * ztt;
    // ztt = pnu x pmu / |pnu x pmu|, where pnu = (0,0,Enu), pmu = muon_momentum_vec.
    double dptt = ztt * proton_momentum_vec_high
      + ztt * proton_momentum_vec_low;
      
    TVector2 mumom_vec_2d(muon_momentum_vec.X(), muon_momentum_vec.Y());
    TVector2 promom_vec_2d_high(proton_momentum_vec_high.X(), proton_momentum_vec_high.Y());
    TVector2 promom_vec_2d_low(proton_momentum_vec_low.X(), proton_momentum_vec_low.Y());
    TVector2 dpt_vec = mumom_vec_2d + promom_vec_2d_high + promom_vec_2d_low;

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

    double dalphat = std::acos((-1. * mumom_vec_2d * dpt_vec) / mumom_vec_2d.Mod() / dpt_vec.Mod());
    double cosdat = (-1. * mumom_vec_2d * dpt_vec) / mumom_vec_2d.Mod() / dpt_vec.Mod();

    if ( muon_correct_flag ) {

      if ( proton_correct_flag ) {

	hist_muon_mom->Fill(muon_momentum, ev.weight);
	hist_proton_mom->Fill(proton_momentum_high, ev.weight);
	hist_proton_mom->Fill(proton_momentum_low, ev.weight);
	hist_proton_mom_high->Fill(proton_momentum_high, ev.weight);
	hist_proton_mom_low->Fill(proton_momentum_low, ev.weight);
	hist_muon_ang->Fill(muon_deg, ev.weight);
	hist_muon_cos->Fill(muon_cosine, ev.weight);
	hist_proton_ang->Fill(proton_deg_high, ev.weight);
	hist_proton_ang->Fill(proton_deg_low, ev.weight);
	hist_proton_ang_high->Fill(proton_deg_high, ev.weight);
	hist_proton_ang_low->Fill(proton_deg_low, ev.weight);
	hist_proton_cos->Fill(proton_cosine_high, ev.weight);
	hist_proton_cos->Fill(proton_cosine_low, ev.weight);
	hist_proton_cos_high->Fill(proton_cosine_high, ev.weight);
	hist_proton_cos_low->Fill(proton_cosine_low, ev.weight);

	hist_mode_muon_mom[mode_id]->Fill(muon_momentum, ev.weight);
	hist_mode_proton_mom[mode_id]->Fill(proton_momentum_high, ev.weight);
	hist_mode_proton_mom[mode_id]->Fill(proton_momentum_low, ev.weight);
	hist_mode_proton_mom_high[mode_id]->Fill(proton_momentum_high, ev.weight);
	hist_mode_proton_mom_low[mode_id]->Fill(proton_momentum_low, ev.weight);
	hist_mode_muon_ang[mode_id]->Fill(muon_deg, ev.weight);
	hist_mode_muon_cos[mode_id]->Fill(muon_cosine, ev.weight);
	hist_mode_proton_ang[mode_id]->Fill(proton_deg_high, ev.weight);
	hist_mode_proton_ang[mode_id]->Fill(proton_deg_low, ev.weight);
	hist_mode_proton_ang_high[mode_id]->Fill(proton_deg_high, ev.weight);
	hist_mode_proton_ang_low[mode_id]->Fill(proton_deg_low, ev.weight);
	hist_mode_proton_cos[mode_id]->Fill(proton_cosine_high, ev.weight);
	hist_mode_proton_cos[mode_id]->Fill(proton_cosine_low, ev.weight);
	hist_mode_proton_cos_high[mode_id]->Fill(proton_cosine_high, ev.weight);
	hist_mode_proton_cos_low[mode_id]->Fill(proton_cosine_low, ev.weight);
	
	hist_muon_mom_ang->Fill(muon_momentum, muon_deg, ev.weight);
	hist_proton_mom_ang->Fill(proton_momentum_high, proton_deg_high, ev.weight);
	hist_proton_mom_ang->Fill(proton_momentum_low, proton_deg_low, ev.weight);
	hist_proton_mom_ang_high->Fill(proton_momentum_high, proton_deg_high, ev.weight);
	hist_proton_mom_ang_low->Fill(proton_momentum_low, proton_deg_low, ev.weight);
	hist_muon_mom_cos->Fill(muon_momentum, muon_cosine, ev.weight);
	hist_proton_mom_cos->Fill(proton_momentum_high, proton_cosine_high, ev.weight);
	hist_proton_mom_cos->Fill(proton_momentum_low, proton_cosine_low, ev.weight);
	hist_proton_mom_cos_high->Fill(proton_momentum_high, proton_cosine_high, ev.weight);
	hist_proton_mom_cos_low->Fill(proton_momentum_low, proton_cosine_low, ev.weight);	 

	hist_open_ang->Fill(open_ang, ev.weight);
	hist_open_cos->Fill(open_cos, ev.weight);
	hist_mom_ratio->Fill(proton_momentum_low / proton_momentum_high, ev.weight);
	hist_mom_vecsum->Fill((proton_momentum_vec_low + proton_momentum_vec_high).Mag(), ev.weight);
	hist_mom_scasum->Fill(proton_momentum_vec_low.Mag() + proton_momentum_vec_high.Mag(), ev.weight);

	hist_mode_open_ang[mode_id]->Fill(open_ang, ev.weight);
	hist_mode_open_cos[mode_id]->Fill(open_cos, ev.weight);
	hist_mode_mom_ratio[mode_id]->Fill(proton_momentum_low / proton_momentum_high, ev.weight);
	hist_mode_mom_vecsum[mode_id]->Fill((proton_momentum_vec_low + proton_momentum_vec_high).Mag(), ev.weight);
	hist_mode_mom_scasum[mode_id]->Fill(proton_momentum_vec_low.Mag() + proton_momentum_vec_high.Mag(), ev.weight);
	
	hist_mom_high_low->Fill(proton_momentum_high, proton_momentum_low, ev.weight);
	hist_ang_high_low->Fill(proton_deg_high, proton_deg_low, ev.weight);
	hist_cos_high_low->Fill(proton_cosine_high, proton_cosine_low, ev.weight);

	hist_dptt->Fill(dptt, ev.weight);
	hist_dpt->Fill(dpt_vec.Mod(), ev.weight);
	hist_pn->Fill(pn, ev.weight);
	hist_dalphat->Fill(dalphat, ev.weight);
	hist_cosdat->Fill(cosdat, ev.weight);

	hist_mode_dptt[mode_id]->Fill(dptt, ev.weight);
	hist_mode_dpt[mode_id]->Fill(dpt_vec.Mod(), ev.weight);
	hist_mode_pn[mode_id]->Fill(pn, ev.weight);
	hist_mode_dalphat[mode_id]->Fill(dalphat, ev.weight);
	hist_mode_cosdat[mode_id]->Fill(cosdat, ev.weight);

	hist_flux_mu_mom_total->Fill(muon_momentum, ev.nu_energy / 1000., ev.weight);
	hist_flux_p_mom_total->Fill(proton_momentum_high, ev.nu_energy / 1000., ev.weight);
	hist_flux_p_mom_total->Fill(proton_momentum_low, ev.nu_energy / 1000., ev.weight);
	hist_flux_p_mom_total_high->Fill(proton_momentum_high, ev.nu_energy / 1000., ev.weight);
	hist_flux_p_mom_total_low->Fill(proton_momentum_low, ev.nu_energy / 1000., ev.weight);
	hist_flux_mu_deg->Fill(muon_deg, ev.nu_energy / 1000., ev.weight);
	hist_flux_mu_cos->Fill(muon_cosine, ev.nu_energy / 1000., ev.weight);
	hist_flux_p_deg->Fill(proton_deg_high, ev.nu_energy / 1000., ev.weight);
	hist_flux_p_deg->Fill(proton_deg_low, ev.nu_energy / 1000., ev.weight);
	hist_flux_p_deg_high->Fill(proton_deg_high, ev.nu_energy / 1000., ev.weight);
	hist_flux_p_deg_low->Fill(proton_deg_low, ev.nu_energy / 1000., ev.weight);
	hist_flux_p_cos->Fill(proton_cosine_high, ev.nu_energy / 1000., ev.weight);
	hist_flux_p_cos->Fill(proton_cosine_low, ev.nu_energy / 1000., ev.weight);
	hist_flux_p_cos_high->Fill(proton_cosine_high, ev.nu_energy / 1000., ev.weight);
	hist_flux_p_cos_low->Fill(proton_cosine_low, ev.nu_energy / 1000., ev.weight);
	hist_flux_dptt->Fill(dptt, ev.nu_energy / 1000., ev.weight);
	hist_flux_dpt->Fill(dpt_vec.Mod(), ev.nu_energy / 1000., ev.weight);
	hist_flux_pn->Fill(pn, ev.nu_energy / 1000., ev.weight);
	hist_flux_dalphat->Fill(dalphat, ev.nu_energy / 1000., ev.weight);
	hist_flux_cosdat->Fill(cosdat, ev.nu_energy / 1000., ev.weight);
	
	
	if ( muon_stop_flag ) {
	  hist_muon_mom_range->Fill(muon_momentum, ev.weight);
	  hist_dptt_range->Fill(dptt, ev.weight);
	  hist_dpt_range->Fill(dpt_vec.Mod(), ev.weight);
	  hist_pn_range->Fill(pn, ev.weight);
	  hist_dalphat_range->Fill(dalphat, ev.weight);
	  hist_cosdat_range->Fill(cosdat, ev.weight);

	  hist_mode_muon_mom_range[mode_id]->Fill(muon_momentum, ev.weight);
	  hist_mode_dptt_range[mode_id]->Fill(dptt, ev.weight);
	  hist_mode_dpt_range[mode_id]->Fill(dpt_vec.Mod(), ev.weight);
	  hist_mode_pn_range[mode_id]->Fill(pn, ev.weight);
	  hist_mode_dalphat_range[mode_id]->Fill(dalphat, ev.weight);
	  hist_mode_cosdat_range[mode_id]->Fill(cosdat, ev.weight);

	  hist_muon_mom_ang_range->Fill(muon_momentum, muon_deg, ev.weight);
	  hist_muon_mom_cos_range->Fill(muon_momentum, muon_cosine, ev.weight);
	  
	  hist_flux_mu_mom_range->Fill(muon_momentum, ev.nu_energy / 1000., ev.weight);
	  hist_flux_dptt_range->Fill(dptt, ev.nu_energy / 1000., ev.weight);
	  hist_flux_dpt_range->Fill(dpt_vec.Mod(), ev.nu_energy / 1000., ev.weight);
	  hist_flux_pn_range->Fill(pn, ev.nu_energy / 1000., ev.weight);
	  hist_flux_dalphat_range->Fill(dalphat, ev.nu_energy / 1000., ev.weight);
	  hist_flux_cosdat_range->Fill(cosdat, ev.nu_energy / 1000., ev.weight);
	}
	else {
	  hist_muon_mom_mcs->Fill(muon_momentum, ev.weight);

	  hist_mode_muon_mom_mcs[mode_id]->Fill(muon_momentum, ev.weight);

	  hist_muon_mom_ang_mcs->Fill(muon_momentum, muon_deg, ev.weight);
	  hist_muon_mom_cos_mcs->Fill(muon_momentum, muon_cosine, ev.weight);
	  
	  hist_flux_mu_mom_mcs->Fill(muon_momentum, ev.nu_energy / 1000., ev.weight);
	}
	
	if ( proton_stop_flag_high ) {
	  hist_proton_mom_range->Fill(proton_momentum_high, ev.weight);
	  hist_proton_mom_range_high->Fill(proton_momentum_high, ev.weight);

	  hist_mode_proton_mom_range[mode_id]->Fill(proton_momentum_high, ev.weight);
	  hist_mode_proton_mom_range_high[mode_id]->Fill(proton_momentum_high, ev.weight);
	  
	  hist_proton_mom_ang_range->Fill(proton_momentum_high, proton_deg_high, ev.weight);
	  hist_proton_mom_cos_range->Fill(proton_momentum_high, proton_cosine_high, ev.weight);
	  hist_proton_mom_ang_range_high->Fill(proton_momentum_high, proton_deg_high, ev.weight);
	  hist_proton_mom_cos_range_high->Fill(proton_momentum_high, proton_cosine_high, ev.weight);

	  hist_flux_p_mom_range->Fill(proton_momentum_high, ev.nu_energy / 1000., ev.weight);
	  hist_flux_p_mom_range_high->Fill(proton_momentum_high, ev.nu_energy / 1000., ev.weight);
	}
	else {
	  hist_proton_mom_mcs->Fill(proton_momentum_high, ev.weight);
	  hist_proton_mom_mcs_high->Fill(proton_momentum_high, ev.weight);

	  hist_mode_proton_mom_mcs[mode_id]->Fill(proton_momentum_high, ev.weight);
	  hist_mode_proton_mom_mcs_high[mode_id]->Fill(proton_momentum_high, ev.weight);

	  hist_proton_mom_ang_mcs->Fill(proton_momentum_high, proton_deg_high, ev.weight);
	  hist_proton_mom_cos_mcs->Fill(proton_momentum_high, proton_cosine_high, ev.weight);
	  hist_proton_mom_ang_mcs_high->Fill(proton_momentum_high, proton_deg_high, ev.weight);
	  hist_proton_mom_cos_mcs_high->Fill(proton_momentum_high, proton_cosine_high, ev.weight);

	  hist_flux_p_mom_mcs->Fill(proton_momentum_high, ev.nu_energy / 1000., ev.weight);
	  hist_flux_p_mom_mcs_high->Fill(proton_momentum_high, ev.nu_energy / 1000., ev.weight);
	}
	
	if ( proton_stop_flag_low ) {
	  hist_proton_mom_range->Fill(proton_momentum_low, ev.weight);
	  hist_proton_mom_range_low->Fill(proton_momentum_low, ev.weight);

	  hist_mode_proton_mom_range[mode_id]->Fill(proton_momentum_low, ev.weight);
	  hist_mode_proton_mom_range_low[mode_id]->Fill(proton_momentum_low, ev.weight);

	  hist_proton_mom_ang_range->Fill(proton_momentum_low, proton_deg_low, ev.weight);
	  hist_proton_mom_cos_range->Fill(proton_momentum_low, proton_cosine_low, ev.weight);
	  hist_proton_mom_ang_range_low->Fill(proton_momentum_low, proton_deg_low, ev.weight);
	  hist_proton_mom_cos_range_low->Fill(proton_momentum_low, proton_cosine_low, ev.weight);

	  hist_flux_p_mom_range->Fill(proton_momentum_low, ev.nu_energy / 1000., ev.weight);
	  hist_flux_p_mom_range_low->Fill(proton_momentum_low, ev.nu_energy / 1000., ev.weight);
	}
	else {
	  hist_proton_mom_mcs->Fill(proton_momentum_low, ev.weight);
	  hist_proton_mom_mcs_low->Fill(proton_momentum_low, ev.weight);

	  hist_mode_proton_mom_mcs[mode_id]->Fill(proton_momentum_low, ev.weight);
	  hist_mode_proton_mom_mcs_low[mode_id]->Fill(proton_momentum_low, ev.weight);

	  hist_proton_mom_ang_mcs->Fill(proton_momentum_low, proton_deg_low, ev.weight);
	  hist_proton_mom_cos_mcs->Fill(proton_momentum_low, proton_cosine_low, ev.weight);
	  hist_proton_mom_ang_mcs_low->Fill(proton_momentum_low, proton_deg_low, ev.weight);
	  hist_proton_mom_cos_mcs_low->Fill(proton_momentum_low, proton_cosine_low, ev.weight);

	  hist_flux_p_mom_mcs->Fill(proton_momentum_low, ev.nu_energy / 1000., ev.weight);
	  hist_flux_p_mom_mcs_low->Fill(proton_momentum_low, ev.nu_energy / 1000., ev.weight);
	}	
	
      }
      else {
	hist_proton_misid_muon_mom->Fill(muon_momentum, ev.weight);
	hist_proton_misid_proton_mom->Fill(proton_momentum_high, ev.weight);
	hist_proton_misid_proton_mom->Fill(proton_momentum_low, ev.weight);
	hist_proton_misid_proton_mom_high->Fill(proton_momentum_high, ev.weight);
	hist_proton_misid_proton_mom_low->Fill(proton_momentum_low, ev.weight);
	hist_proton_misid_muon_ang->Fill(muon_deg, ev.weight);
	hist_proton_misid_muon_cos->Fill(muon_cosine, ev.weight);
	hist_proton_misid_proton_ang->Fill(proton_deg_high, ev.weight);
	hist_proton_misid_proton_ang->Fill(proton_deg_low, ev.weight);
	hist_proton_misid_proton_ang_high->Fill(proton_deg_high, ev.weight);
	hist_proton_misid_proton_ang_low->Fill(proton_deg_low, ev.weight);
	hist_proton_misid_proton_cos->Fill(proton_cosine_high, ev.weight);
	hist_proton_misid_proton_cos->Fill(proton_cosine_low, ev.weight);
	hist_proton_misid_proton_cos_high->Fill(proton_cosine_high, ev.weight);
	hist_proton_misid_proton_cos_low->Fill(proton_cosine_low, ev.weight);
	
	hist_proton_misid_open_ang->Fill(open_ang, ev.weight);
	hist_proton_misid_open_cos->Fill(open_cos, ev.weight);
	hist_proton_misid_mom_ratio->Fill(proton_momentum_low / proton_momentum_high, ev.weight);
	hist_proton_misid_mom_vecsum->Fill((proton_momentum_vec_low + proton_momentum_vec_high).Mag(), ev.weight);
	hist_proton_misid_mom_scasum->Fill(proton_momentum_vec_low.Mag() + proton_momentum_vec_high.Mag(), ev.weight);
	
	hist_proton_misid_dptt->Fill(dptt, ev.weight);
	hist_proton_misid_dpt->Fill(dpt_vec.Mod(), ev.weight);
	hist_proton_misid_pn->Fill(pn, ev.weight);
	hist_proton_misid_dalphat->Fill(dalphat, ev.weight);
	hist_proton_misid_cosdat->Fill(cosdat, ev.weight);
	
	if ( muon_stop_flag ) {
	  hist_proton_misid_muon_mom_range->Fill(muon_momentum, ev.weight);
	  hist_proton_misid_dptt_range->Fill(dptt, ev.weight);
	  hist_proton_misid_dpt_range->Fill(dpt_vec.Mod(), ev.weight);
	  hist_proton_misid_pn_range->Fill(pn, ev.weight);
	  hist_proton_misid_dalphat_range->Fill(dalphat, ev.weight);
	  hist_proton_misid_cosdat_range->Fill(cosdat, ev.weight);
	}
	else {
	  hist_proton_misid_muon_mom_mcs->Fill(muon_momentum, ev.weight);
	}
	
	if ( proton_stop_flag_high ) {
	  hist_proton_misid_proton_mom_range->Fill(proton_momentum_high, ev.weight);
	  hist_proton_misid_proton_mom_range_high->Fill(proton_momentum_high, ev.weight);
	}
	else {
	  hist_proton_misid_proton_mom_mcs->Fill(proton_momentum_high, ev.weight);
	  hist_proton_misid_proton_mom_mcs_high->Fill(proton_momentum_high, ev.weight);
	}
	
	if ( proton_stop_flag_low ) {
	  hist_proton_misid_proton_mom_range->Fill(proton_momentum_low, ev.weight);
	  hist_proton_misid_proton_mom_range_low->Fill(proton_momentum_low, ev.weight);
	}
	else {
	  hist_proton_misid_proton_mom_mcs->Fill(proton_momentum_low, ev.weight);
	  hist_proton_misid_proton_mom_mcs_low->Fill(proton_momentum_low, ev.weight);
	}	
      }
    }
    else {
      hist_muon_misid_muon_mom->Fill(muon_momentum, ev.weight);
      hist_muon_misid_proton_mom->Fill(proton_momentum_high, ev.weight);
      hist_muon_misid_proton_mom->Fill(proton_momentum_low, ev.weight);
      hist_muon_misid_proton_mom_high->Fill(proton_momentum_high, ev.weight);
      hist_muon_misid_proton_mom_low->Fill(proton_momentum_low, ev.weight);
      hist_muon_misid_muon_ang->Fill(muon_deg, ev.weight);
      hist_muon_misid_muon_cos->Fill(muon_cosine, ev.weight);
      hist_muon_misid_proton_ang->Fill(proton_deg_high, ev.weight);
      hist_muon_misid_proton_ang->Fill(proton_deg_low, ev.weight);
      hist_muon_misid_proton_ang_high->Fill(proton_deg_high, ev.weight);
      hist_muon_misid_proton_ang_low->Fill(proton_deg_low, ev.weight);
      hist_muon_misid_proton_cos->Fill(proton_cosine_high, ev.weight);
      hist_muon_misid_proton_cos->Fill(proton_cosine_low, ev.weight);
      hist_muon_misid_proton_cos_high->Fill(proton_cosine_high, ev.weight);
      hist_muon_misid_proton_cos_low->Fill(proton_cosine_low, ev.weight);

      hist_muon_misid_open_ang->Fill(open_ang, ev.weight);
      hist_muon_misid_open_cos->Fill(open_cos, ev.weight);
      hist_muon_misid_mom_ratio->Fill(proton_momentum_low / proton_momentum_high, ev.weight);
      hist_muon_misid_mom_vecsum->Fill((proton_momentum_vec_low + proton_momentum_vec_high).Mag(), ev.weight);
      hist_muon_misid_mom_scasum->Fill(proton_momentum_vec_low.Mag() + proton_momentum_vec_high.Mag(), ev.weight);

      hist_muon_misid_dptt->Fill(dptt, ev.weight);
      hist_muon_misid_dpt->Fill(dpt_vec.Mod(), ev.weight);
      hist_muon_misid_pn->Fill(pn, ev.weight);
      hist_muon_misid_dalphat->Fill(dalphat, ev.weight);
      hist_muon_misid_cosdat->Fill(cosdat, ev.weight);
      
      if ( muon_stop_flag ) {
	hist_muon_misid_muon_mom_range->Fill(muon_momentum, ev.weight);
	hist_muon_misid_dptt_range->Fill(dptt, ev.weight);
	hist_muon_misid_dpt_range->Fill(dpt_vec.Mod(), ev.weight);
	hist_muon_misid_pn_range->Fill(pn, ev.weight);
	hist_muon_misid_dalphat_range->Fill(dalphat, ev.weight);
	hist_muon_misid_cosdat_range->Fill(cosdat, ev.weight);
      }
      else {
	hist_muon_misid_muon_mom_mcs->Fill(muon_momentum, ev.weight);
      }

      if ( proton_stop_flag_high ) {
	hist_muon_misid_proton_mom_range->Fill(proton_momentum_high, ev.weight);
	hist_muon_misid_proton_mom_range_high->Fill(proton_momentum_high, ev.weight);
      }
      else {
	hist_muon_misid_proton_mom_mcs->Fill(proton_momentum_high, ev.weight);
	hist_muon_misid_proton_mom_mcs_high->Fill(proton_momentum_high, ev.weight);
      }

      if ( proton_stop_flag_low ) {
	hist_muon_misid_proton_mom_range->Fill(proton_momentum_low, ev.weight);
	hist_muon_misid_proton_mom_range_low->Fill(proton_momentum_low, ev.weight);
      }
      else {
	hist_muon_misid_proton_mom_mcs->Fill(proton_momentum_low, ev.weight);
	hist_muon_misid_proton_mom_mcs_low->Fill(proton_momentum_low, ev.weight);
      }

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
  hist_muon_mom_range->Write();
  hist_proton_mom_range->Write();
  hist_proton_mom_range_high->Write();
  hist_proton_mom_range_low->Write();
  hist_muon_mom_mcs->Write();
  hist_proton_mom_mcs->Write();
  hist_proton_mom_mcs_high->Write();
  hist_proton_mom_mcs_low->Write();
  hist_muon_mom_ang_mcs->Write();
  hist_muon_mom_cos_mcs->Write();
  hist_proton_mom_ang_mcs->Write();
  hist_proton_mom_ang_mcs_high->Write();
  hist_proton_mom_ang_mcs_low->Write();
  hist_proton_mom_cos_mcs->Write();
  hist_proton_mom_cos_mcs_high->Write();
  hist_proton_mom_cos_mcs_low->Write();
  hist_muon_mom_ang_range->Write();
  hist_muon_mom_cos_range->Write();
  hist_proton_mom_ang_range->Write();
  hist_proton_mom_ang_range_high->Write();
  hist_proton_mom_ang_range_low->Write();
  hist_proton_mom_cos_range->Write();
  hist_proton_mom_cos_range_high->Write();
  hist_proton_mom_cos_range_low->Write();

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
  hist_dptt_range->Write();
  hist_dpt_range->Write();
  hist_pn_range->Write();
  hist_dalphat_range->Write();
  hist_cosdat_range->Write();

  for ( int i = 0; i < num_ninja_mode; i++ ) {
    hist_mode_muon_mom[i]->Write();
    hist_mode_proton_mom[i]->Write();
    hist_mode_proton_mom_high[i]->Write();
    hist_mode_proton_mom_low[i]->Write();
    hist_mode_muon_mom_mcs[i]->Write();
    hist_mode_proton_mom_mcs[i]->Write();
    hist_mode_proton_mom_mcs_high[i]->Write();
    hist_mode_proton_mom_mcs_low[i]->Write();
    hist_mode_muon_mom_range[i]->Write();
    hist_mode_proton_mom_range[i]->Write();
    hist_mode_proton_mom_range_high[i]->Write();
    hist_mode_proton_mom_range_low[i]->Write();
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
    hist_mode_dptt_range[i]->Write();
    hist_mode_dpt_range[i]->Write();
    hist_mode_pn_range[i]->Write();
    hist_mode_dalphat_range[i]->Write();
    hist_mode_cosdat_range[i]->Write();
  }

  hist_muon_misid_muon_mom->Write();
  hist_muon_misid_proton_mom->Write();
  hist_muon_misid_proton_mom_high->Write();
  hist_muon_misid_proton_mom_low->Write();
  hist_muon_misid_muon_ang->Write();
  hist_muon_misid_muon_cos->Write();
  hist_muon_misid_proton_ang->Write();
  hist_muon_misid_proton_ang_high->Write();
  hist_muon_misid_proton_ang_low->Write();
  hist_muon_misid_proton_cos->Write();
  hist_muon_misid_proton_cos_high->Write();
  hist_muon_misid_proton_cos_low->Write();
  hist_muon_misid_muon_mom_mcs->Write();
  hist_muon_misid_muon_mom_range->Write();
  hist_muon_misid_proton_mom_mcs->Write();
  hist_muon_misid_proton_mom_mcs_high->Write();
  hist_muon_misid_proton_mom_mcs_low->Write();
  hist_muon_misid_proton_mom_range->Write();
  hist_muon_misid_proton_mom_range_high->Write();
  hist_muon_misid_proton_mom_range_low->Write();
  hist_muon_misid_open_ang->Write();
  hist_muon_misid_open_cos->Write();
  hist_muon_misid_mom_ratio->Write();
  hist_muon_misid_mom_vecsum->Write();
  hist_muon_misid_mom_scasum->Write();
  hist_muon_misid_dptt->Write();
  hist_muon_misid_dpt->Write();
  hist_muon_misid_pn->Write();
  hist_muon_misid_dalphat->Write();
  hist_muon_misid_cosdat->Write();
  hist_muon_misid_dptt_range->Write();
  hist_muon_misid_dpt_range->Write();
  hist_muon_misid_pn_range->Write();
  hist_muon_misid_dalphat_range->Write();
  hist_muon_misid_cosdat_range->Write();

  hist_proton_misid_muon_mom->Write();
  hist_proton_misid_proton_mom->Write();
  hist_proton_misid_proton_mom_high->Write();
  hist_proton_misid_proton_mom_low->Write();
  hist_proton_misid_muon_ang->Write();
  hist_proton_misid_muon_cos->Write();
  hist_proton_misid_proton_ang->Write();
  hist_proton_misid_proton_ang_high->Write();
  hist_proton_misid_proton_ang_low->Write();
  hist_proton_misid_proton_cos->Write();
  hist_proton_misid_proton_cos_high->Write();
  hist_proton_misid_proton_cos_low->Write();
  hist_proton_misid_muon_mom_mcs->Write();
  hist_proton_misid_muon_mom_range->Write();
  hist_proton_misid_proton_mom_mcs->Write();
  hist_proton_misid_proton_mom_mcs_high->Write();
  hist_proton_misid_proton_mom_mcs_low->Write();
  hist_proton_misid_proton_mom_range->Write();
  hist_proton_misid_proton_mom_range_high->Write();
  hist_proton_misid_proton_mom_range_low->Write();
  hist_proton_misid_open_ang->Write();
  hist_proton_misid_open_cos->Write();
  hist_proton_misid_mom_ratio->Write();
  hist_proton_misid_mom_vecsum->Write();
  hist_proton_misid_mom_scasum->Write();
  hist_proton_misid_dptt->Write();
  hist_proton_misid_dpt->Write();
  hist_proton_misid_pn->Write();
  hist_proton_misid_dalphat->Write();
  hist_proton_misid_cosdat->Write();
  hist_proton_misid_dptt_range->Write();
  hist_proton_misid_dpt_range->Write();
  hist_proton_misid_pn_range->Write();
  hist_proton_misid_dalphat_range->Write();
  hist_proton_misid_cosdat_range->Write();

  hist_flux_mu_mom_total->Write();
  hist_flux_mu_mom_mcs->Write();
  hist_flux_mu_mom_range->Write();
  hist_flux_p_mom_total->Write();
  hist_flux_p_mom_total_high->Write();
  hist_flux_p_mom_total_low->Write();
  hist_flux_p_mom_mcs->Write();
  hist_flux_p_mom_mcs_high->Write();
  hist_flux_p_mom_mcs_low->Write();
  hist_flux_p_mom_range->Write();
  hist_flux_p_mom_range_high->Write();
  hist_flux_p_mom_range_low->Write();
  hist_flux_mu_deg->Write();
  hist_flux_mu_cos->Write();
  hist_flux_p_deg->Write();
  hist_flux_p_deg_high->Write();
  hist_flux_p_deg_low->Write();
  hist_flux_p_cos->Write();
  hist_flux_p_cos_high->Write();
  hist_flux_p_cos_low->Write();
  hist_flux_open_ang->Write();
  hist_flux_open_cos->Write();
  hist_flux_mom_ratio->Write();
  hist_flux_mom_vecsum->Write();
  hist_flux_mom_scasum->Write();
  hist_flux_dptt->Write();
  hist_flux_dpt->Write();
  hist_flux_pn->Write();
  hist_flux_dalphat->Write();
  hist_flux_cosdat->Write();
  hist_flux_dptt_range->Write();
  hist_flux_dpt_range->Write();
  hist_flux_pn_range->Write();
  hist_flux_dalphat_range->Write();
  hist_flux_cosdat_range->Write();

  outputfile->Close();

}
