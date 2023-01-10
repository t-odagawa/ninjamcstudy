#include <iostream>
#include <vector>
#include <cmath>
#include <string>

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/filesystem.hpp>

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>

#include <B2Reader.hh>

#include <McsClass.hpp>

#include "DrawConst.hpp"
#include "HistogramStyle.hpp"

namespace logging = boost::log;
namespace fs = boost::filesystem;

int main (int argc, char* argv[]) {

  logging::core::get()->set_filter
    (
     logging::trivial::severity >= logging::trivial::info
     );

  BOOST_LOG_TRIVIAL(info) << "==========Anti Muon neutrino-water interaction background==========";

  if ( argc != 3 ) {
    BOOST_LOG_TRIVIAL(error) << "Usage : " << argv[0]
			     << " <input momch file> <output file>";
    std::exit(1);
  }

  std::string momchfilename = argv[1];
  auto ev_vec = Momentum_recon::ReadEventInformationBin(momchfilename);

  TFile *ofile = new TFile((TString)argv[2], "recreate");
  
  TH1D *hist_recon_vertex_x = new TH1D("hist_recon_vertex_x", "Reconstructed horizontal position;x [#mum];Entries", 25, 0, 250.e3);
  TH1D *hist_recon_vertex_y = new TH1D("hist_recon_vertex_y", "Reconstructed vertical position;y [#mum];Entries", 25, 0., 250.e3);
  TH1D *hist_recon_vertex_z = new TH1D("hist_recon_vertex_z", "Reconstructed depth;z [#mum];Entries", 25, -250.e3, 0.);
  TH2D *hist_recon_vertex_xy = new TH2D("hist_recon_vertex_xy", "Reconstructed film position;x [#mum];y [#mum]", 250, 0., 250.e3, 250, 0., 250.e3);

  TH1D *hist_vertex_pl = new TH1D("hist_vertex_pl", "", 133, 0.5, 133.5);

  TH1D *hist_multi = new TH1D("hist_multi", "Multiplicity;# of tracks;Entries", total_multi_bin_size - 1, total_multi_bins);
  TH2D *hist_multi_z = new TH2D("hist_multi_z", ";# of tracks;Reconstructed depth [#mum]", 250, -250.e3, 0., 5, 0.5, 5.5);
  TH1D *hist_plate_diff = new TH1D("hist_plate_diff", ";Plate difference;Plate difference", 20,-10,10);
  TH1D *hist_charge = new TH1D("hist_charge", ";Reconstructed charge;Entries",3, -1.5, 1.5);

  TH1D *hist_multi_p = new TH1D("hist_multi_p", "", hadron_multi_bin_size - 1, hadron_multi_bins);
  TH1D *hist_multi_pi = new TH1D("hist_multi_pi", "", hadron_multi_bin_size - 1, hadron_multi_bins);
  TH1D *hist_multi_other = new TH1D("hist_multi_other", "", hadron_multi_bin_size - 1, hadron_multi_bins);

  TH1D *hist_anu_bg_mu_mom = new TH1D("hist_anu_bg_mu_mom", "", muon_mom_bin_size - 1, muon_mom_bins);
  TH1D *hist_anu_bg_mu_mom_mcs = new TH1D("hist_anu_bg_mu_mom_mcs", "", muon_mom_bin_size - 1, muon_mom_bins);
  TH1D *hist_anu_bg_mu_mom_range = new TH1D("hist_anu_bg_mu_mom_range", "", muon_mom_bin_size - 1, muon_mom_bins);

  TH1D *hist_anu_bg_mu_mom_mcs_all = new TH1D("hist_anu_bg_mu_mom_mcs_all", "", muon_mom_bin_size - 1, muon_mom_bins);
  
  TH1D *hist_anu_bg_mu_mom_true = new TH1D("hist_anu_bg_mu_mom_true", "", muon_mom_bin_size - 1, muon_mom_bins);
  TH1D *hist_anu_bg_mu_mom_true_stop = new TH1D("hist_anu_mu_mom_true_stop", "", muon_mom_bin_size - 1, muon_mom_bins);
  TH1D *hist_anu_bg_mu_mom_true_nonstop = new TH1D("hist_anu_mu_mom_true_nonstop", "", muon_mom_bin_size - 1, muon_mom_bins);

  TH1D *hist_anu_bg_mu_cos = new TH1D("hist_anu_bg_mu_cos", "", muon_cos_bin_size - 1, muon_cos_bins);
  TH1D *hist_anu_bg_mu_deg = new TH1D("hist_anu_bg_mu_deg", "", muon_deg_bin_size - 1, muon_deg_bins);

  TH1D *hist_anu_bg_pi_mom = new TH1D("hist_anu_bg_pi_mom", ";Reconstructed momentum [MeV/c];Entries", hadron_mom_bin_size - 1, hadron_mom_bins);
  TH1D *hist_anu_bg_pi_mom_mcs = new TH1D("hist_anu_bg_pi_mom_mcs", ";Reconstructed momentum [MeV/c];Entries", hadron_mom_bin_size - 1, hadron_mom_bins);
  TH1D *hist_anu_bg_pi_mom_range = new TH1D("hist_anu_bg_pi_mom_range", ";Reconstructed momentum [MeV/c];Entries",
					    hadron_mom_bin_size - 1, hadron_mom_bins);
  TH1D *hist_anu_bg_p_mom = new TH1D("hist_anu_bg_p_mom", ";Reconstructed momentum [MeV/c];Entries", hadron_mom_bin_size - 1, hadron_mom_bins);
  TH1D *hist_anu_bg_p_mom_mcs = new TH1D("hist_anu_bg_p_mom_mcs", ";Reconstructed momentum [MeV/c];Entries", hadron_mom_bin_size - 1, hadron_mom_bins);
  TH1D *hist_anu_bg_p_mom_range = new TH1D("hist_anu_bg_p_mom_range", ";Reconstructed momentum [MeV/c];Entries",
					   hadron_mom_bin_size - 1, hadron_mom_bins);
  TH1D *hist_anu_bg_pi_cos = new TH1D("hist_anu_bg_pi_cos", ";Reconstructed angle;Entries", hadron_cos_bin_size - 1, hadron_cos_bins);
  TH1D *hist_anu_bg_pi_deg = new TH1D("hist_anu_bg_pi_deg", ";Reconstructed angle [degree];Entries", hadron_deg_bin_size - 1, hadron_deg_bins);
  TH1D *hist_anu_bg_p_cos = new TH1D("hist_anu_bg_p_cos", ";Reconstructed angle;Entries", hadron_cos_bin_size - 1, hadron_cos_bins);
  TH1D *hist_anu_bg_p_deg = new TH1D("hist_anu_bg_p_deg", ";Reconstructed angle [degree];Entries", hadron_deg_bin_size - 1, hadron_deg_bins);

  TH1D *hist_anu_bg_other_mom_mcs = new TH1D("hist_anu_bg_other_mom_mcs", "", hadron_mom_bin_size - 1, hadron_mom_bins);
  TH1D *hist_anu_bg_other_cos = new TH1D("hist_anu_bg_other_cos", "", hadron_cos_bin_size - 1, hadron_cos_bins);
  TH1D *hist_anu_bg_other_deg = new TH1D("hist_anu_bg_other_deg", "", hadron_deg_bin_size - 1, hadron_deg_bins);

  TH1D *hist_side_recon_vertex_x = new TH1D("hist_side_recon_vertex_x", "Reconstructed horizontal position;x [#mum];Entries", 25, 0, 250.e3);
  TH1D *hist_side_recon_vertex_y = new TH1D("hist_side_recon_vertex_y", "Reconstructed vertical position;y [#mum];Entries", 25, 0., 250.e3);
  TH1D *hist_side_recon_vertex_z = new TH1D("hist_side_recon_vertex_z", "Reconstructed depth;z [#mum];Entries", 25, -250.e3, 0.);
  TH2D *hist_side_recon_vertex_xy = new TH2D("hist_side_recon_vertex_xy", "Reconstructed film position;x [#mum];y [#mum]", 250, 0., 250.e3, 250, 0., 250.e3);

  TH1D *hist_side_multi = new TH1D("hist_side_multi", "Multiplicity;# of tracks;Entries", 5, 0.5, 5.5);
  TH2D *hist_side_multi_z = new TH2D("hist_side_multi_z", ";# of tracks;Reconstructed depth [#mum]", 250, -250.e3, 0., 5, 0.5, 5.5);

  TH1D *hist_pene_recon_vertex_x = new TH1D("hist_pene_recon_vertex_x", "Reconstructed horizontal position;x [#mum];Entries", 25, 0, 250.e3);
  TH1D *hist_pene_recon_vertex_y = new TH1D("hist_pene_recon_vertex_y", "Reconstructed vertical position;y [#mum];Entries", 25, 0., 250.e3);
  TH1D *hist_pene_recon_vertex_z = new TH1D("hist_pene_recon_vertex_z", "Reconstructed depth;z [#mum];Entries", 25, -250.e3, 0.);
  TH2D *hist_pene_recon_vertex_xy = new TH2D("hist_pene_recon_vertex_xy", "Reconstructed film position;x [#mum];y [#mum]", 250, 0., 250.e3, 250, 0., 250.e3);

  TH1D *hist_pene_multi = new TH1D("hist_pene_multi", "Multiplicity;# of tracks;Entries", 5, 0.5, 5.5);
  TH2D *hist_pene_multi_z = new TH2D("hist_pene_multi_z", ";# of tracks;Reconstructed depth [#mum]", 250, -250.e3, 0., 5, 0.5, 5.5);

  // 0pi1p TKI
  TH1D *hist_dpt = new TH1D("hist_dpt", "", dpt_bin_size-1, dpt_bins);
  TH1D *hist_dalphat = new TH1D("hist_dalphat", "", dalphat_bin_size-1, dalphat_bins);
  TH1D *hist_dphit = new TH1D("hist_dphit", "", dphit_bin_size-1, dphit_bins);
  TH1D *hist_dptx = new TH1D("hist_dptx", "", dptx_bin_size-1, dptx_bins);
  TH1D *hist_dpty = new TH1D("hist_dpty", "", dpty_bin_size-1, dpty_bins);

  // 0pi2p protons
  TH1D *hist_open_cos = new TH1D("hist_open_cos", "", open_cos_bin_size-1, open_cos_bins);

  for ( auto ev : ev_vec ) {
    
    if ( ev.chains.empty() ) continue;

    if ( ev.vertex_material == 0 ) {
      hist_recon_vertex_x->Fill(ev.recon_vertex_position[0], ev.weight);
      hist_recon_vertex_y->Fill(ev.recon_vertex_position[1], ev.weight);
      hist_recon_vertex_z->Fill(ev.recon_vertex_position[2], ev.weight);
      hist_recon_vertex_xy->Fill(ev.recon_vertex_position[0], ev.recon_vertex_position[1], ev.weight);
      hist_vertex_pl->Fill(ev.vertex_pl % 1000, ev.weight);
      hist_multi->Fill(ev.chains.size(), ev.weight);
      hist_multi_z->Fill(ev.recon_vertex_position[2], ev.chains.size(), ev.weight);
      hist_plate_diff->Fill(ev.vertex_pl / 1000 - ev.vertex_pl % 1000, ev.weight);

      int num_proton = 0;
      int num_pion = 0;
      for ( auto chain : ev.chains ) {
	int recon_particle_id = chain.particle_flag % 10000;
	if ( recon_particle_id == 2212 ) num_proton++;
	else if ( recon_particle_id == 211 ) num_pion++;
      }

      if ( ev.chains.size() != num_proton + 1 ) continue;

      hist_multi_p->Fill((double)num_proton, ev.weight);
      hist_multi_pi->Fill((double)num_pion, ev.weight);
      hist_multi_other->Fill(ev.chains.size() - num_proton - num_pion - 1, ev.weight);

      // ループが分かれているのは突貫工事だから
      for ( auto chain : ev.chains ) {

	int recon_particle_id = chain.particle_flag % 10000;

	double momentum = -1.;
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
	double tangent = chain.direction * std::hypot(std::tan(thetax), std::tan(thetay));
	double theta = std::atan(tangent);
	double theta_deg = theta * TMath::RadToDeg();
	if ( chain.direction == -1 ) {
	  theta_deg += 180.;
	}
	double cosine = std::cos(theta_deg * TMath::DegToRad());

	if ( recon_particle_id == 13 ) {

	  double mu_mom_true = -1.;
	  for ( auto true_chain : ev.true_chains ) {
	    if ( chain.chainid == true_chain.chainid ) {
	      mu_mom_true = true_chain.bm_range_mom;
	      break;
	    }
	  }

	  if ( chain.stop_flag == 1 ) {
	    momentum = chain.bm_range_mom;
	    hist_anu_bg_mu_mom_range->Fill(momentum, ev.weight);
	    hist_anu_bg_mu_mom_true_stop->Fill(mu_mom_true, ev.weight);
	  }
	  else if ( chain.stop_flag == 0 ) {
	    momentum = chain.ecc_mcs_mom[0];
	    hist_anu_bg_mu_mom_mcs->Fill(momentum, ev.weight);
	    hist_anu_bg_mu_mom_true_nonstop->Fill(mu_mom_true, ev.weight);
	  }
	  hist_anu_bg_mu_mom->Fill(momentum, ev.weight);
	  hist_anu_bg_mu_deg->Fill(theta_deg, ev.weight);
	  hist_anu_bg_mu_cos->Fill(cosine, ev.weight);
	  hist_charge->Fill(chain.charge_sign, ev.weight);

	  hist_anu_bg_mu_mom_mcs_all->Fill(chain.ecc_mcs_mom[0], ev.weight);

	  hist_anu_bg_mu_mom_true->Fill(mu_mom_true, ev.weight);
	}
	else if ( recon_particle_id == 2212 ) {
	  if ( chain.stop_flag == 0 ) {
	    momentum = chain.ecc_mcs_mom[1];
	    hist_anu_bg_p_mom_mcs->Fill(momentum, ev.weight);
	  }
	  else if ( chain.stop_flag == 2 ) {
	    momentum = chain.ecc_range_mom[1];
	    hist_anu_bg_p_mom_range->Fill(momentum, ev.weight);
	  }
	  hist_anu_bg_p_mom->Fill(momentum, ev.weight);
	  hist_anu_bg_p_cos->Fill(cosine, ev.weight);
	  hist_anu_bg_p_deg->Fill(theta_deg, ev.weight);
	}
	else if ( recon_particle_id == 211 ) {
	  if ( chain.stop_flag == 0 ) {
	    momentum = chain.bm_curvature_mom;
	    hist_anu_bg_pi_mom->Fill(momentum, ev.weight);
	    hist_anu_bg_pi_mom_mcs->Fill(momentum, ev.weight);
	  }
	  else if ( chain.stop_flag == 2 ) {
	    momentum = chain.bm_curvature_mom;
	    hist_anu_bg_pi_mom->Fill(momentum, ev.weight);
	    hist_anu_bg_pi_mom_mcs->Fill(momentum, ev.weight);
	    // momentum = chain.ecc_range_mom[0];
	    // hist_anu_bg_pi_mom->Fill(momentum, ev.weight);
	    // hist_anu_bg_pi_mom_range->Fill(momentum, ev.weight);
	  }
	  hist_anu_bg_pi_cos->Fill(cosine, ev.weight);
	  hist_anu_bg_pi_deg->Fill(theta_deg, ev.weight);
	}
	else {
	  momentum = chain.ecc_mcs_mom[0];
	  hist_anu_bg_other_mom_mcs->Fill(momentum, ev.weight);
	  hist_anu_bg_other_cos->Fill(cosine, ev.weight);
	  hist_anu_bg_other_deg->Fill(theta_deg, ev.weight);
	}
      }
      
      // 0pi1p
      if ( num_proton == 1 && ev.chains.size() == 2 ) {
	double muon_momentum;
	double proton_momentum;
	TVector3 muon_tangent;
	TVector3 proton_tangent;
	int proton_direction;
	TVector3 muon_momentum_vec;
	TVector3 proton_momentum_vec;
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

	  if ( chain.particle_flag % 10000 == 13 ) {
	    muon_tangent.SetXYZ(ax, ay, 1.);
	    if ( chain.stop_flag == 1 ) {
	      muon_momentum = chain.bm_range_mom;
	    }
	    else if ( chain.stop_flag == 0 ) {
	      muon_momentum = chain.ecc_mcs_mom[0];
	    }
	    muon_momentum_vec = (muon_momentum / muon_tangent.Mag()) * muon_tangent;
	  }
	  else if ( chain.particle_flag % 10000 == 2212 ) {
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
	    }
	    else {
	      proton_momentum = chain.ecc_mcs_mom[1];
	    }
	    proton_momentum_vec = (proton_momentum * chain.direction / proton_tangent.Mag()) * proton_tangent;	
	  }
	}

	TVector2 mumom_vec_2d(muon_momentum_vec.X(), muon_momentum_vec.Y());
	TVector2 promom_vec_2d(proton_momentum_vec.X(), proton_momentum_vec.Y());
	TVector2 dpt_vec = mumom_vec_2d + promom_vec_2d;
	Double_t dalphat = std::acos((-1. * mumom_vec_2d * dpt_vec) / mumom_vec_2d.Mod() / dpt_vec.Mod()) * TMath::RadToDeg();
	double cosdat = (-1. * mumom_vec_2d * dpt_vec) / mumom_vec_2d.Mod() / dpt_vec.Mod();
	double dphit = std::acos((-1. * mumom_vec_2d * promom_vec_2d) / mumom_vec_2d.Mod() / promom_vec_2d.Mod()) * TMath::RadToDeg();
	double cosdphit = (-1. * mumom_vec_2d * promom_vec_2d) / mumom_vec_2d.Mod() / promom_vec_2d.Mod();
	double dptx = TMath::Sign(1., mumom_vec_2d.X()) * dpt_vec.Mod() * std::sin(dalphat);
	double dpty = dpt_vec.Mod() * cosdat;

	hist_dpt->Fill(dpt_vec.Mod(), ev.weight);
	hist_dalphat->Fill(dalphat, ev.weight);
	hist_dphit->Fill(dphit, ev.weight);
	hist_dptx->Fill(dptx, ev.weight);
	hist_dpty->Fill(dpty, ev.weight);

      }

      // 0pi2p
      if ( num_proton == 2 && ev.chains.size() == 3 ) {

	TVector3 proton_tangent_high;
	TVector3 proton_tangent_low;
	int proton_direction_high;
	int proton_direction_low;
	double proton_momentum_high;
	double proton_momentum_low;
	TVector3 proton_momentum_vec_high;
	TVector3 proton_momentum_vec_low;
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
	  if ( chain.particle_flag % 10000 == 2212 ) {
	    num_proton++;
	    if ( num_proton == 1 ) {
	      proton_tangent_high.SetXYZ(ax, ay, 1);
	      proton_direction_high = chain.direction;
	      if ( chain.stop_flag == 2 ) {
		proton_momentum_high = chain.ecc_range_mom[1];
	      }
	      else {
		proton_momentum_high = chain.ecc_mcs_mom[1];
	      }
	      proton_momentum_vec_high = (proton_momentum_high * chain.direction / proton_tangent_high.Mag()) * proton_tangent_high;
	    }
	    else if ( num_proton == 2 ) {
	      proton_tangent_low.SetXYZ(ax, ay, 1);
	      proton_direction_low = chain.direction;
	      if ( chain.stop_flag == 2 ) {
		proton_momentum_low = chain.ecc_range_mom[1];
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
	  std::swap(proton_momentum_low, proton_momentum_high);
	  std::swap(proton_momentum_vec_low, proton_momentum_vec_high);
	}

	// Calculate
	double open_cos = (proton_momentum_vec_high * proton_momentum_vec_low)
	  / proton_momentum_vec_high.Mag() / proton_momentum_vec_low.Mag();
	hist_open_cos->Fill(open_cos, ev.weight);
      }

    }
    else if ( ev.vertex_material == -3 ) {
      hist_side_recon_vertex_x->Fill(ev.recon_vertex_position[0], ev.weight);
      hist_side_recon_vertex_y->Fill(ev.recon_vertex_position[1], ev.weight);
      hist_side_recon_vertex_z->Fill(ev.recon_vertex_position[2], ev.weight);
      hist_side_recon_vertex_xy->Fill(ev.recon_vertex_position[0], ev.recon_vertex_position[1], ev.weight);
      hist_side_multi->Fill(ev.chains.size(), ev.weight);
      hist_side_multi_z->Fill(ev.recon_vertex_position[2], ev.chains.size(), ev.weight);
    }
    else if ( ev.vertex_material == -2 ) {
      hist_pene_recon_vertex_x->Fill(ev.recon_vertex_position[0], ev.weight);
      hist_pene_recon_vertex_y->Fill(ev.recon_vertex_position[1], ev.weight);
      hist_pene_recon_vertex_z->Fill(ev.recon_vertex_position[2], ev.weight);
      hist_pene_recon_vertex_xy->Fill(ev.recon_vertex_position[0], ev.recon_vertex_position[1], ev.weight);
      hist_pene_multi->Fill(ev.chains.size(), ev.weight);
      hist_pene_multi_z->Fill(ev.recon_vertex_position[2], ev.chains.size(), ev.weight);
    }
    
  }
  
  ofile->cd();
  hist_recon_vertex_x->Write();
  hist_recon_vertex_y->Write();
  hist_recon_vertex_z->Write();
  hist_recon_vertex_xy->Write();
  hist_vertex_pl->Write();
  hist_multi->Write();
  hist_multi_z->Write();
  hist_plate_diff->Write();

  hist_multi_p->Write();
  hist_multi_pi->Write();
  hist_multi_other->Write();

  hist_anu_bg_mu_mom->Write();
  hist_anu_bg_mu_mom_mcs->Write();
  hist_anu_bg_mu_mom_range->Write();
  hist_anu_bg_mu_cos->Write();
  hist_anu_bg_mu_deg->Write();

  hist_anu_bg_mu_mom_mcs_all->Write();
  
  hist_anu_bg_mu_mom_true->Write();
  hist_anu_bg_mu_mom_true_stop->Write();
  hist_anu_bg_mu_mom_true_nonstop->Write();
  
  hist_anu_bg_pi_mom->Write();
  hist_anu_bg_pi_mom_mcs->Write();
  hist_anu_bg_pi_mom_range->Write();
  hist_anu_bg_p_mom->Write();
  hist_anu_bg_p_mom_mcs->Write();
  hist_anu_bg_p_mom_range->Write();
  hist_anu_bg_pi_cos->Write();
  hist_anu_bg_pi_deg->Write();
  hist_anu_bg_p_cos->Write();
  hist_anu_bg_p_deg->Write();

  hist_anu_bg_other_mom_mcs->Write();
  hist_anu_bg_other_cos->Write();
  hist_anu_bg_other_deg->Write();

  hist_side_recon_vertex_x->Write();
  hist_side_recon_vertex_y->Write();
  hist_side_recon_vertex_z->Write();
  hist_side_recon_vertex_xy->Write();
  hist_side_multi->Write();
  hist_side_multi_z->Write();

  hist_pene_recon_vertex_x->Write();
  hist_pene_recon_vertex_y->Write();
  hist_pene_recon_vertex_z->Write();
  hist_pene_recon_vertex_xy->Write();
  hist_pene_multi->Write();
  hist_pene_multi_z->Write();

  hist_dpt->Write();
  hist_dalphat->Write();
  hist_dphit->Write();
  hist_dptx->Write();
  hist_dpty->Write();

  hist_open_cos->Write();

  ofile->Close();

  std::exit(0);

}
