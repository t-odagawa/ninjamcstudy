#include <cmath>
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/filesystem.hpp>

#include <McsClass.hpp>

int main(int argc, char* argv[]) {

  if ( argc != 3 ) {
    BOOST_LOG_TRIVIAL(error) << "Usage : " << argv[0]
			     << " <input momch file> <output file>";
  }

  std::string momchfilename = argv[1];
  if ( !fs::exists(momchfilename) ) {
    throw std::runtime_error("File not found : " + momchfilename);
  }
  auto ev_vec = Momentum_recon::ReadEventInformationBin(momchfilename);

  TFile *outputfile = new TFile((TString)argv[2], "recreate");

  TH1D *hist_muon_mom = new TH1D("hist_muon_mom", "", 20, 0., 2000.);
  TH1D *hist_muon_ang = new TH1D("hist_muon_ang", "", 18, 0., 90.);
  TH1D *hist_muon_cos = new TH1D("hist_muon_cos", "", 20, 0., 1.);

  TH1D *hist_proton_mom = new TH1D("hist_proton_mom", "", 15, 0., 1500.);
  TH1D *hist_proton_ang = new TH1D("hist_proton_ang", "", 36, 0., 180.);
  TH1D *hist_proton_cos = new TH1D("hist_proton_cos", "", 40, -1., 1.);

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

  TH1D *hist_muon_mom_2 = new TH1D("hist_muon_mom_2", "", 20, 0., 2000.);
  TH1D *hist_muon_ang_2 = new TH1D("hist_muon_ang_2", "", 18, 0., 90.);
  TH1D *hist_muon_cos_2 = new TH1D("hist_muon_cos_2", "", 20, 0., 1.);

  TH1D *hist_proton_mom_2 = new TH1D("hist_proton_mom_2", "", 15, 0., 1500.);
  TH1D *hist_proton_ang_2 = new TH1D("hist_proton_ang_2", "", 36, 0., 180.);
  TH1D *hist_proton_cos_2 = new TH1D("hist_proton_cos_2", "", 40, -1., 1.);

  TH1D *hist_proton_mom_high = new TH1D("hist_proton_mom_high",
					"",
					15, 0., 1500.);
  TH1D *hist_proton_mom_low = new TH1D("hist_proton_mom_low",
				       "",
				       15, 0., 1500.);
  TH1D *hist_proton_ang_high = new TH1D("hist_proton_ang_high",
					"",
					36, 0., 180.);
  TH1D *hist_proton_ang_low = new TH1D("hist_proton_ang_low",
				       "",
				       36, 0., 180.);
  TH1D *hist_proton_cos_high = new TH1D("hist_proton_cos_high",
					"", 
					40, -1., 1.);
  TH1D *hist_proton_cos_low = new TH1D("hist_proton_cos_low",
				       "",
				       40, -1., 1.);

  TH1D *hist_open_ang = new TH1D("hist_open_ang",
				 "",
				 36, 0., 180.);
  TH1D *hist_open_cos = new TH1D("hist_open_cos",
				 "",
				 40, -1., 1.);
  TH1D *hist_mom_ratio = new TH1D("hist_mom_ratio",
				  "",
				  10, 0., 0.5);
  TH1D *hist_mom_vecsum = new TH1D("hist_mom_vecsum",
				   "",
				   100, 0., 2000.);
  TH1D *hist_mom_scasum = new TH1D("hist_mom_scasum",
				   "",				  
				   100, 0., 2000.);

  TH1D *hist_dptt = new TH1D("hist_dptt",
			     "",
			     100, -500., 500.);
  TH1D *hist_dpt_2 = new TH1D("hist_dpt_2",
			      "",
			    200, 0., 2000.);
  TH1D *hist_pn = new TH1D("hist_pn",
			   "",
			   200, 0., 2000.);
  TH1D *hist_dalphat_2 = new TH1D("hist_dalphat_2",
				  "",
				  36, 0., 180.);
  TH1D *hist_cosdat_2 = new TH1D("hist_cosdat_2",
			       "",
			       40, -1., 1.);

  for ( auto ev : ev_vec ) {

    if ( ev.vertex_material != B2Material::kWater ) continue;

    if ( ev.chains.empty() ) continue;

    int num_muon = 0;
    int num_proton = 0;
    int num_pion = 0;
    
    for ( auto chain : ev.chains ) {
      int recon_particle_id = chain.particle_flag % 10000;
      int true_particle_id = chain.particle_flag / 10000;
      if ( recon_particle_id == 13 ) num_muon++;
      else if ( recon_particle_id == 2212 ) num_proton++;
      else if ( recon_particle_id == 211 ) num_pion++;
    }

    if ( ev.chains.size() == 2 &&
	 num_muon == 1 && num_proton == 1 && num_pion == 0 ) {

      TVector3 muon_tangent;
      double muon_momentum;
      TVector3 proton_tangent;
      double proton_momentum;
      TVector3 proton_momentum_vec;
      int proton_direction;

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

	ax = std::tan(thetax); ay = std::atn(thetay);

	int recon_particle_id = chain.particle_flag % 10000;
	if ( recon_particle_id == 13 ) {
	  if ( chain.stop_flag == 0 ) 
	    muon_momentum = chain.ecc_mcs_mom[0];
	  else if ( chain.stop_flag == 1 )
	    muon_momentum = chain.bm_range_mom;

	  muon_tangent.SetXYZ(ax, ay, 1.);
	  muon_momentum_vec = (muon_momentum / muon_tangent.Mag()) * muon_tangent;
	}
	else if ( recon_particle_id == 2212 ) {
	  if ( chain.direction == 1 ) {
	    proton_tangent.SetXYZ(ax, ay, 1);
	    proton_direction = chain.direction;
	  }
	  else if ( chain.direction == -1 ) {
	    proton_tangent.SetXYZ(ax, ay, 1);
	    proton_direction = chain.direction;
	  }

	  if ( chain.stop_flag == 2 )
	    proton_momentum = chain.ecc_range_mom[1];	   
	  else if ( chain.stop_flag == 0 ) 
	    proton_momentum = chain.ecc_mcs_mom[1];

	  proton_momentum_vec = (proton_momentum * chain.direction / proton_tangent.Mag()) * proton_tangent;
	}
      }

      double muon_rad = std::atan(std::hypot(muon_tangent.X(), muon_tangent.Y()));
      double muon_ang = muon_rad * TMath::RadToDeg();

      hist_muon_mom->Fill(muon_momentum, ev.weight);
      hist_muon_ang->Fill(muon_ang, ev.weight);
      hist_muon_cos->Fill(std::cos(muon_ang), ev.weight);

      double proton_rad = std::atan(proton_direction * std::hypot(proton_tangent.X(), proton_tangent.Y()));
      double proton_ang = proton_rad * TMath::RadToDeg();
      if ( proton_direction == -1 )
	proton_ang += 180.;
      
      hist_proton_mom->Fill(proton_momentum, ev.weight);
      hist_proton_ang->Fill(proton_ang, ev.weight);
      hist_proton_cos->Fill(std::cos(proton_ang * TMath::DegToRad()), ev.weight);

      TVector2 mumom_vec_2d(muon_momentum_vec.X(), muon_momentum_vec.Y());
      TVector2 promom_vec_2d(proton_momentum_vec.X(), proton_momentum_vec.Y());
      TVector2 dpt_vec = mumom_vec_2d + promom_vec_2d;
      double dalphat = std::acos((-1. * mumom_vec_2d * dpt_vec) / mumom_vec_2d.Mod() / dpt_vec.Mod());
      hist_dpt->Fiill(dpt_vec.Mod(), ev.weight);
      hist_dalphat->Fill(dalphat * TMath::RadToDeg(), ev.weight);
      hist_cosdat->Fill((-1. * mumom_vec_2d * dpt_vec) / mumom_vec_2d.Mod() / dpt_vec.Mod(), ev.weight);
      hist_dphit->Fill(std::acos((-1. * mumom_vec_2d * promom_vec_2d) / mumom_vec_2d.Mod() / promom_vec_2d.Mod()) * TMath::RadToDeg(), ev.weight);
      hist_cosdphit->Fill((-1. * mumom_vec_2d * promom_vec_2d) / mumom_vec_2d.Mod() / promom_vec_2d.Mod(), ev.weight);
      hist_dptx->Fill(TMath::Sign(1., mumom_vec_2d.X()) * dpt_vec.Mod() * std::sin(dalphat), ev.weight);
      hist_dpty->Fill(dpt_vec.Mod() * std::cos(dalphat), ev.weight);

    }
    else if ( ev.chains.size() == 3 &&
	      num_muon == 1 && num_proton == 2 && num_pion == 0 ) {

      TVector3 muon_tangent;
      double muon_momentum;
      TVector3 muon_tangent_vec;
      TVector3 proton_tangent_tmp;
      TVector3 proton_tangent_high;
      TVector3 proton_tangent_low;
      double proton_momentum_tmp;
      double proton_momentum_high;
      double proton_momentum_low;
      TVector3 proton_momentum_vec_tmp;
      TVector3 proton_momentum_vec_high;
      TVector3 proton_momentum_vec_low;
      int proton_direction_tmp;
      int proton_direction_high;
      int proton_direction_low;

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

	int recon_particle_id = chain.particle_flag % 10000;

	int num_proton = 0;

	if ( recon_particle_id == 13 ) {
	  if ( chain.stop_flag == 1 ) 
	    muon_momentum = chain.bm_range_mom;
	  else if ( chain.stop_flag == 0 ) 
	    muon_momentum = chain.ecc_mcs_mom[0];

	  muon_tangent.SetXYZ(ax, ay, 1.);
	  muon_momentum_vec = (muon_momentum / muon_tangent.Mag()) * muon_tangent;	  
	}
	else if ( recon_particle_id == 2212 ) {
	  num_proton++;
	  if ( chain.stop_flag == 2 )
	    proton_momentum_tmp = chain.ecc_range_mom[1];
	  else if ( chain.stop_flag == 0 )
	    proton_momentum_tmp = chain.ecc_mcs_mom[1];

	  if ( chain.direction == 1 ) {
	    proton_tangent_tmp.SetXYZ(ax, ay, 1.);
	    proton_direction_tmp = 1;
	  }
	  else if ( chain.direction == -1 ) {
	    proton_tangent_tmp.SetXYZ(ax, ay, 1);
	    proton_direction_tmp = -1;
	  }

	  proton_momentum_vec_tmp = (proton_momentum_tmp * proton_direction_tmp / proton_tangent_tmp.Mag()) * proton_tangent_tmp;

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
      }

      if ( proton_momentum_low > proton_momentum_high ) {
	std::swap(proton_tangent_high, proton_tangent_low);
	std::swap(proton_direction_high, proton_direction_low);
	std::swap(proton_momentum_high, proton_momentum_low);
	std::swap(proton_momentum_vec_high, proton_momentum_vec_low);
      }

      double muon_rad = std::atan(std::hypot(muon_tangent.X(), muon_tangent.Y()));
      double muon_ang = muon_rad * TMath::RadToDeg();
      
      hist_muon_mom_2->Fill(muon_momentum, ev.weight);
      hist_muon_ang_2->Fill(muon_ang, ev.weight);
      hist_muon_cos_2->Fill(std::cos(muon_rad), ev.weight);

      double proton_rad_high = std::atan(proton_direction_high
					 * std::hypot(proton_tangent_high.X(), proton_tangent_high.Y()));
      double proton_ang_high = proton_rad_high * TMath::RadToDeg();
      if ( proton_direction_high == -1 ) proton_ang_high += 180.;

      double proton_rad_low = std::atan(proton_direction_low
					* std::hypot(proton_tangent_low.X(), proton_tangent_low.Y()));
      double proton_ang_low = proton_rad_low * TMath::RadToDeg();
      if ( proton_direction_low == -1 ) proton_ang_low += 180.;

      hist_proton_mom_2->Fill(proton_mom_high, ev.weight);
      hist_proton_mom_2->Fill(proton_mom_low, ev.weight);
      hist_proton_mom_high->Fill(proton_mom_high, ev.weight);
      hist_proton_mom_low->Fill(proton_mom_low, ev.weight);
      
      hist_proton_ang_2->Fill(proton_ang_high, ev.weight);
      hist_proton_ang_2->Fill(proton_ang_low, ev.weight);
      hist_proton_ang_high->Fill(proton_ang_high, ev.weight);
      hist_proton_ang_low->Fill(proton_ang_low, ev.weight);

      hist_proton_cos_2->Fill(std::cos(proton_ang_high * TMath::DegToRad()), ev.weight);
      hist_proton_cos_2->Fill(std::cos(proton_ang_low * TMath::DegToRad()), ev.weight);
      hist_proton_cos_high->Fill(std::cos(proton_ang_high * TMath::DegToRad()), ev,weight);
      hist_proton_cos_low->Fill(std::cos(proton_ang_low * TMath::DegToRad()), ev.weight);

      double open_cos = (proton_momentum_vec_high * proton_momentum_vec_low)
	/ proton_momentum_vec_high.Mag() / proton_momentum_vec_low.Mag();
      double open_ang = std::acos(open_cos) * TMath::RadToDeg();

      hist_open_ang->Fill(open_ang, ev.weight);
      hist_open_cos->Fill(open_cos, ev.weight);
      hist_mom_ratio->Fill(proton_momentum_low / proton_momentum_high, ev.weight);
      hist_mom_vecsum->Fill((proton_momentum_vec_high + proton_momentum_vec_low).Mag(), ev.weight);
      hist_mom_scasum->Fill(proton_momentum_vec_high.Mag() + proton_momentum_vec_low.Mag(), ev.weight);

      TVector3 ztt = (-1. * muon_momentum_vec.Y(), muon_momentum_Vec.X(), 0.);
      ztt = (1. / std::hypot(muon_momentum_vec.X(), muon_momentum_vec.Y())) * ztt;
      double dptt = ztt * proton_momentum_vec_high + ztt * proton_momentum_vec_low;
      hist_dptt->Fill(dptt, ev.weight);

      TVector2 mumom_vec_2d(muon_momentum_vec.X(), muon_momentum_vec.Y());
      TVector2 promom_vec_2d_high(proton_momentum_vec_high.X(), proton_momentum_vec_high.Y());
      TVector2 promom_vec_2d_low(proton_momentum_vec_low.X(), proton_momentum_vec_low.Y());
      TVector2 dpt_vec = mumom_vec_2d + promom_vec_2d_high + promom_vec_2d_low;
      hist_dpt_2->Fill(dpt_vec.Mod(), ev.weight);

      double muon_energy = std::sqrt(muon_momentum * muon_momentum + muon_mass * muon_mass);
      double proton_energy_high = std::sqrt(proton_momentum_high * proton_momentum_high
					    + proton_mass * proton_mass);
      double proton_energy_low = std::sqrt(proton_momentum_low * proton_momentum_low
					   + proton_mass * proton_mass);
      double residual_mass = 23.;
      double pl = 0.5 * (proton_mass + muon_momentum_vec.Z() + proton_momentum_vec_high.Z() + proton_momentum_vec_low.Z()
			 - muon_energy - proton_energy_high - proton_energy_low)
	- 0.5 * (dpt_vec.Mod2() + residual_mass * residual_mass)
	/ (proton_mass + muon_momentum_vec.Z() + proton_momentum_vec_high.Z() + proton_momentum_vec_low.Z()
	   - muon_energy - proton_energy_high - proton_energy_low);
      double pn = std::sqrt(dpt_vec.Mod2() + pl * pl);
      hist_pn->Fill(pn, ev.weight);

      double dalphat = std::acos((-1. * mumom_vec_2d * dpt_vec) / mumom_vec_2d.Mod() / dpt_vec.Mod());
      hist_dalphat_2->Fill(dalphat, ev.weight);
      hist_cosdat_2->Fill((-1. * mumom_vec_2d * dpt_vec) / mumom_vec_2d.Mod() / dpt_vec.Mod(), ev.weight);

    }
    
  }


  outputfile->cd();
  hist_muon_mom->Write();
  hist_muon_ang->Write();
  hist_muon_cos->Write();
  
  hist_proton_mom->Write();
  hist_proton_ang->Write();
  hist_proton_cos->Write();

  hist_dpt->Write();
  hist_dalphat->Write();
  hist_cosdat->Write();
  hist_dphit->Write();
  hist_cosdphit->Write();
  hist_dptx->Write();
  hist_dpty->Write();

  hist_muon_mom_2->Write();
  hist_muon_ang_2->Write();
  hist_muon_cos_2->Write();
  hist_proton_mom_2->Write();
  hist_proton_ang_2->Write();
  hist_proton_cos_2->Write();
  hist_proton_mom_high->Write();
  hist_proton_mom_low->Write();
  hist_proton_ang_high->Write();
  hist_proton_ang_low->Write();
  hist_proton_cos_high->Write();
  hist_proton_cos_low->Write();

  hist_open_ang->Write();
  hist_open_cos->Write();
  hist_mom_ratio->Write();
  hist_mom_vecsum->Write();
  hist_mom_scasum->Write();
  hist_dptt->Write();
  hist_dpt_2->Write();
  hist_pn->Write();
  hist_dalphat_2->Write();
  hist_cosdat_2->Write();

  ofile->Close();

}
