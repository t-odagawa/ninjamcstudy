#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/filesystem.hpp>

#include <TFile.h>

#include <McsClass.hpp>

namespace logging = boost::log;
namespace fs = boost::filesystem;

void Analyze0pi1p(std::string b2filename,
		  std::string momchfilename,
		  std::string outpufilename) {

  BOOST_LOG_TRIVIAL(infl) << "==========0pi1p mode==========";

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
				 50, 0., 1500.);
  TH1D *hist_proton_mom = new TH1D("hist_proton_mom",
				   "Proton reconstructed momentum;p_{p} [MeV/c];Entries",
				   50, 0., 1500.);
  /*
  TH1D *hist_muon_ang;
  TH1D *hist_muon_cos;
  TH1D *hist_proton_ang;
  TH1D *hist_proton_cos;
  TH2D *hist_muon_mom_ang = new TH2D("hist_muon_mom_ang",
				     "Muon momentum;p_{#mu, true} [MeV/c];#theta_{#mu} [degree]",
				     15, 0., 1500., 15, 0., 1500.);
  TH2D *hist_proton_mom_ang = new TH1D("hist_proton_mom_ang");
  TH1D *hist_nu_ene_bias;
  TH1D *hist_nu_ene_recon;
  TH2D *hist_nu_ene_recon_true;
  */
  TH1D *hist_dpt = new TH1D("hist_dpt",
			    "#deltap_{T};#deltap_{T} [MeV/c];Entries",
			    200, 0., 2.);
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
  TH1D *hist_mode_dpt[num_ninja_mode];
  TH1D *hist_mode_dalphat[num_ninja_mode];
  TH1D *hist_mode_cosdat[num_ninja_mode];
  TH1D *hist_mode_dphit[num_ninja_mode];
  TH1D *hist_mode_cosdphit[num_ninja_mode];
  TH1D *hist_mode_dptx[num_ninja_mode];
  TH1D *hist_mode_dpty[num_ninja_mode];
  for ( int i = 0; i < num_ninja_mode; i++ ) {
    hist_mode_dpt[i] = new TH1D(Form("hist_dpt_%d", i), "", 200, 0., 2.);
    hist_mode_dpt[i]->SetFillColor(mode_color[i]);
    hist_mode_dpt[i]->SetFillStyle(mode_style[i]);
    hist_mode_dalphat[i] = new TH1D(Form("hist_dalphat_%d", i), 90, 0., 180.);
    hist_mode_dalphat[i]->SetFillColor(mode_color[i]);
    hist_mode_dalphat[i]->SetFillStyle(mode_style[i]);
    hist_mode_cosdat[i] = new TH1D(Form("hist_cosdat_%d", i), 100, -1., 1.);
    hist_mode_cosdat[i]->SetFillColor(mode_color[i]);
    hist_mode_cosdat[i]->SetFillStyle(mode_style[i]);
    hist_mode_dphit[i] = new TH1D(Form("hist_dphit_%d", i), 90, 0., 180.);
    hist_mode_dphit[i]->SetFillColor(mode_color[i]);
    hist_mode_dphit[i]->SetFillStyle(mode_style[i]);
    hist_mode_cosdphit[i] = new TH1D(Form("hist_cosdphit_%d", i), 100, -1., 1.);
    hist_mode_cosdphit[i]->SetFillColor(mode_color[i]);
    hist_mode_cosdphit[i]->SetFillStyle(mode_style[i]);
    hist_mode_dptx[i] = new TH1D(Form("hist_dptx_%d", i), 200, -2., 2.);
    hist_mode_dptx[i]->SetFillColor(mode_color[i]);
    hist_mode_dptx[i]->SetFillStyle(mode_style[i]);
    hist_mode_dpty[i] = new TH1D(Form("hist_dpty_%d", i), 200, -2., 2.);
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
      
      TVector3 muon_tangent;
      TVector3 proton_tangent;
      double muon_momentum;
      double proton_momentum;
      TVector3 muon_momentum_vec;
      TVector3 proton_momentum_vec;

      for ( auto chain : ev.chains ) {
	int particle_id = chain.particle_flag % 10000;
	if ( particle_id == 13 ) {
	  num_muon++;
	  muon_tangent.SetXYZ(chain.base.back().ax,
			      chain.base.back().ay,
			      1.);
	  muon_momentum = chain.ecc_mcs_mom[0]; // Baby MIND range?
	  muon_momentum_vec = (muon_momentum / muon_tangent.Mag()) * muon_tangent;
	}
	else if ( particle_id == 2212 ) {
	  num_proton++;
	  if ( chain.direction == 1 ) {
	    proton_tangent.SetXYZ(chain.base.back().ax,
				  chain.base.back().ay,
				  1.);
	  }
	  else if ( chain.direction == -1 ) {
	    proton_tangent.SetXYZ(chain.base.front().ax,
				  chain.base.front().ay,
				  1.);
	  }
	  if ( chain.stop_flag == 2 ) {
	    proton_momentum = chain.ecc_range_mom[1];
	  }
	  else {
	    proton_momentum = chain.ecc_mcs_mom[1];
	  }
	  proton_momentum_vec = (proton_momentum * proton_tangent.Z() / proton_tangent.Mag()) * proton_tangent;
	}
      }

      if ( num_muon != 1 && num_proton != 1 ) continue;

      TVector2 dpt_vec(muon_momentum_vec.X() + proton_momentum_vec.X(),
		       muon_momentum_vec.Y() + proton_momentum_vec.Y());
      TVector2 mumom_vec_2d(muon_momentum_vec.X(), muon_momentum_vec.Y());
      TVector2 promom_vec_2d(proton_momentum_vec.X(), proton_momentum_vec.Y());
      Double_t dalphat = std::acos((-mumom_vec_2d * dpt_vec) / mumom_vec_2d.Mag() / dpt_vec.Mag());
      hist_dpt->Fill(dpt_vec.Mag(), ev.weight);
      hist_mode_dpt[mode_id]->Fill(dpt_vec.Mag(), ev.weight);
      hist_dalphat->Fill(dalphat * TMath::RadToDeg(),
			 ev.weight);
      hist_mode_dalphat[mode_id]->Fill(dalphat * TMath::RadToDeg(),
				  ev.weight);
      hist_cosdat->Fill((-mumom_vec_2d * dpt_vec) / mumom_vec_2d.Mag() / dpt_vec.Mag(), ev.weight);
      hist_mode_cosdat[mode_id]->Fill((-mumom_vec_2d * dpt_vec) / mumom_vec_2d.Mag() / dpt_vec.Mag(), ev.weight);
      hist_dphit->Fill(std::acos((-mumom_vec_2d * promom_vec_2d) / mumom_vec_2d.Mag() / promom_vec_2d.Mag()) * TMath::RadToDeg(),
		       ev.weight);
      hist_mode_dphit[mode_id]->Fill(std::acos((-mumom_vec_2d * promom_vec_2d) / mumom_vec_2d.Mag() / promom_vec_2d.Mag()) * TMath::RadToDeg(),
				     ev.weight);
      hist_cosdphit->Fill((-mumom_vec_2d * promom_vec_2d) / mumom_vec_2d.Mag() / promom_vec_2d.Mag(),
			  ev.weight);
      hist_mode_cosdphit[mode_id]->Fill((-mumom_vec_2d * promom_vec_2d) / mumom_vec_2d.Mag() / promom_vec_2d.Mag(),
					ev.weight);
      hist_dptx->Fill(TMath::Sign(1., mumom_vec_2d.X()) * dpt_vec.Mag() * std::sin(dalphat),
		      ev.weight);
      hist_mode_dptx[mode_id]->Fill(TMath::Sign(1., mumom_vec_2d.X()) * dpt_vec.Mag() * std::sin(dalphat),
				    ev.weight);
      hist_dpty->Fill(dpt_vec.Mag() * std::cos(dalphat), ev.weight);
      hist_mode_dpty[mode_id]->Fill(dpt_vec.Mag() * std::cos(dalphat), ev.weight);
    }
 
  }

  outputfile->cd();
  hist_muon_mom->Write();
  hist_proton_mom->Write();
  /*
  hist_muon_ang->Write();
  hist_muon_cos->Write();
  hist_proton_ang->Write();
  hist_proton_cos->Write();
  hist_muon_mom_ang->Write();
  hist_proton_mom_ang->Write();
  hist_nu_ene_bias->Write();
  hist_nu_ene_recon->Write();
  hist_nu_ene_recon_true->Write();
  */
  hist_dpt->Write();
  hist_dalphat->Write();
  hist_cosdat->Write();
  hist_dphit->Write();
  hist_cosdphit->Write();
  hist_dptx->Write();
  hist_dpty->Write();
  for ( int i = 0; i < num_ninja_mode; i++ ) {
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
