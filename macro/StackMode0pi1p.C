#include "/home/t2k/odagawa/NinjaMCStudy/lib/HistogramStyle.hpp"

void StackMode0pi1p() {

  TString filename = "/hsm/nu/ninja/pra_tmp/mc_tmp_20220620/output/output_mode6.root";
  TFile *file = new TFile(filename, "read");

  THStack *hs_muon_mom = new THStack("hs_muon_mom", "Muon momentum;p_{#mu} [MeV/c];Entries");
  THStack *hs_proton_mom = new THStack("hs_proton_mom", "Proton momentum;p_{p} [MeV/c];Entries");  
  THStack *hs_muon_ang = new THStack("hs_muon_ang", "Muon angle;#theta_{#mu} [degree];Entries");
  THStack *hs_proton_ang = new THStack("hs_proton_ang", "Proton angle;#theta_{p} [degree];Entries");
  THStack *hs_muon_cos = new THStack("hs_muon_cos", "Muon angle;cos#theta_{#mu};Entries");
  THStack *hs_proton_cos = new THStack("hs_proton_cos", "Proton angle;cos#theta_{p}; Entries");

  THStack *hs_q2 = new THStack("hs_q2", "Q^{2};Q^{2} [MeV^{2}/c^{2}];Entries");
  THStack *hs_nu_ene_bias = new THStack("hs_nu_ene_bias", "Energy bias;(E_{#nu,recon} - E_{#nu,true}) / E_{#nu,true};Entries");
  THStack *hs_nu_ene_recon = new THStack("hs_nu_ene_recon", "Reconstructed neutrino energy;E_{#nu,recon} [MeV];Entries");

  THStack *hs_dpt = new THStack("hs_dpt", "#deltap_{T};#deltap_{T} [MeV/c];Entries");
  THStack *hs_dalphat = new THStack("hs_dalphat", "#delta#alpha_{T};#delta#alpha_{T} [degree];Entries");
  THStack *hs_cosdat = new THStack("hs_cosdat", "cos#delta#alpha_{T};cos#delta#alpha_{T};Entries");
  THStack *hs_dphit = new THStack("hs_dphit", "#delta#phi_{T};#delta#phi_{T} [degree];Entries");
  THStack *hs_cosdphit = new THStack("hs_cosdphit", "cos#delta#phi_{T};cos#delta#phi_{T};Entries");
  THStack *hs_dptx = new THStack("hs_dptx", "#deltap_{Tx};#deltap_{Tx} [MeV/c];Entries");
  THStack *hs_dpty = new THStack("hs_dpty", "#deltap_{Ty};#deltap_{Ty} [MeV/c];Entries");

  TLegend *l = new TLegend(0.6, 0.5, 0.85, 0.85);
  l->SetName("l");

  TH1D *hist_muon_mom[num_ninja_mode];
  TH1D *hist_proton_mom[num_ninja_mode];
  TH1D *hist_muon_ang[num_ninja_mode];
  TH1D *hist_proton_ang[num_ninja_mode];
  TH1D *hist_muon_cos[num_ninja_mode];
  TH1D *hist_proton_cos[num_ninja_mode];

  TH1D *hist_q2[num_ninja_mode];
  TH1D *hist_nu_ene_bias[num_ninja_mode];
  TH1D *hist_nu_ene_recon[num_ninja_mode];

  TH1D *hist_dpt[num_ninja_mode];
  TH1D *hist_dalphat[num_ninja_mode];
  TH1D *hist_cosdat[num_ninja_mode];
  TH1D *hist_dphit[num_ninja_mode];
  TH1D *hist_cosdphit[num_ninja_mode];
  TH1D *hist_dptx[num_ninja_mode];
  TH1D *hist_dpty[num_ninja_mode];

  for ( int i = 0; i < num_ninja_mode; i++ ) {
    hist_muon_mom[i] = (TH1D*)file->Get(Form("hist_muon_mom_%d", i));
    hist_proton_mom[i] = (TH1D*)file->Get(Form("hist_proton_mom_%d", i));
    hist_muon_ang[i] = (TH1D*)file->Get(Form("hist_muon_ang_%d", i));
    hist_proton_ang[i] = (TH1D*)file->Get(Form("hist_proton_ang_%d", i));
    hist_muon_cos[i] = (TH1D*)file->Get(Form("hist_muon_cos_%d", i));
    hist_proton_cos[i] = (TH1D*)file->Get(Form("hist_proton_cos_%d", i));

    hist_q2[i] = (TH1D*)file->Get(Form("hist_q2_%d", i));
    hist_nu_ene_bias[i] = (TH1D*)file->Get(Form("hist_nu_ene_bias_%d", i));
    hist_nu_ene_recon[i] = (TH1D*)file->Get(Form("hist_nu_ene_recon_%d", i));

    hist_dpt[i] = (TH1D*)file->Get(Form("hist_dpt_%d", i));
    hist_dalphat[i] = (TH1D*)file->Get(Form("hist_dalphat_%d", i));
    hist_cosdat[i] = (TH1D*)file->Get(Form("hist_cosdat_%d", i));
    hist_dphit[i] = (TH1D*)file->Get(Form("hist_dphit_%d", i));
    hist_cosdphit[i] = (TH1D*)file->Get(Form("hist_cosdphit_%d", i));
    hist_dptx[i] = (TH1D*)file->Get(Form("hist_dptx_%d", i));
    hist_dpty[i] = (TH1D*)file->Get(Form("hist_dpty_%d", i));

    l->AddEntry(hist_muon_mom[i], mode_name[i], "f");
  }

  for ( int i = 0; i < num_ninja_mode; i++ ) {
    hs_muon_mom->Add(hist_muon_mom[mode_stack_order[i]]);
    hs_proton_mom->Add(hist_proton_mom[mode_stack_order[i]]);
    hs_muon_ang->Add(hist_muon_ang[mode_stack_order[i]]);
    hs_proton_ang->Add(hist_proton_ang[mode_stack_order[i]]);
    hs_muon_cos->Add(hist_muon_cos[mode_stack_order[i]]);
    hs_proton_cos->Add(hist_proton_cos[mode_stack_order[i]]);

    hs_q2->Add(hist_q2[mode_stack_order[i]]);
    hs_nu_ene_bias->Add(hist_nu_ene_bias[mode_stack_order[i]]);
    hs_nu_ene_recon->Add(hist_nu_ene_recon[mode_stack_order[i]]);

    hs_dpt->Add(hist_dpt[mode_stack_order[i]]);
    hs_dalphat->Add(hist_dalphat[mode_stack_order[i]]);
    hs_cosdat->Add(hist_cosdat[mode_stack_order[i]]);
    hs_dphit->Add(hist_dphit[mode_stack_order[i]]);
    hs_cosdphit->Add(hist_cosdphit[mode_stack_order[i]]);
    hs_dptx->Add(hist_dptx[mode_stack_order[i]]);
    hs_dpty->Add(hist_dpty[mode_stack_order[i]]);
  }


  TFile *ofile = new TFile("~/stack_0pi1p.root", "recreate");
  
  ofile->cd();
  hs_muon_mom->Write();
  hs_proton_mom->Write();
  hs_muon_ang->Write();
  hs_proton_ang->Write();
  hs_muon_cos->Write();
  hs_proton_cos->Write();

  hs_q2->Write();
  hs_nu_ene_bias->Write();
  hs_nu_ene_recon->Write();

  hs_dpt->Write();
  hs_dalphat->Write();
  hs_cosdat->Write();
  hs_dphit->Write();
  hs_cosdphit->Write();
  hs_dptx->Write();
  hs_dpty->Write();

  l->Write();
  ofile->Close();
}
