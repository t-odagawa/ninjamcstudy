#include "/home/t2k/odagawa/NinjaMCStudy/lib/HistogramStyle.hpp"
void StackModeAng() {

  // Legend 
  TLegend *l_muon_ang = new TLegend(0.5, 0.45, 0.85, 0.85);
  l_muon_ang->SetName("l_muon_ang");
  TLegend *l_proton_ang = new TLegend(0.5, 0.45, 0.85, 0.85);
  l_proton_ang->SetName("l_proton_ang");

  // Signal + packing/mis-pid background
  TString filename = "/hsm/nu/ninja/pra_tmp/CC0pi_20221213/signal/output/output_mode3.root";
  TFile *file = new TFile(filename, "read");
  double scale = 1. / 976. / 33.156 * 0.47 * 0.99;

  // Signal
  TH1D *hist_muon_deg[num_ninja_mode];
  TH1D *hist_muon_cos[num_ninja_mode];
  TH1D *hist_pion_deg[num_ninja_mode];
  TH1D *hist_pion_cos[num_ninja_mode];
  TH1D *hist_proton_deg[num_ninja_mode];
  TH1D *hist_proton_cos[num_ninja_mode];
  for ( int i = 0; i < num_ninja_mode - 1; i++ ) {
    hist_muon_deg[i] = (TH1D*)file->Get(Form("hist_muon_ang_deg_%d", i));
    hist_muon_cos[i] = (TH1D*)file->Get(Form("hist_muon_ang_cos_%d", i));
    hist_pion_deg[i] = (TH1D*)file->Get(Form("hist_pion_ang_deg_%d", i));
    hist_pion_cos[i] = (TH1D*)file->Get(Form("hist_pion_ang_cos_%d", i));
    hist_proton_deg[i] = (TH1D*)file->Get(Form("hist_proton_ang_deg_%d", i));
    hist_proton_cos[i] = (TH1D*)file->Get(Form("hist_proton_ang_cos_%d", i));
    hist_muon_deg[i]->Scale(scale);
    hist_muon_cos[i]->Scale(scale);
    hist_pion_deg[i]->Scale(scale);
    hist_pion_cos[i]->Scale(scale);
    hist_proton_deg[i]->Scale(scale);
    hist_proton_cos[i]->Scale(scale);
    l_muon_ang->AddEntry(hist_muon_deg[i], mode_name[i], "f");
    l_proton_ang->AddEntry(hist_proton_deg[i], mode_name[i], "f");
  }

  // Packing background
  auto hist_muon_deg_pack_bg = (TH1D*)file->Get("hist_muon_ang_deg_single");
  auto hist_muon_cos_pack_bg = (TH1D*)file->Get("hist_muon_ang_cos_single");
  hist_muon_deg_pack_bg->Scale(scale * 0.09);
  hist_muon_cos_pack_bg->Scale(scale * 0.09);

  // Muon mis-id CC + NC
  auto hist_muon_misid_ang_deg = (TH1D*)file->Get("hist_muon_misid_ang_deg");
  auto hist_muon_misid_ang_cos = (TH1D*)file->Get("hist_muon_misid_ang_cos");
  auto hist_pion_ang_deg_muon_misid = (TH1D*)file->Get("hist_pion_ang_deg_muon_misid");
  auto hist_pion_ang_cos_muon_misid = (TH1D*)file->Get("hist_pion_ang_cos_muon_misid");
  auto hist_proton_ang_deg_muon_misid = (TH1D*)file->Get("hist_proton_ang_deg_muon_misid");
  auto hist_proton_ang_cos_muon_misid = (TH1D*)file->Get("hist_proton_ang_cos_muon_misid");
  hist_muon_misid_ang_deg->Scale(scale);
  hist_muon_misid_ang_cos->Scale(scale);
  hist_pion_ang_deg_muon_misid->Scale(scale);
  hist_pion_ang_cos_muon_misid->Scale(scale);
  hist_proton_ang_deg_muon_misid->Scale(scale);
  hist_proton_ang_cos_muon_misid->Scale(scale);

  // Partner mis-id
  auto hist_pion_misid_ang_deg = (TH1D*)file->Get("hist_pion_misid_ang_deg");
  auto hist_pion_misid_ang_cos = (TH1D*)file->Get("hist_pion_misid_ang_cos");
  auto hist_proton_misid_ang_deg = (TH1D*)file->Get("hist_proton_misid_ang_deg");
  auto hist_proton_misid_ang_cos = (TH1D*)file->Get("hist_proton_misid_ang_cos");
  hist_pion_misid_ang_deg->Scale(scale);
  hist_pion_misid_ang_deg->SetFillColor(kOrange);
  hist_pion_misid_ang_cos->Scale(scale);
  hist_pion_misid_ang_cos->SetFillColor(kOrange);
  hist_proton_misid_ang_deg->Scale(scale);
  hist_proton_misid_ang_deg->SetFillColor(kOrange);
  hist_proton_misid_ang_cos->Scale(scale);
  hist_proton_misid_ang_cos->SetFillColor(kOrange);


  // Iron (+emulsion) interaction background
  TString fefilename = "/hsm/nu/ninja/pra_tmp/CC0pi_20221213/fe/output/bg_dist.root";
  TFile *fefile = new TFile(fefilename, "read");
  auto hist_fe_mu_ang_deg = (TH1D*)fefile->Get("hist_fe_bg_mu_deg");
  auto hist_fe_mu_ang_cos = (TH1D*)fefile->Get("hist_fe_bg_mu_cos");
  auto hist_fe_pi_ang_deg = (TH1D*)fefile->Get("hist_fe_bg_pi_deg");
  auto hist_fe_pi_ang_cos = (TH1D*)fefile->Get("hist_fe_bg_pi_cos");
  auto hist_fe_p_ang_deg = (TH1D*)fefile->Get("hist_fe_bg_p_deg");
  auto hist_fe_p_ang_cos = (TH1D*)fefile->Get("hist_fe_bg_p_cos");
  double fescale = 1. / 999. / 33.156 * 0.47 * 0.99;
  fescale *= (70. * 0.5 * 8.03) / (58. * 2.3);
  fescale *= 1.2;
  hist_fe_mu_ang_deg->Scale(fescale);
  hist_fe_mu_ang_cos->Scale(fescale);
  hist_fe_pi_ang_deg->Scale(fescale);
  hist_fe_pi_ang_cos->Scale(fescale);
  hist_fe_p_ang_deg->Scale(fescale);
  hist_fe_p_ang_cos->Scale(fescale);

  // Anti neutrino interaction background
  TString anufilename = "/hsm/nu/ninja/pra_tmp/CC0pi_20221213/anu/output/bg_dist.root";
  TFile *anufile = new TFile(anufilename, "read");
  auto hist_anu_mu_ang_deg = (TH1D*)anufile->Get("hist_anu_bg_mu_deg");
  auto hist_anu_mu_ang_cos = (TH1D*)anufile->Get("hist_anu_bg_mu_cos");
  auto hist_anu_pi_ang_deg = (TH1D*)anufile->Get("hist_anu_bg_pi_deg");
  auto hist_anu_pi_ang_cos = (TH1D*)anufile->Get("hist_anu_bg_pi_cos");
  auto hist_anu_p_ang_deg = (TH1D*)anufile->Get("hist_anu_bg_p_deg");
  auto hist_anu_p_ang_cos = (TH1D*)anufile->Get("hist_anu_bg_p_cos");
  double anuscale = 1. / 977./ 33.156 * 0.47 * 0.99;
  hist_anu_mu_ang_deg->Scale(anuscale);
  hist_anu_mu_ang_cos->Scale(anuscale);
  hist_anu_pi_ang_deg->Scale(anuscale);
  hist_anu_pi_ang_cos->Scale(anuscale);
  hist_anu_p_ang_deg->Scale(anuscale);
  hist_anu_p_ang_cos->Scale(anuscale);

  TH1D *hist_bg_int_mu_deg = new TH1D("hist_bg_int_mu_deg", "", muon_deg_bin_size-1, muon_deg_bins);
  TH1D *hist_bg_int_mu_cos = new TH1D("hist_bg_int_mu_cos", "", muon_cos_bin_size-1, muon_cos_bins);
  TH1D *hist_bg_int_pi_deg = new TH1D("hist_bg_int_pi_deg", "", hadron_deg_bin_size-1, hadron_deg_bins);
  TH1D *hist_bg_int_pi_cos = new TH1D("hist_bg_int_pi_cos", "", hadron_cos_bin_size-1, hadron_cos_bins);
  TH1D *hist_bg_int_p_deg = new TH1D("hist_bg_int_p_deg", "", hadron_deg_bin_size-1, hadron_deg_bins);
  TH1D *hist_bg_int_p_cos = new TH1D("hist_bg_int_p_cos", "", hadron_cos_bin_size-1, hadron_cos_bins);
  TList *int_list_mu_deg = new TList();
  TList *int_list_mu_cos = new TList();
  TList *int_list_pi_deg = new TList();
  TList *int_list_pi_cos = new TList();
  TList *int_list_p_deg = new TList();
  TList *int_list_p_cos = new TList();

  int_list_mu_deg->Add(hist_fe_mu_ang_deg);
  int_list_mu_deg->Add(hist_anu_mu_ang_deg);
  int_list_mu_deg->Add(hist_muon_deg_pack_bg);
  int_list_mu_deg->Add(hist_muon_misid_ang_deg);
  hist_bg_int_mu_deg->Merge(int_list_mu_deg);
  hist_bg_int_mu_deg->SetFillStyle(3021);
  hist_bg_int_mu_deg->SetFillColor(kMagenta);

  int_list_mu_cos->Add(hist_fe_mu_ang_cos);
  int_list_mu_cos->Add(hist_anu_mu_ang_cos);
  int_list_mu_cos->Add(hist_muon_cos_pack_bg);
  int_list_mu_cos->Add(hist_muon_misid_ang_cos);
  hist_bg_int_mu_cos->Merge(int_list_mu_cos);
  hist_bg_int_mu_cos->SetFillStyle(3021);
  hist_bg_int_mu_cos->SetFillColor(kMagenta);

  int_list_pi_deg->Add(hist_fe_pi_ang_deg);
  int_list_pi_deg->Add(hist_anu_pi_ang_deg);
  int_list_pi_deg->Add(hist_pion_ang_deg_muon_misid);
  hist_bg_int_pi_deg->Merge(int_list_pi_deg);
  hist_bg_int_pi_deg->SetFillStyle(3021);
  hist_bg_int_pi_deg->SetFillColor(kMagenta);

  int_list_pi_cos->Add(hist_fe_pi_ang_cos);
  int_list_pi_cos->Add(hist_anu_pi_ang_cos);
  int_list_pi_cos->Add(hist_pion_ang_cos_muon_misid);
  hist_bg_int_pi_cos->Merge(int_list_pi_cos);
  hist_bg_int_pi_cos->SetFillStyle(3021);
  hist_bg_int_pi_cos->SetFillColor(kMagenta);

  int_list_p_deg->Add(hist_fe_p_ang_deg);
  int_list_p_deg->Add(hist_anu_p_ang_deg);
  int_list_p_deg->Add(hist_proton_ang_deg_muon_misid);
  hist_bg_int_p_deg->Merge(int_list_p_deg);
  hist_bg_int_p_deg->SetFillStyle(3021);
  hist_bg_int_p_deg->SetFillColor(kMagenta);

  int_list_p_cos->Add(hist_fe_p_ang_cos);
  int_list_p_cos->Add(hist_anu_p_ang_cos);
  int_list_p_cos->Add(hist_proton_ang_cos_muon_misid);
  hist_bg_int_p_cos->Merge(int_list_p_cos);
  hist_bg_int_p_cos->SetFillStyle(3021);
  hist_bg_int_p_cos->SetFillColor(kMagenta);

  l_muon_ang->AddEntry(hist_bg_int_mu_deg, "Internal Beam-related Background", "f");
  // l_muon_ang->AddEntry(hist_bg_ext_mu_deg, "External Beam-related Background", "f");
  l_proton_ang->AddEntry(hist_bg_int_p_deg, "Internal Beam-related Background", "f");
  // l_proton_ang->AddEntry(hist_bg_ext_p_deg, "External Beam-related Background", "f");

  l_proton_ang->AddEntry(hist_proton_misid_ang_deg, "Partner mis-id", "f");

  // Chance coincidence
  TH1D *hist_bg_cc_deg = new TH1D("hist_bg_cc_deg", "", muon_deg_bin_size-1, muon_deg_bins);
  TH1D *hist_bg_cc_cos = new TH1D("hist_bg_cc_cos", "", muon_cos_bin_size-1, muon_cos_bins);

  TH1D *hist_bg_cc_deg_p = new TH1D("hist_bg_cc_deg_p", "", hadron_deg_bin_size-1, hadron_deg_bins);
  TH1D *hist_bg_cc_cos_p = new TH1D("hist_bg_cc_cos_p", "", hadron_cos_bin_size-1, hadron_cos_bins);

  // l_muon_ang->AddEntry(hist_bg_cc_deg, "Beam-unrelated Background", "f");
  // l_proton_ang->AddEntry(hist_bg_cc_deg_p, "Beam-unrelated Background", "f");

  // Data
  TH1D *hist_deg = new TH1D("hist_deg", "", muon_deg_bin_size-1, muon_deg_bins);
  TH1D *hist_cos = new TH1D("hist_cos", "", muon_cos_bin_size-1, muon_cos_bins);
  TH1D *hist_deg_p = new TH1D("hist_deg_p", "", hadron_deg_bin_size-1, hadron_deg_bins);
  TH1D *hist_cos_p = new TH1D("hist_cos_p", "", hadron_cos_bin_size-1, hadron_cos_bins);

  l_muon_ang->AddEntry(hist_deg, "Data", "lp");
  l_proton_ang->AddEntry(hist_deg_p, "Data", "lp");

  // Stack
  THStack *hs_muon_deg = new THStack("hs_muon_deg", "Muon angle;#theta_{#mu} [deg];Entries");
  THStack *hs_muon_cos = new THStack("hs_muon_cos", "Muon angle;#theta_{#mu} [deg];Entries");
  THStack *hs_pion_deg = new THStack("hs_pion_deg", "Pion angle;#theta_{#pi} [deg];Entries");
  THStack *hs_pion_cos = new THStack("hs_pion_cos", "Pion angle;#theta_{#pi} [deg];Entries");
  THStack *hs_proton_deg = new THStack("hs_proton_deg", "Proton angle;#theta_{p} [deg];Entries");
  THStack *hs_proton_cos = new THStack("hs_proton_cos", "Proton angle;#theta_{p} [deg];Entries");

  hs_pion_deg->Add(hist_pion_misid_ang_deg);
  hs_pion_cos->Add(hist_pion_misid_ang_cos);
  hs_proton_deg->Add(hist_proton_misid_ang_deg);
  hs_proton_cos->Add(hist_proton_misid_ang_cos);

  // hs_muon_deg->Add(hist_bg_ext_mu_deg);
  hs_muon_deg->Add(hist_bg_int_mu_deg);
  // hs_muon_cos->Add(hist_bg_ext_mu_cos);
  hs_muon_cos->Add(hist_bg_int_mu_cos);

  // hs_pion_deg->Add(hist_bg_ext_pi_deg);
  hs_pion_deg->Add(hist_bg_int_pi_deg);
  // hs_pion_cos->Add(hist_bg_ext_pi_cos);
  hs_pion_cos->Add(hist_bg_int_pi_cos);

  // hs_proton_deg->Add(hist_bg_ext_p_deg);
  hs_proton_deg->Add(hist_bg_int_p_deg);
  // hs_proton_cos->Add(hist_bg_ext_p_cos);
  hs_proton_cos->Add(hist_bg_int_p_cos);

  for ( int i = 1; i < num_ninja_mode; i++ ) {

    hs_muon_deg->Add(hist_muon_deg[mode_stack_order[i]]);
    hs_pion_deg->Add(hist_pion_deg[mode_stack_order[i]]);
    hs_proton_deg->Add(hist_proton_deg[mode_stack_order[i]]);
    hs_muon_cos->Add(hist_muon_cos[mode_stack_order[i]]);
    hs_pion_cos->Add(hist_pion_cos[mode_stack_order[i]]);
    hs_proton_cos->Add(hist_proton_cos[mode_stack_order[i]]);
  }

  TH1D *hist_muon_deg_norm = new TH1D(*((TH1D*)(hs_muon_deg->GetStack()->Last())));
  hist_muon_deg_norm->SetName("hist_muon_deg_norm");

  TFile *ofile = new TFile("~/stack_ang_0pi.root", "recreate");

  ofile->cd();
  hs_muon_deg->Write();
  hs_pion_deg->Write();
  hs_proton_deg->Write();
  hs_muon_cos->Write();
  hs_pion_cos->Write();
  hs_proton_cos->Write();

  hist_muon_deg_norm->Write();

  l_muon_ang->Write();
  l_proton_ang->Write();
  ofile->Close();

}
