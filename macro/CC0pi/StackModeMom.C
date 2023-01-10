#include "/home/t2k/odagawa/NinjaMCStudy/lib/HistogramStyle.hpp"

void StackModeMom() {

  // Legend 
  TLegend *l_muon_mom = new TLegend(0.6, 0.5, 0.85, 0.85);
  l_muon_mom->SetName("l_muon_mom");
  TLegend *l_proton_mom = new TLegend(0.6, 0.5, 0.85, 0.85);
  l_proton_mom->SetName("l_proton_mom");
  TLegend *l_method = new TLegend(0.6, 0.5, 0.85, 0.85);
  l_method->SetName("l_method");

  // Signal + packing/mis-pid background
  TString filename = "/hsm/nu/ninja/pra_tmp/CC0pi_20221213/signal/output/output_mode2.root";
  TFile *file = new TFile(filename, "read");
  double scale = 1. / 976. / 33.156 * 0.47 * 0.99;

  // Signal method only
  auto hist_muon_mcs = (TH1D*)file->Get("hist_muon_mom_mcs");
  auto hist_muon_range = (TH1D*)file->Get("hist_muon_mom_range");
  auto hist_pion_mcs = (TH1D*)file->Get("hist_pion_mom_mcs");
  auto hist_pion_range = (TH1D*)file->Get("hist_pion_mom_range");
  auto hist_proton_mcs = (TH1D*)file->Get("hist_proton_mom_mcs");
  auto hist_proton_range = (TH1D*)file->Get("hist_proton_mom_range");

  hist_muon_mcs->SetFillColor(kRed);
  hist_muon_range->SetFillColor(kBlue);
  hist_pion_mcs->SetFillColor(kRed);
  hist_pion_range->SetFillColor(kBlue);
  hist_proton_mcs->SetFillColor(kRed);
  hist_proton_range->SetFillColor(kBlue);

  hist_muon_mcs->Scale(scale);
  hist_muon_range->Scale(scale);
  hist_pion_mcs->Scale(scale);
  hist_pion_range->Scale(scale);
  hist_proton_mcs->Scale(scale);
  hist_proton_range->Scale(scale);

  // Signal  
  TH1D *hist_muon_mom[num_ninja_mode];
  TH1D *hist_muon_mom_mcs[num_ninja_mode];
  TH1D *hist_muon_mom_range[num_ninja_mode];
  TH1D *hist_pion_mom[num_ninja_mode];
  TH1D *hist_pion_mom_mcs[num_ninja_mode];
  TH1D *hist_pion_mom_range[num_ninja_mode];
  TH1D *hist_proton_mom[num_ninja_mode];
  TH1D *hist_proton_mom_mcs[num_ninja_mode];
  TH1D *hist_proton_mom_range[num_ninja_mode];

  // TH1D *hist_muon_mom_mcs_all[num_ninja_mode];

  for ( int i = 0; i < num_ninja_mode - 1; i++ ) {
    hist_muon_mom[i] = (TH1D*)file->Get(Form("hist_muon_mom_%d", i));
    hist_muon_mom_mcs[i] = (TH1D*)file->Get(Form("hist_muon_mom_mcs_%d", i));
    hist_muon_mom_range[i] = (TH1D*)file->Get(Form("hist_muon_mom_range_%d", i));
    hist_muon_mom[i]->Scale(scale);
    hist_muon_mom_mcs[i]->Scale(scale);
    hist_muon_mom_range[i]->Scale(scale);
    hist_pion_mom[i] = (TH1D*)file->Get(Form("hist_pion_mom_%d", i));
    hist_pion_mom_mcs[i] = (TH1D*)file->Get(Form("hist_pion_mom_mcs_%d", i));
    hist_pion_mom_range[i] = (TH1D*)file->Get(Form("hist_pion_mom_range_%d", i));
    hist_pion_mom[i]->Scale(scale);
    hist_pion_mom_mcs[i]->Scale(scale);
    hist_pion_mom_range[i]->Scale(scale);
    hist_proton_mom[i] = (TH1D*)file->Get(Form("hist_proton_mom_%d", i));
    hist_proton_mom_mcs[i] = (TH1D*)file->Get(Form("hist_proton_mom_mcs_%d", i));
    hist_proton_mom_range[i] = (TH1D*)file->Get(Form("hist_proton_mom_range_%d", i));
    hist_proton_mom[i]->Scale(scale);
    hist_proton_mom_mcs[i]->Scale(scale);
    hist_proton_mom_range[i]->Scale(scale);

    // hist_muon_mom_mcs_all[i] = (TH1D*)file->Get(Form("hist_muon_mom_mcs_all_%d", i));
    // hist_muon_mom_mcs_all[i]->Scale(scale);

    l_muon_mom->AddEntry(hist_muon_mom[i], mode_name[i], "f");    
    l_proton_mom->AddEntry(hist_proton_mom[i], mode_name[i], "f");
  }

  // Packing background
  auto hist_muon_pack_bg = (TH1D*)file->Get("hist_muon_mom_single");
  auto hist_muon_mcs_pack_bg = (TH1D*)file->Get("hist_muon_mom_single_mcs");
  auto hist_muon_range_pack_bg = (TH1D*)file->Get("hist_muon_mom_single_range");
  hist_muon_pack_bg->Scale(scale * 0.09);
  hist_muon_mcs_pack_bg->Scale(scale * 0.09);
  hist_muon_range_pack_bg->Scale(scale * 0.09);

  // auto hist_muon_pack_bg_mcs_all = (TH1D*)file->Get("hist_muon_mom_single_mcs_all");
  // hist_muon_pack_bg_mcs_all->Scale(scale * 0.09);

  // Muon mis-id CC + NC
  auto hist_muon_misid_mom = (TH1D*)file->Get("hist_muon_misid_mom");
  auto hist_muon_misid_mom_mcs = (TH1D*)file->Get("hist_muon_misid_mom_mcs");
  auto hist_muon_misid_mom_range = (TH1D*)file->Get("hist_muon_misid_mom_range");
  auto hist_pion_mom_muon_misid = (TH1D*)file->Get("hist_pion_mom_muon_misid");
  auto hist_pion_mom_mcs_muon_misid = (TH1D*)file->Get("hist_pion_mom_mcs_muon_misid");
  auto hist_pion_mom_range_muon_misid = (TH1D*)file->Get("hist_pion_mom_range_muon_misid");
  auto hist_proton_mom_muon_misid = (TH1D*)file->Get("hist_proton_mom_muon_misid");
  auto hist_proton_mom_mcs_muon_misid = (TH1D*)file->Get("hist_proton_mom_mcs_muon_misid");
  auto hist_proton_mom_range_muon_misid = (TH1D*)file->Get("hist_proton_mom_range_muon_misid");
  hist_muon_misid_mom->Scale(scale);
  hist_muon_misid_mom_mcs->Scale(scale);
  hist_muon_misid_mom_range->Scale(scale);
  hist_pion_mom_muon_misid->Scale(scale);
  hist_pion_mom_mcs_muon_misid->Scale(scale);
  hist_pion_mom_range_muon_misid->Scale(scale);
  hist_proton_mom_muon_misid->Scale(scale);
  hist_proton_mom_mcs_muon_misid->Scale(scale);
  hist_proton_mom_range_muon_misid->Scale(scale);

  // auto hist_muon_misid_mom_mcs_all = (TH1D*)file->Get("hist_muon_misid_mom_mcs_all");
  // hist_muon_misid_mom_mcs_all->Scale(scale);

  // Partner mis-id
  auto hist_pion_misid_mom = (TH1D*)file->Get("hist_pion_misid_mom");
  auto hist_pion_misid_mom_mcs = (TH1D*)file->Get("hist_pion_misid_mom_mcs");
  auto hist_pion_misid_mom_range = (TH1D*)file->Get("hist_pion_misid_mom_range");
  auto hist_proton_misid_mom = (TH1D*)file->Get("hist_proton_misid_mom");
  auto hist_proton_misid_mom_mcs = (TH1D*)file->Get("hist_proton_misid_mom_mcs");
  auto hist_proton_misid_mom_range = (TH1D*)file->Get("hist_proton_misid_mom_range");
  hist_pion_misid_mom->Scale(scale);
  hist_pion_misid_mom->SetFillColor(kOrange);
  hist_pion_misid_mom_mcs->Scale(scale);
  hist_pion_misid_mom_mcs->SetFillColor(kOrange);
  hist_pion_misid_mom_range->Scale(scale);
  hist_pion_misid_mom_range->SetFillColor(kOrange);
  hist_proton_misid_mom->Scale(scale);
  hist_proton_misid_mom->SetFillColor(kOrange);
  hist_proton_misid_mom_mcs->Scale(scale);
  hist_proton_misid_mom_mcs->SetFillColor(kOrange);
  hist_proton_misid_mom_range->Scale(scale);
  hist_proton_misid_mom_range->SetFillColor(kOrange);

  // Iron (+ emulsion?) interaction background
  // TString fefilename = "/hsm/nu/ninja/pra_tmp/mc_tmp_fe_20220712/output/bg_dist.root";
  TString fefilename = "/group/nu/ninja/work/odagawa/20220930-iron-mc-new-matching/output/bg_dist.root";
  TFile *fefile = new TFile(fefilename, "read");
  auto hist_fe_mu_mom = (TH1D*)fefile->Get("hist_fe_bg_mu_mom");
  auto hist_fe_mu_mom_mcs = (TH1D*)fefile->Get("hist_fe_bg_mu_mom_mcs");
  auto hist_fe_mu_mom_range = (TH1D*)fefile->Get("hist_fe_bg_mu_mom_range");
  auto hist_fe_p_mom = (TH1D*)fefile->Get("hist_fe_bg_p_mom");
  auto hist_fe_p_mom_mcs = (TH1D*)fefile->Get("hist_fe_bg_p_mom_mcs");
  auto hist_fe_p_mom_range = (TH1D*)fefile->Get("hist_fe_bg_p_mom_range");
  auto hist_fe_pi_mom = (TH1D*)fefile->Get("hist_fe_bg_pi_mom");
  auto hist_fe_pi_mom_mcs = (TH1D*)fefile->Get("hist_fe_bg_pi_mom_mcs");
  auto hist_fe_pi_mom_range = (TH1D*)fefile->Get("hist_fe_bg_pi_mom_range");
  double fescale = 1. / 999. / 33.156 * 0.47 * 0.99;
  fescale *= (70. * 0.5 * 8.03) / (58. * 2.3);
  fescale *= 1.2;
  hist_fe_mu_mom->Scale(fescale);
  hist_fe_mu_mom_mcs->Scale(fescale);
  hist_fe_mu_mom_range->Scale(fescale);
  hist_fe_pi_mom->Scale(fescale);
  hist_fe_pi_mom_mcs->Scale(fescale);
  hist_fe_pi_mom_range->Scale(fescale);
  hist_fe_p_mom->Scale(fescale);
  hist_fe_p_mom_mcs->Scale(fescale);
  hist_fe_p_mom_range->Scale(fescale);

  // auto hist_fe_mu_mom_mcs_all = (TH1D*)fefile->Get("hist_fe_bg_mu_mom_mcs_all");
  // hist_fe_mu_mom_mcs_all->Scale(fescale);

  // Anti neutrino interaction background
  // TString anufilename = "/hsm/nu/ninja/pra_tmp/mc_tmp_anu_20220714/output/bg_dist.root";
  TString anufilename = "/group/nu/ninja/work/odagawa/20220930-bkg-mc-new-matching/anu/output/bg_dist.root";
  TFile *anufile = new TFile(anufilename, "read");
  auto hist_anu_mu_mom = (TH1D*)anufile->Get("hist_anu_bg_mu_mom");
  auto hist_anu_mu_mom_mcs = (TH1D*)anufile->Get("hist_anu_bg_mu_mom_mcs");
  auto hist_anu_mu_mom_range = (TH1D*)anufile->Get("hist_anu_bg_mu_mom_range");
  auto hist_anu_p_mom = (TH1D*)anufile->Get("hist_anu_bg_p_mom");
  auto hist_anu_p_mom_mcs = (TH1D*)anufile->Get("hist_anu_bg_p_mom_mcs");
  auto hist_anu_p_mom_range = (TH1D*)anufile->Get("hist_anu_bg_p_mom_range");
  auto hist_anu_pi_mom = (TH1D*)anufile->Get("hist_anu_bg_pi_mom");
  auto hist_anu_pi_mom_mcs = (TH1D*)anufile->Get("hist_anu_bg_pi_mom_mcs");
  auto hist_anu_pi_mom_range = (TH1D*)anufile->Get("hist_anu_bg_pi_mom_range");
  double anuscale = 1. / 977. / 33.156 * 0.47 * 0.99;
  hist_anu_mu_mom->Scale(anuscale);
  hist_anu_mu_mom_mcs->Scale(anuscale);
  hist_anu_mu_mom_range->Scale(anuscale);
  hist_anu_pi_mom->Scale(anuscale);
  hist_anu_pi_mom_mcs->Scale(anuscale);
  hist_anu_pi_mom_range->Scale(anuscale);
  hist_anu_p_mom->Scale(anuscale);
  hist_anu_p_mom_mcs->Scale(anuscale);
  hist_anu_p_mom_range->Scale(anuscale);

  // auto hist_anu_mu_mom_mcs_all = (TH1D*)anufile->Get("hist_anu_bg_mu_mom_mcs_all");
  // hist_anu_mu_mom_mcs_all->Scale(anuscale);

  TH1D *hist_bg_int_mu_mom = new TH1D("hist_bg_int_mu_mom", "", muon_mom_bin_size-1, muon_mom_bins);
  TH1D *hist_bg_int_mu_mom_mcs = new TH1D("hist_bg_int_mu_mom_mcs", "", muon_mom_bin_size-1, muon_mom_bins);
  TH1D *hist_bg_int_mu_mom_range = new TH1D("hist_bg_int_mu_mom_range", "", muon_mom_bin_size-1, muon_mom_bins);
  TH1D *hist_bg_int_pi_mom = new TH1D("hist_bg_int_pi_mom", "", hadron_mom_bin_size-1, hadron_mom_bins);
  TH1D *hist_bg_int_pi_mom_mcs = new TH1D("hist_bg_int_pi_mom_mcs", "", hadron_mom_bin_size-1, hadron_mom_bins);
  TH1D *hist_bg_int_pi_mom_range = new TH1D("hist_bg_int_pi_mom_range", "", hadron_mom_bin_size-1, hadron_mom_bins);
  TH1D *hist_bg_int_p_mom = new TH1D("hist_bg_int_p_mom", "", hadron_mom_bin_size-1, hadron_mom_bins);
  TH1D *hist_bg_int_p_mom_mcs = new TH1D("hist_bg_int_p_mom_mcs", "", hadron_mom_bin_size-1, hadron_mom_bins);
  TH1D *hist_bg_int_p_mom_range = new TH1D("hist_bg_int_p_mom_range", "", hadron_mom_bin_size-1, hadron_mom_bins);
  TList *int_list_mu_mom = new TList();
  TList *int_list_mu_mom_mcs = new TList();
  TList *int_list_mu_mom_range = new TList();
  TList *int_list_pi_mom = new TList();
  TList *int_list_pi_mom_mcs = new TList();
  TList *int_list_pi_mom_range = new TList();
  TList *int_list_p_mom = new TList();
  TList *int_list_p_mom_mcs = new TList();
  TList *int_list_p_mom_range = new TList();

  // TH1D *hist_bg_int_mu_mom_mcs_all = new TH1D("hist_bg_int_mu_mom_mcs_all","", 20,0,2000);
  // TList *int_list_mu_mom_mcs_all = new TList();

  int_list_mu_mom->Add(hist_fe_mu_mom);
  int_list_mu_mom->Add(hist_anu_mu_mom);
  int_list_mu_mom->Add(hist_muon_pack_bg);
  int_list_mu_mom->Add(hist_muon_misid_mom);
  hist_bg_int_mu_mom->Merge(int_list_mu_mom);
  hist_bg_int_mu_mom->SetFillStyle(3021);
  hist_bg_int_mu_mom->SetFillColor(kMagenta);

  int_list_mu_mom_mcs->Add(hist_fe_mu_mom_mcs);
  int_list_mu_mom_mcs->Add(hist_anu_mu_mom_mcs);
  int_list_mu_mom_mcs->Add(hist_muon_mcs_pack_bg);
  int_list_mu_mom_mcs->Add(hist_muon_misid_mom_mcs);
  hist_bg_int_mu_mom_mcs->Merge(int_list_mu_mom_mcs);
  hist_bg_int_mu_mom_mcs->SetFillStyle(3021);
  hist_bg_int_mu_mom_mcs->SetFillColor(kMagenta);

  int_list_mu_mom_range->Add(hist_fe_mu_mom_range);
  int_list_mu_mom_range->Add(hist_anu_mu_mom_range);
  int_list_mu_mom_range->Add(hist_muon_range_pack_bg);
  int_list_mu_mom_range->Add(hist_muon_misid_mom_range);
  hist_bg_int_mu_mom_range->Merge(int_list_mu_mom_range);
  hist_bg_int_mu_mom_range->SetFillStyle(3021);
  hist_bg_int_mu_mom_range->SetFillColor(kMagenta);

  int_list_pi_mom->Add(hist_fe_pi_mom);
  int_list_pi_mom->Add(hist_anu_pi_mom);
  int_list_pi_mom->Add(hist_pion_mom_muon_misid);
  hist_bg_int_pi_mom->Merge(int_list_pi_mom);
  hist_bg_int_pi_mom->SetFillStyle(3021);
  hist_bg_int_pi_mom->SetFillColor(kMagenta);

  int_list_pi_mom_mcs->Add(hist_fe_pi_mom_mcs);
  int_list_pi_mom_mcs->Add(hist_anu_pi_mom_mcs);
  int_list_pi_mom_mcs->Add(hist_pion_mom_mcs_muon_misid);
  hist_bg_int_pi_mom_mcs->Merge(int_list_pi_mom_mcs);
  hist_bg_int_pi_mom_mcs->SetFillStyle(3021);
  hist_bg_int_pi_mom_mcs->SetFillColor(kMagenta);

  int_list_pi_mom_range->Add(hist_fe_pi_mom_range);
  int_list_pi_mom_range->Add(hist_anu_pi_mom_range);
  int_list_pi_mom_range->Add(hist_pion_mom_range_muon_misid);
  hist_bg_int_pi_mom_range->Merge(int_list_pi_mom_range);
  hist_bg_int_pi_mom_range->SetFillStyle(3021);
  hist_bg_int_pi_mom_range->SetFillColor(kMagenta);

  int_list_p_mom->Add(hist_fe_p_mom);
  int_list_p_mom->Add(hist_anu_p_mom);
  int_list_p_mom->Add(hist_proton_mom_muon_misid);
  hist_bg_int_p_mom->Merge(int_list_p_mom);
  hist_bg_int_p_mom->SetFillStyle(3021);
  hist_bg_int_p_mom->SetFillColor(kMagenta);

  int_list_p_mom_mcs->Add(hist_fe_p_mom_mcs);
  int_list_p_mom_mcs->Add(hist_anu_p_mom_mcs);
  int_list_p_mom_mcs->Add(hist_proton_mom_mcs_muon_misid);
  hist_bg_int_p_mom_mcs->Merge(int_list_p_mom_mcs);
  hist_bg_int_p_mom_mcs->SetFillStyle(3021);
  hist_bg_int_p_mom_mcs->SetFillColor(kMagenta);

  int_list_p_mom_range->Add(hist_fe_p_mom_range);
  int_list_p_mom_range->Add(hist_anu_p_mom_range);
  int_list_p_mom_range->Add(hist_proton_mom_range_muon_misid);
  hist_bg_int_p_mom_range->Merge(int_list_p_mom_range);
  hist_bg_int_p_mom_range->SetFillStyle(3021);
  hist_bg_int_p_mom_range->SetFillColor(kMagenta);
  /*
  int_list_mu_mom_mcs_all->Add(hist_fe_mu_mom_mcs_all);
  int_list_mu_mom_mcs_all->Add(hist_anu_mu_mom_mcs_all);
  int_list_mu_mom_mcs_all->Add(hist_muon_pack_bg_mcs_all);
  int_list_mu_mom_mcs_all->Add(hist_muon_misid_mom_mcs_all);
  hist_bg_int_mu_mom_mcs_all->Merge(int_list_mu_mom_mcs_all);
  hist_bg_int_mu_mom_mcs_all->SetFillStyle(3021);
  hist_bg_int_mu_mom_mcs_all->SetFillColor(kMagenta);
  */
  l_muon_mom->AddEntry(hist_bg_int_mu_mom, "Internal Beam-related Background", "f");
  // l_muon_mom->AddEntry(hist_bg_ext_mu_mom, "External Beam-related Background", "f");  
  l_proton_mom->AddEntry(hist_bg_int_p_mom, "Internal Beam-related Background", "f");
  // l_proton_mom->AddEntry(hist_bg_ext_p_mom, "External Beam-related Background", "f");

  l_proton_mom->AddEntry(hist_proton_misid_mom, "Partner mis-id", "f");

  // Chance coincidence
  TH1D *hist_bg_cc = new TH1D("hist_bg_cc", "", muon_mom_bin_size-1, muon_mom_bins);
  TH1D *hist_bg_cc_p = new TH1D("hist_bg_cc_p", "", hadron_mom_bin_size-1, hadron_mom_bins);
  // l_muon_mom->AddEntry(hist_bg_cc, "Beam-unrelated Background", "f");
  // l_proton_mom->AddEntry(hist_bg_cc_p, "Beam-unrelated Background", "f");
  
  // Data
  TH1D *hist = new TH1D("hist", "", 20, 0., 2000.);
  TH1D *hist_p = new TH1D("hist_p", "", 15, 0., 1500.);
  l_muon_mom->AddEntry(hist, "Data", "lp");
  l_proton_mom->AddEntry(hist_p, "Data", "lp");
  
  // stack
  THStack *hs_muon_mom = new THStack("hs_muon_mom", "Muon momentum;Momentum [MeV/c];Entries");
  THStack *hs_muon_mom_mcs = new THStack("hs_muon_mom_mcs", "Muon MCS momentum;Momentum [MeV/c];Entries");
  THStack *hs_muon_mom_range = new THStack("hs_muon_mom_range", "Muon range momentum;Momentum [MeV/c];Entries");
  THStack *hs_pion_mom = new THStack("hs_pion_mom", "Pion momentum;Momentum [MeV/c];Entries");
  THStack *hs_pion_mom_mcs = new THStack("hs_pion_mom_mcs", "Pion MCS momentum;Momentum [MeV/c];Entries");
  THStack *hs_pion_mom_range = new THStack("hs_pion_mom_range", "Pion range momentum;Momentum [MeV/c];Entries");
  THStack *hs_proton_mom = new THStack("hs_proton_mom", "Proton momentum;Momentum [MeV/c];Entries");
  THStack *hs_proton_mom_mcs = new THStack("hs_proton_mom_mcs", "Proton MCS momentum;Momentum [MeV/c];Entries");
  THStack *hs_proton_mom_range = new THStack("hs_proton_mom_range", "Proton range momentum;Momentum [MeV/c];Entries");

  THStack *hs_muon_mom_method = new THStack("hs_muon_mom_method", "Muon momentum (signal only);Momentum [MeV/c];Entries");
  THStack *hs_pion_mom_method = new THStack("hs_pion_mom_method", "Pion momentum (signal only);Momentum [MeV/c];Entries");
  THStack *hs_proton_mom_method = new THStack("hs_proton_mom_method", "Proton momentum (signal only);Momentum [MeV/c];Entries");

  // THStack *hs_muon_mom_mcs_all = new THStack("hs_muon_mom_mcs_all", "Muon MCS momentum;Momentum [MeV/c];Entries");

  hs_pion_mom->Add(hist_pion_misid_mom);
  hs_pion_mom_mcs->Add(hist_pion_misid_mom_mcs);
  hs_pion_mom_range->Add(hist_pion_misid_mom_range);
  hs_proton_mom->Add(hist_proton_misid_mom);
  hs_proton_mom_mcs->Add(hist_proton_misid_mom_mcs);
  hs_proton_mom_range->Add(hist_proton_misid_mom_range);

  // hs_muon_mom->Add(hist_bg_ext_mu_mom);
  // hs_muon_mom_mcs->Add(hist_bg_ext_mu_mom_mcs);
  // hs_muon_mom_range->Add(hist_bg_ext_mu_mom_range);
  hs_muon_mom->Add(hist_bg_int_mu_mom);
  hs_muon_mom_mcs->Add(hist_bg_int_mu_mom_mcs);
  hs_muon_mom_range->Add(hist_bg_int_mu_mom_range);

  /*
  hs_muon_mom_mcs_all->Add(hist_bg_ext_mu_mom_mcs_all);
  hs_muon_mom_mcs_all->Add(hist_bg_int_mu_mom_mcs_all);
  */

  // hs_pion_mom->Add(hist_bg_ext_pi_mom);
  // hs_pion_mom_mcs->Add(hist_bg_ext_pi_mom_mcs);
  // hs_pion_mom_range->Add(hist_bg_ext_pi_mom_range);
  hs_pion_mom->Add(hist_bg_int_pi_mom);
  hs_pion_mom_mcs->Add(hist_bg_int_pi_mom_mcs);
  hs_pion_mom_range->Add(hist_bg_int_pi_mom_range);

  // hs_proton_mom->Add(hist_bg_ext_p_mom);
  // hs_proton_mom_mcs->Add(hist_bg_ext_p_mom_mcs);
  // hs_proton_mom_range->Add(hist_bg_ext_p_mom_range);
  hs_proton_mom->Add(hist_bg_int_p_mom);
  hs_proton_mom_mcs->Add(hist_bg_int_p_mom_mcs);
  hs_proton_mom_range->Add(hist_bg_int_p_mom_range);

  for ( int i = 1; i < num_ninja_mode; i++ ) {

    hs_muon_mom->Add(hist_muon_mom[mode_stack_order[i]]);
    hs_muon_mom_mcs->Add(hist_muon_mom_mcs[mode_stack_order[i]]);
    hs_muon_mom_range->Add(hist_muon_mom_range[mode_stack_order[i]]);
    hs_pion_mom->Add(hist_pion_mom[mode_stack_order[i]]);
    hs_pion_mom_mcs->Add(hist_pion_mom_mcs[mode_stack_order[i]]);
    hs_pion_mom_range->Add(hist_pion_mom_range[mode_stack_order[i]]);
    hs_proton_mom->Add(hist_proton_mom[mode_stack_order[i]]);
    hs_proton_mom_mcs->Add(hist_proton_mom_mcs[mode_stack_order[i]]);
    hs_proton_mom_range->Add(hist_proton_mom_range[mode_stack_order[i]]);

    // hs_muon_mom_mcs_all->Add(hist_muon_mom_mcs_all[mode_stack_order[i]]);

  }

  hs_muon_mom_method->Add(hist_muon_range);
  hs_muon_mom_method->Add(hist_muon_mcs);
  hs_pion_mom_method->Add(hist_pion_range);
  hs_pion_mom_method->Add(hist_pion_mcs);
  hs_proton_mom_method->Add(hist_proton_range);
  hs_proton_mom_method->Add(hist_proton_mcs);
  l_method->AddEntry(hist_pion_mcs, "MCS", "f");
  l_method->AddEntry(hist_pion_range, "Range", "f");

  // TH1D *hist_mcs_all_norm = new TH1D(*((TH1D*)(hs_muon_mom_mcs_all->GetStack()->Last())));
  // hist_mcs_all_norm->SetName("hist_mcs_all_norm");
  TH1D *hist_mcs_norm = new TH1D(*((TH1D*)(hs_muon_mom_mcs->GetStack()->Last())));
  hist_mcs_norm->SetName("hist_mcs_norm");
  TH1D *hist_range_norm = new TH1D(*((TH1D*)(hs_muon_mom_range->GetStack()->Last())));
  hist_range_norm->SetName("hist_range_norm");

  TFile *ofile = new TFile("~/stack_mom_0pi.root", "recreate");
  // TFile *ofile = new TFile("~/stack_mom_syst.root", "recreate");

  ofile->cd();
  hs_muon_mom->Write();
  hs_muon_mom_mcs->Write();
  hs_muon_mom_range->Write();
  hs_pion_mom->Write();
  hs_pion_mom_mcs->Write();
  hs_pion_mom_range->Write();
  hs_proton_mom->Write();
  hs_proton_mom_mcs->Write();
  hs_proton_mom_range->Write();
  hs_muon_mom_method->Write();
  hs_pion_mom_method->Write();
  hs_proton_mom_method->Write();

  // hs_muon_mom_mcs_all->Write();
  // hist_mcs_all_norm->Write();
  hist_mcs_norm->Write();
  hist_range_norm->Write();
  l_muon_mom->Write();
  l_proton_mom->Write();
  l_method->Write();
  ofile->Close();
  
}
