#include "/home/t2k/odagawa/NinjaMCStudy/lib/HistogramStyle.hpp"

void Stack0pi1p() {

  // Legend
  TLegend *l_dat = new TLegend(0.11, 0.45, 0.55, 0.89);
  l_dat->SetName("l_dat");

  // Signal + mis-pid background
  TString filename = "/group/nu/ninja/work/odagawa/20221020-phd-thesis-preliminary/signal/output/output_mode6.root";
  TFile *file = new TFile(filename, "read");
  double scale = 1. / 976. / 33.156 * 0.47 * 0.99;

  TH1D *hist_dat[num_ninja_mode];
  for ( int i = 0; i < num_ninja_mode - 1; i++ ) {
    hist_dat[i] = (TH1D*)file->Get(Form("hist_dalphat_%d", i));
    hist_dat[i]->Scale(scale);
    l_dat->AddEntry(hist_dat[i], mode_name[i], "f");
  }

  // Muon mis-id CC+NC
  auto hist_dat_muon_misid = (TH1D*)file->Get("hist_muon_misid_dalphat");
  hist_dat_muon_misid->Scale(scale);

  // Partner mis-id/non-id
  auto hist_dat_proton_misid = (TH1D*)file->Get("hist_proton_misid_dalphat");
  hist_dat_proton_misid->Scale(scale);
  hist_dat_proton_misid->SetFillColor(kOrange);

  // Iron (+emulsion?) interaction background
  TString fefilename = "/hsm/nu/ninja/pra_tmp/CC0pi_20221213/fe/output/bg_dist.root";
  TFile *fefile = new TFile(fefilename, "read");
  auto fehist_dat = (TH1D*)fefile->Get("hist_dalphat");
  double fescale = 1. / 999. / 33.156 * 0.47 * 0.99;
  fescale *= (70. * 0.05 * 8.03) / (58. * 2.3);
  fescale *= 1.2;
  fehist_dat->Scale(fescale);

  // Anti neutrino interaction background
  TString anufilename = "/hsm/nu/ninja/pra_tmp/CC0pi_20221213/anu/output/bg_dist.root";
  TFile *anufile = new TFile(anufilename, "read");
  auto anuhist_dat = (TH1D*)anufile->Get("hist_dalphat");
  double anuscale = 1. / 977. / 33.156 * 0.47 * 0.99;
  anuhist_dat->Scale(anuscale);

  TH1D *hist_bg_int_dat = new TH1D("hist_bg_int_dat", "", dalphat_bin_size-1, dalphat_bins);
  TList *int_dat_list = new TList();
  int_dat_list->Add(fehist_dat);
  int_dat_list->Add(anuhist_dat);
  int_dat_list->Add(hist_dat_muon_misid);
  hist_bg_int_dat->Merge(int_dat_list);
  hist_bg_int_dat->SetFillStyle(3021);
  hist_bg_int_dat->SetFillColor(kMagenta);

  l_dat->AddEntry(hist_bg_int_dat, "Internal Beam-related Background", "f");

  l_dat->AddEntry(hist_dat_proton_misid, "Partner mis-id", "f");

  TH1D *hist = new TH1D("hist", "", 1, 0, 1);

  l_dat->AddEntry(hist, "Data", "lp");

  THStack *hs_dat = new THStack("hs_dat", ";#delta#alpha_{T} [deg];Entries/(20 deg)");

  hs_dat->Add(hist_dat_proton_misid);
  hs_dat->Add(hist_bg_int_dat);
  for( int i = 1; i < num_ninja_mode; i++ ) {
    hs_dat->Add(hist_dat[mode_stack_order[i]]);
  }

  TFile *ofile = new TFile("~/stack_0pi1p.root", "recreate");

  ofile->cd();
  hs_dat->Write();
  l_dat->Write();
  ofile->Close();

}
