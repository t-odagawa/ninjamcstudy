#include "/home/t2k/odagawa/NinjaMCStudy/lib/HistogramStyle.hpp"

void Stack0pi2p() {

  // Legend
  TLegend *l_open_cos = new TLegend(0.11, 0.45, 0.55, 0.89);
  l_open_cos->SetName("l_open_cos");

  // Signal + mis-pid background
  TString filename = "/group/nu/ninja/work/odagawa/20221020-phd-thesis-preliminary/signal/output/output_mode7.root";
  TFile *file = new TFile(filename, "read");
  double scale = 1. / 976. / 33.156 * 0.47 * 0.99;

  TH1D *hist_open_cos[num_ninja_mode];
  for ( int i = 0; i < num_ninja_mode - 1; i++ ) {
    hist_open_cos[i] = (TH1D*)file->Get(Form("hist_open_cos_%d", i));
    hist_open_cos[i]->Scale(scale);
    l_open_cos->AddEntry(hist_open_cos[i], mode_name[i], "f");
  }

  // Muon mis-id CC+NC
  auto hist_open_cos_muon_misid = (TH1D*)file->Get("hist_muon_misid_open_cos");
  hist_open_cos_muon_misid->Scale(scale);

  // Partner mis-id/non-id
  auto hist_open_cos_proton_misid = (TH1D*)file->Get("hist_proton_misid_open_cos");
  hist_open_cos_proton_misid->Scale(scale);
  hist_open_cos_proton_misid->SetFillColor(kOrange);

  // Iron (+emulsion?) interaction background
  TString fefilename = "/hsm/nu/ninja/pra_tmp/CC0pi_20221213/fe/output/bg_dist.root";
  TFile *fefile = new TFile(fefilename, "read");
  auto fehist_open_cos = (TH1D*)fefile->Get("hist_open_cos");
  double fescale = 1. / 999. / 33.156 * 0.47 * 0.99;
  fescale *= (70. * 0.05 * 8.03) / (58. * 2.3);
  fescale *= 1.2;
  fehist_open_cos->Scale(fescale);

  // Anti neutrino interaction background
  TString anufilename = "/hsm/nu/ninja/pra_tmp/CC0pi_20221213/anu/output/bg_dist.root";
  TFile *anufile = new TFile(anufilename, "read");
  auto anuhist_open_cos = (TH1D*)anufile->Get("hist_open_cos");
  double anuscale = 1. / 977. / 33.156 * 0.47 * 0.99;
  anuhist_open_cos->Scale(anuscale);

  TH1D *hist_bg_int_open_cos = new TH1D("hist_bg_int_open_cos", "", open_cos_bin_size-1, open_cos_bins);
  TList *int_open_cos_list = new TList();
  int_open_cos_list->Add(fehist_open_cos);
  int_open_cos_list->Add(anuhist_open_cos);
  int_open_cos_list->Add(hist_open_cos_muon_misid);
  hist_bg_int_open_cos->Merge(int_open_cos_list);
  hist_bg_int_open_cos->SetFillStyle(3021);
  hist_bg_int_open_cos->SetFillColor(kMagenta);

  l_open_cos->AddEntry(hist_bg_int_open_cos, "Internal Beam-related Background", "f");

  l_open_cos->AddEntry(hist_open_cos_proton_misid, "Partner mis-id", "f");

  TH1D *hist = new TH1D("hist", "", 1, 0, 1);

  l_open_cos->AddEntry(hist, "Data", "lp");

  THStack *hs_open_cos = new THStack("hs_open_cos", ";cos#gamma;Entries/(0.05)");

  hs_open_cos->Add(hist_open_cos_proton_misid);
  hs_open_cos->Add(hist_bg_int_open_cos);
  for( int i = 1; i < num_ninja_mode; i++ ) {
    hs_open_cos->Add(hist_open_cos[mode_stack_order[i]]);
  }

  TFile *ofile = new TFile("~/stack_0pi2p.root", "recreate");

  ofile->cd();
  hs_open_cos->Write();
  l_open_cos->Write();
  ofile->Close();

}
