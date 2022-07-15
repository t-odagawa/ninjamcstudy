#include "/home/t2k/odagawa/NinjaMCStudy/lib/HistogramStyle.hpp"

void StackModeMom() {
  
  TString filename = "/hsm/nu/ninja/pra_tmp/mc_tmp_20220620/output/output_mode2.root";
  TFile *file = new TFile(filename, "read");

  THStack *hs_muon_mom = new THStack("hs_muon_mom", "Muon momentum;p_{#mu} [MeV/c];Entries");
  THStack *hs_pion_mom = new THStack("hs_pion_mom", "Pion momentum;p_{#pi} [MeV/c];Entries");
  THStack *hs_proton_mom = new THStack("hs_proton_mom", "Proton momentum;p_{p} [MeV/c];Entries");

  THStack *hs_muon_mom_method = new THStack("hs_muon_mom_method", "Muon momentum;p_{#mu} [MeV/c];Entries");
  THStack *hs_pion_mom_method = new THStack("hs_pion_mom_method", "Pion momentum;p_{#pi} [MeV/c];Entries");
  THStack *hs_proton_mom_method = new THStack("hs_proton_mom_method", "Proton momentum;p_{p} [MeV/c];Entries");

  TLegend *l_muon_mom = new TLegend(0.6, 0.5, 0.85, 0.85);
  l_muon_mom->SetName("l_muon_mom");

  TLegend *l_method = new TLegend(0.6, 0.5, 0.85, 0.85);
  l_method->SetName("l_method");

  double scale = 1. / 976. / 8.98 * 0.47 * 0.99;

  TH1D *hist_muon_mom[num_ninja_mode];
  TH1D *hist_pion_mom[num_ninja_mode];
  TH1D *hist_proton_mom[num_ninja_mode];
  for ( int i = 0; i < num_ninja_mode; i++ ) {
    hist_muon_mom[i] = (TH1D*)file->Get(Form("hist_muon_mom_%d", i));
    hist_muon_mom[i]->Scale(scale);
    hist_pion_mom[i] = (TH1D*)file->Get(Form("hist_pion_mom_%d", i));
    hist_pion_mom[i]->Scale(scale);
    hist_proton_mom[i] = (TH1D*)file->Get(Form("hist_proton_mom_%d", i));
    hist_proton_mom[i]->Scale(scale);
    l_muon_mom->AddEntry(hist_muon_mom[i], mode_name[i], "f");    
  }

  for ( int i = 0; i < num_ninja_mode; i++ ) {
    hs_muon_mom->Add(hist_muon_mom[mode_stack_order[i]]);
    hs_pion_mom->Add(hist_pion_mom[mode_stack_order[i]]);
    hs_proton_mom->Add(hist_proton_mom[mode_stack_order[i]]);
  }

  TH1D *hist_muon_mcs = (TH1D*)file->Get("hist_muon_mom_mcs");
  TH1D *hist_muon_range = (TH1D*)file->Get("hist_muon_mom_range");
  TH1D *hist_pion_mcs = (TH1D*)file->Get("hist_pion_mom_mcs");
  TH1D *hist_pion_range = (TH1D*)file->Get("hist_pion_mom_range");
  TH1D *hist_proton_mcs = (TH1D*)file->Get("hist_proton_mom_mcs");
  TH1D *hist_proton_range = (TH1D*)file->Get("hist_proton_mom_range");

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

  hs_muon_mom_method->Add(hist_muon_range);
  hs_muon_mom_method->Add(hist_muon_mcs);
  hs_pion_mom_method->Add(hist_pion_range);
  hs_pion_mom_method->Add(hist_pion_mcs);
  hs_proton_mom_method->Add(hist_proton_range);
  hs_proton_mom_method->Add(hist_proton_mcs);
  l_method->AddEntry(hist_pion_mcs, "MCS", "f");
  l_method->AddEntry(hist_pion_range, "Range", "f");

  TFile *ofile = new TFile("~/stack_mom.root", "recreate");

  ofile->cd();
  hs_muon_mom->Write();
  hs_pion_mom->Write();
  hs_proton_mom->Write();
  hs_muon_mom_method->Write();
  hs_pion_mom_method->Write();
  hs_proton_mom_method->Write();
  l_muon_mom->Write();
  l_method->Write();
  ofile->Close();
  
}
