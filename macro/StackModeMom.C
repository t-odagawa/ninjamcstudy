#include "/home/t2k/odagawa/NinjaMCStudy/lib/HistogramStyle.hpp"

void StackModeMom() {
  
  TString filename = "/hsm/nu/ninja/pra_tmp/mc_tmp_20220505/output/output_mode2.root";
  TFile *file = new TFile(filename, "read");

  THStack *hs_muon_mom = new THStack("hs_muon_mom", "Muon momentum;p_{#mu} [MeV/c];Entries");

  TLegend *l_muon_mom = new TLegend(0.6, 0.5, 0.85, 0.85);
  l_muon_mom->SetName("l_muon_mom");

  TH1D *hist_muon_mom[num_ninja_mode];
  for ( int i = 0; i < num_ninja_mode; i++ ) {
    hist_muon_mom[i] = (TH1D*)file->Get(Form("hist_muon_mom_%d", i));
    l_muon_mom->AddEntry(hist_muon_mom[i], mode_name[i], "f");
  }

  for ( int i = 0; i < num_ninja_mode; i++ ) {
    hs_muon_mom->Add(hist_muon_mom[mode_stack_order[i]]);
  }

  TFile *ofile = new TFile("~/stack_mom.root", "recreate");

  ofile->cd();
  hs_muon_mom->Write();
  l_muon_mom->Write();
  ofile->Close();
  
}
