#include "/home/t2k/odagawa/NinjaMCStudy/lib/HistogramStyle.hpp"

void StackModeAng() {

  TString filename = "/hsm/nu/ninja/pra_tmp/mc_tmp_20220505/output/output_mode3.root";
  TFile *file = new TFile(filename, "read");

  THStack *hs_muon_ang = new THStack("hs_muon_ang", "Muon angle;#theta_{#mu} [deg];Entries");
  THStack *hs_pion_ang = new THStack("hs_pion_ang", "Pion angle;#theta_{#pi} [deg];Entries");
  THStack *hs_proton_ang = new THStack("hs_proton_ang", "Proton angle;#theta_{p} [deg];Entries");

  TLegend *l_muon_ang = new TLegend(0.6, 0.5, 0.85, 0.85);
  l_muon_ang->SetName("l_muon_ang");

  TH1D *hist_muon_ang[num_ninja_mode];
  TH1D *hist_pion_ang[num_ninja_mode];
  TH1D *hist_proton_ang[num_ninja_mode];
  for ( int i = 0; i < num_ninja_mode; i++ ) {
    hist_muon_ang[i] = (TH1D*)file->Get(Form("hist_muon_ang_deg_%d", i));
    hist_pion_ang[i] = (TH1D*)file->Get(Form("hist_pion_ang_deg_%d", i));
    hist_proton_ang[i] = (TH1D*)file->Get(Form("hist_proton_ang_deg_%d", i));
    l_muon_ang->AddEntry(hist_muon_ang[i], mode_name[i], "f");
  }

  for ( int i = 0; i < num_ninja_mode; i++ ) {
    hs_muon_ang->Add(hist_muon_ang[mode_stack_order[i]]);
    hs_pion_ang->Add(hist_pion_ang[mode_stack_order[i]]);
    hs_proton_ang->Add(hist_proton_ang[mode_stack_order[i]]);
  }

  TFile *ofile = new TFile("~/stack_ang.root", "recreate");

  ofile->cd();
  hs_muon_ang->Write();
  hs_pion_ang->Write();
  hs_proton_ang->Write();
  l_muon_ang->Write();
  ofile->Close();

}
