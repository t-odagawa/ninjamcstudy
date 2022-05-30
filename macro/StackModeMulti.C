#include "/home/t2k/odagawa/NinjaMCStudy/lib/HistogramStyle.hpp"

void StackModeMulti() {
  
  TString filename = "/hsm/nu/ninja/pra_tmp/mc_tmp_20220505/output/output_mode0.root";
  TFile *file = new TFile(filename, "read");

  THStack *hs_water_multi = new THStack("hs_water_multi","Multiplicity (Area normalized);# of tracks;Entries");
  THStack *hs_water_proton_multi = new THStack("hs_water_proton_multi", "Proton multiplicity;# of protons;Entries");
  THStack *hs_water_pion_multi = new THStack("hs_water_pion_multi", "Pion multiplicity;# of pions;Entries");

  TLegend *l_water_multi = new TLegend(0.6, 0.5, 0.85, 0.85);
  l_water_multi->SetName("l_water_multi");

  TH1D *hist = new TH1D("hist", "", 10, 0.5, 10.5);
  hist->Fill(1,25);
  hist->Fill(2,37);
  hist->Fill(3,22);
  hist->Fill(4,5);
  hist->Fill(5,2);
  hist->Fill(6,3);
  hist->Fill(7,1);
  hist->Sumw2(0);  

  auto hist_water_total_multi = (TH1D*)file->Get("hist_water_total_multi");
  double scale = hist->Integral() / hist_water_total_multi->Integral();

  TH1D *hist_water_multi[num_ninja_mode];
  TH1D *hist_water_proton_multi[num_ninja_mode];
  TH1D *hist_water_pion_multi[num_ninja_mode];
  for ( int i = 0; i < num_ninja_mode; i++ ) {
    hist_water_multi[i] = (TH1D*)file->Get(Form("hist_water_mode_multi_%d",i));
    hist_water_proton_multi[i] = (TH1D*)file->Get(Form("hist_water_mode_proton_multi_%d",i));
    hist_water_pion_multi[i] = (TH1D*)file->Get(Form("hist_water_mode_pion_multi_%d",i));
    l_water_multi->AddEntry(hist_water_multi[i], mode_name[i], "f");
  }

  l_water_multi->AddEntry(hist, "Data", "lp");

  for ( int i = 0; i < num_ninja_mode; i++ ) {
    // hist_water_multi[mode_stack_order[i]]->Scale(0.47e-3);
    hist_water_multi[mode_stack_order[i]]->Scale(scale);
    hs_water_multi->Add(hist_water_multi[mode_stack_order[i]]);
    hs_water_proton_multi->Add(hist_water_proton_multi[mode_stack_order[i]]);
    hs_water_pion_multi->Add(hist_water_pion_multi[mode_stack_order[i]]);
  }


  TFile *ofile = new TFile("~/stack_mode.root", "recreate");

  ofile->cd();
  hs_water_multi->Write();
  hs_water_proton_multi->Write();
  hs_water_pion_multi->Write();
  l_water_multi->Write();
  hist->Write();
  ofile->Close();

}
