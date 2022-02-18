#include "/home/t2k/odagawa/NinjaMCStudy/lib/HistogramStyle.hpp"

void StackMode() {

  TString filename = "~/NinjaMCStudy/build/testadd.root";
  TFile *file = new TFile(filename, "read");

  // Stacked histograms
  THStack *hs_mom = new THStack("hs_mom", "Momentum;p [GeV/c];Entries/0.1 GeV/c");
  THStack *hs_cos = new THStack("hs_cos", "Angle;cos#theta;Entries/0.01");
  THStack *hs_deg = new THStack("hs_deg", "Angle;#theta [degree];Entries/2 deg");

  // Legends
  TLegend *l_mom = new TLegend(0.6, 0.5, 0.85, 0.85);
  TLegend *l_cos = new TLegend(0.6, 0.5, 0.85, 0.85);
  TLegend *l_deg = new TLegend(0.6, 0.5, 0.85, 0.85);

  // Canvases
  TCanvas *c_mom = new TCanvas("c_mom","", 800, 600);
  TCanvas *c_cos = new TCanvas("c_cos","", 800, 600);
  TCanvas *c_deg = new TCanvas("c_deg","", 800, 600);
  c_mom->SetGrid(1); c_cos->SetGrid(1); c_deg->SetGrid(1);

  TH1D *hist_mom[num_ninja_mode];
  TH1D *hist_cos[num_ninja_mode];
  TH1D *hist_deg[num_ninja_mode];
  for ( Int_t imode = 0; imode < num_ninja_mode; imode++ ) {
    hist_mom[imode] = (TH1D*)file->Get(Form("hist_mom_%d", imode));
    hist_cos[imode] = (TH1D*)file->Get(Form("hist_cos_%d", imode));
    hist_deg[imode] = (TH1D*)file->Get(Form("hist_deg_%d", imode));
    hist_mom[imode]->SetLineColor(mode_color[imode]);
    hist_cos[imode]->SetLineColor(mode_color[imode]);
    hist_deg[imode]->SetLineColor(mode_color[imode]);
    hist_mom[imode]->SetFillColor(mode_color[imode]);
    hist_cos[imode]->SetFillColor(mode_color[imode]);
    hist_deg[imode]->SetFillColor(mode_color[imode]);
    hist_mom[imode]->SetFillStyle(mode_style[imode]);
    hist_cos[imode]->SetFillStyle(mode_style[imode]);
    hist_deg[imode]->SetFillStyle(mode_style[imode]);
    l_mom->AddEntry(hist_mom[imode], mode_name[imode], "f");
    l_cos->AddEntry(hist_cos[imode], mode_name[imode], "f");
    l_deg->AddEntry(hist_deg[imode], mode_name[imode], "f");
  }

  for ( Int_t imode = 0; imode < num_ninja_mode; imode++ ) {
    hs_mom->Add(hist_mom[mode_stack_order[imode]]);
    hs_cos->Add(hist_cos[mode_stack_order[imode]]);
    hs_deg->Add(hist_deg[mode_stack_order[imode]]);
  }

  TFile *ofile = new TFile("stack_mode.root", "recreate");

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  c_mom->cd();
  hs_mom->Draw("hist e");
  l_mom->Draw();
  c_mom->Update();
  ofile->cd();
  c_mom->Write("hs_mom");

  c_cos->cd();
  hs_cos->Draw("hist e");
  l_cos->Draw();
  c_cos->Update();
  ofile->cd();
  c_cos->Write("hs_cos");

  c_deg->cd();
  hs_deg->Draw("hist e");
  l_deg->Draw();
  c_deg->Update();
  ofile->cd();
  c_deg->Write("hs_deg");

  ofile->Close();

}
