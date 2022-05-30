void StackPid() {

  TString filename = "/hsm/nu/ninja/pra_tmp/mc_tmp_20220505/output/output_mode4.root";
  TFile *file = new TFile(filename, "read");

  THStack *hs_lr = new THStack("hs_lr", "Likelihood ratio;Likelihood ratio;Entries");

  TLegend *l_lr = new TLegend(0.6, 0.5, 0.85, 0.85);
  l_lr->SetName("l_lr");

  auto hist_lr_p = (TH1D*)file->Get("hist_likelihood_ratio_p");
  hist_lr_p->SetFillColor(624);
  l_lr->AddEntry(hist_lr_p, "Proton", "f");
  auto hist_lr_pi = (TH1D*)file->Get("hist_likelihood_ratio_pi");
  hist_lr_pi->SetFillColor(394);
  l_lr->AddEntry(hist_lr_pi, "Pion", "f");

  double scale = 0.47e-3;

  hs_lr->Add(hist_lr_p);
  hs_lr->Add(hist_lr_pi);

  TFile *ofile = new TFile("~/stack_pid.root", "recreate");

  ofile->cd();
  hs_lr->Write();
  l_lr->Write();
  ofile->Close();

}
