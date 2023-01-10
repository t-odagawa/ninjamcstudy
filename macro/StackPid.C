void StackPid() {

  // TString filename = "/hsm/nu/ninja/pra_tmp/mc_tmp_20220620/output/output_mode4.root";
  TString filename = "/group/nu/ninja/work/odagawa/20221020-phd-thesis-preliminary/signal/output/output_mode4.root";
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

  double scale = 1 / 976. / 33.156 * 0.47 * 0.99;

  hist_lr_p->Scale(scale);
  hist_lr_pi->Scale(scale);

  double proton_eff = hist_lr_p->Integral(1,10);
  double proton_mis = hist_lr_p->Integral(11,20);
  double pion_eff = hist_lr_pi->Integral(11,20);
  double pion_mis = hist_lr_pi->Integral(1,10);

  std::cout << "Proton efficiency : " << proton_eff / (proton_eff + proton_mis) << std::endl;
  std::cout << "Proton purity : " << proton_eff / (proton_eff + pion_mis) << std::endl;
  std::cout << "Pion efficiency : " << pion_eff / (pion_eff + pion_mis) << std::endl;
  std::cout << "Pion purity : " << pion_eff / (pion_eff + proton_mis) << std:: endl;

  hs_lr->Add(hist_lr_p);
  hs_lr->Add(hist_lr_pi);

  TFile *ofile = new TFile("~/stack_pid.root", "recreate");

  ofile->cd();
  hs_lr->Write();
  l_lr->Write();
  ofile->Close();

}
