#include "/home/t2k/odagawa/NinjaMCStudy/lib/HistogramStyle.hpp"

void StackModeMulti() {

  // Legend
  TLegend *l_water_multi = new TLegend(0.6, 0.5, 0.85, 0.85);
  l_water_multi->SetName("l_water_multi");  
  TLegend *l_water_multi_p = new TLegend(0.6, 0.5, 0.85, 0.85);
  l_water_multi_p->SetName("l_water_multi_p");

  // Signal + packing/mis-pid background
  TString filename = "/hsm/nu/ninja/pra_tmp/mc_tmp_20220620/output/output_mode0.root";
  TFile *file = new TFile(filename, "read");
  auto hist_water_total_multi = (TH1D*)file->Get("hist_water_total_multi");
  double scale = 1. / 976. / 33.156 * 0.47 * 0.99;

  TH1D *hist_water_multi[num_ninja_mode];
  TH1D *hist_water_proton_multi[num_ninja_mode];
  TH1D *hist_water_pion_multi[num_ninja_mode];
  for ( int i = 0; i < num_ninja_mode - 1; i++ ) {
    hist_water_multi[i] = (TH1D*)file->Get(Form("hist_water_mode_multi_%d",i));
    hist_water_multi[i]->Scale(scale);
    hist_water_proton_multi[i] = (TH1D*)file->Get(Form("hist_water_mode_proton_multi_%d",i));
    hist_water_proton_multi[i]->Scale(scale);
    hist_water_pion_multi[i] = (TH1D*)file->Get(Form("hist_water_mode_pion_multi_%d",i));
    hist_water_pion_multi[i]->Scale(scale);
    l_water_multi->AddEntry(hist_water_multi[i], mode_name[i], "f");
    l_water_multi_p->AddEntry(hist_water_proton_multi[i], mode_name[i], "f");
  }


  // Packing background
  auto hist_pack_bg = new TH1D("hist_pack_bg", "", 10, 0.5, 10.5);
  hist_pack_bg->Fill(1, hist_water_total_multi->Integral(0,1) * 0.09);
  hist_pack_bg->Scale(scale);

  // Muon mis-id CC + NC
  auto hist_water_multi_muon_misid = (TH1D*)file->Get("hist_water_multi_muon_misid");
  auto hist_water_proton_multi_muon_misid = (TH1D*)file->Get("hist_water_proton_multi_muon_misid");
  auto hist_water_pion_multi_muon_misid = (TH1D*)file->Get("hist_water_pion_multi_muon_misid");
  hist_water_multi_muon_misid->Scale(scale);
  hist_water_proton_multi_muon_misid->Scale(scale);
  hist_water_pion_multi_muon_misid->Scale(scale);

  // Partner mis-id/non-id
  auto hist_water_proton_misid_multi = (TH1D*)file->Get("hist_water_proton_misid_multi");
  auto hist_water_pion_misid_multi = (TH1D*)file->Get("hist_water_pion_misid_multi");
  hist_water_proton_misid_multi->Scale(scale);
  hist_water_proton_misid_multi->SetFillColor(kOrange);
  hist_water_pion_misid_multi->Scale(scale);
  hist_water_pion_misid_multi->SetFillColor(kOrange);

  // Proton Module background
  TString pmfilename = "/hsm/nu/ninja/pra_tmp/mc_tmp_pm_20220627/output/bg_dist.root";
  TFile *pmfile = new TFile(pmfilename, "read");
  auto pmhist = (TH1D*)pmfile->Get("hist_multi");
  auto pmhist_p = (TH1D*)pmfile->Get("hist_multi_p");
  auto pmhist_pi = (TH1D*)pmfile->Get("hist_multi_pi");
  double pmscale = 1. / 999. * 0.47 * 0.99;
  pmhist->Scale(pmscale);
  pmhist_p->Scale(pmscale);
  pmhist_pi->Scale(pmscale);

  // Upstream WAGASCI background
  TString uwgfilename = "/hsm/nu/ninja/pra_tmp/mc_tmp_uwg_20220722/output/bg_dist.root";
  TFile *uwgfile = new TFile(uwgfilename, "read");
  auto uwghist = (TH1D*)uwgfile->Get("hist_multi");
  auto uwghist_p = (TH1D*)uwgfile->Get("hist_multi_p");
  auto uwghist_pi = (TH1D*)uwgfile->Get("hist_multi_pi");
  double uwgscale = 1. / 999. * 0.47 * 0.99;
  uwghist->Scale(uwgscale);
  uwghist_p->Scale(uwgscale);
  uwghist_pi->Scale(uwgscale);

  // Downstream WAGASCI backgorund
  TString dwgfilename = "/hsm/nu/ninja/pra_tmp/mc_tmp_dwg_20220712/output/bg_dist.root";
  TFile *dwgfile = new TFile(dwgfilename, "read");
  auto dwghist = (TH1D*)dwgfile->Get("hist_multi");
  auto dwghist_p = (TH1D*)dwgfile->Get("hist_multi_p");
  auto dwghist_pi = (TH1D*)dwgfile->Get("hist_multi_pi");
  double dwgscale = 1. / 999. * 0.47 * 0.99;
  dwghist->Scale(dwgscale);
  dwghist_p->Scale(dwgscale);
  dwghist_pi->Scale(dwgscale);

  // Baby MIND background (negligible)
  // TString bmfilename = "/hsm/nu/ninja/pra_tmp/mc_tmp_bm_20220711/output/bg_dist.root";
  // TFile *bmfile = new TFile(bmfilename, "read");
  // auto bmhist = (TH1D*)bmfile->Get("hist_multi");
  // auto bmhist_p = (TH1D*)bmfile->Get("hist_multi_p");
  // auto bmhist_pi = (TH1D*)bmfile->Get("hist_multi_pi");
  // double bmscale = 1. / 999. * 0.47 * 0.99;
  // bmhist->Scale(bmscale);

  // Other ECC background
  TString otherfilename = "/hsm/nu/ninja/pra_tmp/mc_tmp_20220505/output/bg_dist.root";
  TFile *otherfile = new TFile(otherfilename, "read");
  auto otherhist = (TH1D*)otherfile->Get("hist_multi");
  auto otherhist_p = (TH1D*)otherfile->Get("hist_multi_p");
  auto otherhist_pi = (TH1D*)otherfile->Get("hist_multi_pi");
  double otherscale = 1. / 990. * 0.47 * 0.99;
  otherhist->Scale(otherscale);
  otherhist_p->Scale(otherscale);
  otherhist_pi->Scale(otherscale);

  // Wall background
  TString wallfilename = "/hsm/nu/ninja/pra_tmp/wall_mc_20220616/output/bg_dist.root";
  TFile *wallfile = new TFile(wallfilename, "read");
  auto wallhist = (TH1D*)wallfile->Get("hist_multi");
  auto wallhist_p = (TH1D*)wallfile->Get("hist_multi_p");
  auto wallhist_pi = (TH1D*)wallfile->Get("hist_multi_pi");
  double wallscale = 1. / 9999. * 0.47 * 0.99;
  wallhist->Scale(wallscale);
  wallhist_p->Scale(wallscale);
  wallhist_pi->Scale(wallscale);

  TH1D *hist_bg_ext = new TH1D("hist_bg_ext", "", 10, 0.5, 10.5);
  TList *ext_list = new TList();
  ext_list->Add(pmhist);
  ext_list->Add(uwghist);
  ext_list->Add(dwghist);
  // ext_list->Add(bmhist);
  ext_list->Add(otherhist);
  ext_list->Add(wallhist);
  hist_bg_ext->Merge(ext_list);
  hist_bg_ext->SetFillStyle(3011);
  hist_bg_ext->SetFillColor(kGreen);

  TH1D *hist_bg_ext_p = new TH1D("hist_bg_ext_p", "", 10, -0.5, 9.5);
  TList *ext_list_p = new TList();
  ext_list_p->Add(pmhist_p);
  ext_list_p->Add(uwghist_p);
  ext_list_p->Add(dwghist_p);
  ext_list_p->Add(otherhist_p);
  ext_list_p->Add(wallhist_p);
  hist_bg_ext_p->Merge(ext_list_p);
  hist_bg_ext_p->SetFillStyle(3011);
  hist_bg_ext_p->SetFillColor(kGreen);

  TH1D *hist_bg_ext_pi = new TH1D("hist_bg_ext_pi", "", 10, -0.5, 9.5);
  TList *ext_list_pi = new TList();
  ext_list_pi->Add(pmhist_pi);
  ext_list_pi->Add(uwghist_pi);
  ext_list_pi->Add(dwghist_pi);
  ext_list_pi->Add(otherhist_pi);
  ext_list_pi->Add(wallhist_pi);
  hist_bg_ext_pi->Merge(ext_list_pi);
  hist_bg_ext_pi->SetFillStyle(3011);
  hist_bg_ext_pi->SetFillColor(kGreen);


  // Iron (+ emulsion?) interaction background
  TString fefilename = "/hsm/nu/ninja/pra_tmp/mc_tmp_fe_20220712/output/bg_dist.root";
  TFile *fefile = new TFile(fefilename, "read");
  auto fehist = (TH1D*)fefile->Get("hist_multi");
  auto fehist_p = (TH1D*)fefile->Get("hist_multi_p");
  auto fehist_pi = (TH1D*)fefile->Get("hist_multi_pi");
  double fescale = 1. / 999. / 33.156 * 0.47 * 0.99;
  fescale *= 1.2; // + emulsion
  fehist->Scale(fescale);
  fehist_p->Scale(fescale);
  fehist_pi->Scale(fescale);

  // Anti neutrino interaction background
  TString anufilename = "/hsm/nu/ninja/pra_tmp/mc_tmp_anu_20220714/output/bg_dist.root";
  TFile *anufile = new TFile(anufilename, "read");
  auto anuhist = (TH1D*)anufile->Get("hist_multi");
  auto anuhist_p = (TH1D*)anufile->Get("hist_multi_p");
  auto anuhist_pi = (TH1D*)anufile->Get("hist_multi_pi");
  double anuscale = 1. / 977. / 33.156 * 0.47 * 0.99;
  anuhist->Scale(anuscale);
  anuhist_p->Scale(anuscale);
  anuhist_pi->Scale(anuscale);


  TH1D *hist_bg_int = new TH1D("hist_bg_int", "", 10, 0.5, 10.5);
  TList *int_list = new TList();
  int_list->Add(fehist);
  int_list->Add(anuhist);
  int_list->Add(hist_pack_bg);
  int_list->Add(hist_water_multi_muon_misid);
  hist_bg_int->Merge(int_list);
  hist_bg_int->SetFillStyle(3021);
  hist_bg_int->SetFillColor(kMagenta);

  TH1D *hist_bg_int_p = new TH1D("hist_bg_int_p", "", 10, -0.5, 9.5);
  TList *int_list_p = new TList();
  int_list_p->Add(fehist_p);
  int_list_p->Add(anuhist_p);
  int_list_p->Add(hist_water_proton_multi_muon_misid);
  hist_bg_int_p->Merge(int_list_p);
  hist_bg_int_p->SetFillStyle(3021);
  hist_bg_int_p->SetFillColor(kMagenta);

  TH1D *hist_bg_int_pi = new TH1D("hist_bg_int_pi", "", 10, -0.5, 9.5);
  TList *int_list_pi = new TList();
  int_list_pi->Add(fehist_pi);
  int_list_pi->Add(anuhist_pi);
  int_list_pi->Add(hist_water_pion_multi_muon_misid);
  hist_bg_int_pi->Merge(int_list_pi);
  hist_bg_int_pi->SetFillStyle(3021);
  hist_bg_int_pi->SetFillColor(kMagenta);


  l_water_multi->AddEntry(hist_bg_int, "Internal Beam-related Background", "f");
  l_water_multi->AddEntry(hist_bg_ext, "External Beam-related Background", "f");
  l_water_multi_p->AddEntry(hist_bg_int_p, "Internal Beam-related Background", "f");
  l_water_multi_p->AddEntry(hist_bg_ext_p, "External Beam-related Background", "f");

  // Chance coincidence
  // signal を shift する
  TH1D *hist_bg_cc = new TH1D("hist_bg_cc", "", 10, 0.5, 10.5);
  l_water_multi->AddEntry(hist_bg_cc, "Beam-unrelated Background", "f");
  l_water_multi_p->AddEntry(hist_bg_cc, "Beam-unrelated Background", "f");

  l_water_multi_p->AddEntry(hist_water_proton_misid_multi, "Partner mis-id","f");

  // Data
  TH1D *hist = new TH1D("hist", "", 10, 0.5, 10.5);
  hist->Fill(1,25);
  hist->Fill(2,37);
  hist->Fill(3,22);
  hist->Fill(4,5);
  hist->Fill(5,2);
  hist->Fill(6,3);
  hist->Fill(7,1);
  hist->Sumw2(0);

  TH1D *hist_p = new TH1D("hist_p", "", 10, -0.5, 9.5);
  hist_p->Fill(0., 45);
  hist_p->Fill(1, 41);
  hist_p->Fill(2, 6);
  hist_p->Fill(3, 1);
  hist_p->Fill(4, 2);
  hist_p->Sumw2(0);

  TH1D *hist_pi = new TH1D("hist_pi", "", 10, -0.5, 9.5);
  hist_pi->Fill(0., 53);
  hist_pi->Fill(1, 28);
  hist_pi->Fill(2, 10);
  hist_pi->Fill(3, 3);
  hist_pi->Fill(4, 1);
  hist_pi->Sumw2(0);
  
  l_water_multi->AddEntry(hist, "Data", "lp");
  l_water_multi_p->AddEntry(hist_p, "Data", "lp");


  THStack *hs_water_multi = new THStack("hs_water_multi","Multiplicity;# of tracks;Entries");
  THStack *hs_water_proton_multi = new THStack("hs_water_proton_multi", "Proton multiplicity;# of protons;Entries");
  THStack *hs_water_pion_multi = new THStack("hs_water_pion_multi", "Pion multiplicity;# of pions;Entries");

  hs_water_proton_multi->Add(hist_water_proton_misid_multi);
  hs_water_pion_multi->Add(hist_water_pion_misid_multi);

  hs_water_multi->Add(hist_bg_ext);
  hs_water_multi->Add(hist_bg_int);
  hs_water_proton_multi->Add(hist_bg_ext_p);
  hs_water_proton_multi->Add(hist_bg_int_p);
  hs_water_pion_multi->Add(hist_bg_ext_pi);
  hs_water_pion_multi->Add(hist_bg_int_pi);
  for ( int i = 1; i < num_ninja_mode; i++ ) {
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
  l_water_multi_p->Write();
  hist->Write();
  hist_p->Write();
  hist_pi->Write();
  ofile->Close();

}
