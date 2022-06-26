#include "/home/t2k/odagawa/NinjaMCStudy/lib/HistogramStyle.hpp"

void StackMode0pi2p() {

  TString filename = "/hsm/nu/ninja/pra_tmp/mc_tmp_20220505/output/output_mode7.root";
  TFile *file = new TFile(filename, "read");

  THStack *hs_muon_mom = new THStack("hs_muon_mom", "Muon momentum;p_{#mu} [MeV/c];Entries");
  THStack *hs_proton_mom = new THStack("hs_proton_mom", "Proton momentum;p_{p} [MeV/c];Entries");
  THStack *hs_proton_mom_high = new THStack("hs_proton_mom_high","Higher proton momentum;p_{p1} [MeV/c];Entries");
  THStack *hs_proton_mom_low = new THStack("hs_proton_mom_low", "Lower proton momentum;p_{p2} [MeV/c];Entries");
  THStack *hs_muon_ang = new THStack("hs_muon_ang", "Muon angle;#theta_{#mu} [deg];Entries");
  THStack *hs_proton_ang = new THStack("hs_proton_ang", "Proton angle;#theta_{p} [deg];Entries");
  THStack *hs_proton_ang_high = new THStack("hs_proton_ang_high","Higher proton angle;#theta_{p1} [deg];Entries");
  THStack *hs_proton_ang_low = new THStack("hs_proton_ang_low", "Lower proton angle;#theta_{p2} [deg];Entries");
  THStack *hs_muon_cos = new THStack("hs_muon_cos", "Muon angle;cos#theta_{#mu};Entries");
  THStack *hs_proton_cos = new THStack("hs_proton_cos", "Proton angle;cos#theta_{p};Entries");
  THStack *hs_proton_cos_high = new THStack("hs_proton_cos_high","Higher proton angle;cos#theta_{p1};Entries");
  THStack *hs_proton_cos_low = new THStack("hs_proton_cos_low", "Lower proton angle;cos#theta_{p2};Entries");
  THStack *hs_open_ang = new THStack("hs_open_ang", "Proton-proton opening angle;#theta_{pp} [deg];Entries");
  THStack *hs_open_cos = new THStack("hs_open_cos", "Proton-proton opening angle;cos#theta_{pp};Entries");
  THStack *hs_mom_ratio = new THStack("hs_mom_ratio", "Proton-proton momentum ratio;p_{p2}/p_{p1};Entries");
  THStack *hs_mom_vecsum = new THStack("hs_mom_vecsum", "Proton momentum vector sum;|#vec{p}_{p1} + #vec{p}_{p2}| [MeV/c];Entries");
  THStack *hs_mom_scasum = new THStack("hs_mom_scasum", "Proton momentum scalar sum;p_{p1} + p_{p2} [MeV/c];Entries");
  THStack *hs_dptt = new THStack("hs_dptt", "#deltap_{TT};#deltap_{TT} [MeV/c];Entries");
  THStack *hs_dpt = new THStack("hs_dpt", "#deltap_{T};#deltap_{T} [MeV/c];Entries");
  THStack *hs_pn = new THStack("hs_pn", "p_{N};p_{N} [MeV/c];Entries");
  THStack *hs_dalphat = new THStack("hs_dalphat", "#delta#alpha_{T};#delta#alpha_{T} [deg];Entries");
  THStack *hs_cosdat = new THStack("hs_cosdat", "cos#delta#alpha_{T};cos#delta#alpha_{T};Entries");

  THStack *hs_proton_mom_high_low = new THStack("hs_proton_mom_high_low", "Proton momentum;p_{p} [MeV/c];Entries");
  THStack *hs_proton_ang_high_low = new THStack("hs_proton_ang_high_low", "Proton angle;#theta_{p} [deg];Entries");
  THStack *hs_proton_cos_high_low = new THStack("hs_proton_cos_high_low", "Proton angle;cos#theta_{p};Entries");

  TLegend *l_muon_mom = new TLegend(0.6, 0.5, 0.85, 0.85);
  l_muon_mom->SetName("l_muon_mom");

  TLegend *l_high_low = new TLegend(0.6, 0.5, 0.85, 0.85);
  l_high_low->SetName("l_high_low");

  TH1D *hist_muon_mom[num_ninja_mode];
  TH1D *hist_proton_mom[num_ninja_mode];
  TH1D *hist_proton_mom_high[num_ninja_mode];
  TH1D *hist_proton_mom_low[num_ninja_mode];
  TH1D *hist_muon_ang[num_ninja_mode];
  TH1D *hist_proton_ang[num_ninja_mode];
  TH1D *hist_proton_ang_high[num_ninja_mode];
  TH1D *hist_proton_ang_low[num_ninja_mode];
  TH1D *hist_muon_cos[num_ninja_mode];
  TH1D *hist_proton_cos[num_ninja_mode];
  TH1D *hist_proton_cos_high[num_ninja_mode];
  TH1D *hist_proton_cos_low[num_ninja_mode];

  TH1D *hist_open_ang[num_ninja_mode];
  TH1D *hist_open_cos[num_ninja_mode];
  TH1D *hist_mom_ratio[num_ninja_mode];
  TH1D *hist_mom_vecsum[num_ninja_mode];
  TH1D *hist_mom_scasum[num_ninja_mode];

  TH1D *hist_dptt[num_ninja_mode];
  TH1D *hist_dpt[num_ninja_mode];
  TH1D *hist_pn[num_ninja_mode];
  TH1D *hist_dalphat[num_ninja_mode];
  TH1D *hist_cosdat[num_ninja_mode];

  for ( int i = 0; i < num_ninja_mode; i++ ) {
    hist_muon_mom[i] = (TH1D*)file->Get(Form("hist_muon_mom_%d", i));
    hist_proton_mom[i] = (TH1D*)file->Get(Form("hist_proton_mom_%d", i));
    hist_proton_mom_high[i] = (TH1D*)file->Get(Form("hist_proton_mom_high_%d", i));
    hist_proton_mom_low[i] = (TH1D*)file->Get(Form("hist_proton_mom_low_%d", i));
    hist_muon_ang[i] = (TH1D*)file->Get(Form("hist_muon_ang_%d", i));
    hist_proton_ang[i] = (TH1D*)file->Get(Form("hist_proton_ang_%d", i));
    hist_proton_ang_high[i] = (TH1D*)file->Get(Form("hist_proton_ang_high_%d", i));
    hist_proton_ang_low[i] = (TH1D*)file->Get(Form("hist_proton_ang_low_%d", i));
    hist_muon_cos[i] = (TH1D*)file->Get(Form("hist_muon_cos_%d", i));
    hist_proton_cos[i] = (TH1D*)file->Get(Form("hist_proton_cos_%d", i));
    hist_proton_cos_high[i] = (TH1D*)file->Get(Form("hist_proton_cos_high_%d", i));
    hist_proton_cos_low[i] = (TH1D*)file->Get(Form("hist_proton_cos_low_%d", i));

    hist_open_ang[i] = (TH1D*)file->Get(Form("hist_open_ang_%d", i));
    hist_open_cos[i] = (TH1D*)file->Get(Form("hist_open_cos_%d", i));
    hist_mom_ratio[i] = (TH1D*)file->Get(Form("hist_mom_ratio_%d", i));
    hist_mom_vecsum[i] = (TH1D*)file->Get(Form("hist_mom_vecsum_%d", i));
    hist_mom_scasum[i] = (TH1D*)file->Get(Form("hist_mom_scasum_%d", i));

    hist_dptt[i] = (TH1D*)file->Get(Form("hist_dptt_%d", i));
    hist_dpt[i] = (TH1D*)file->Get(Form("hist_dpt_%d", i));
    hist_pn[i] = (TH1D*)file->Get(Form("hist_pn_%d", i));
    hist_dalphat[i] = (TH1D*)file->Get(Form("hist_dalphat_%d", i));
    hist_cosdat[i] = (TH1D*)file->Get(Form("hist_cosdat_%d", i));
    
    l_muon_mom->AddEntry(hist_muon_mom[i], mode_name[i], "f");
  }
  
  for ( int i = 0; i < num_ninja_mode; i++ ) {
    hs_muon_mom->Add(hist_muon_mom[mode_stack_order[i]]);
    hs_proton_mom->Add(hist_proton_mom[mode_stack_order[i]]);
    hs_proton_mom_high->Add(hist_proton_mom_high[mode_stack_order[i]]);
    hs_proton_mom_low->Add(hist_proton_mom_low[mode_stack_order[i]]);
    hs_muon_ang->Add(hist_muon_ang[mode_stack_order[i]]);
    hs_proton_ang->Add(hist_proton_ang[mode_stack_order[i]]);
    hs_proton_ang_high->Add(hist_proton_ang_high[mode_stack_order[i]]);
    hs_proton_ang_low->Add(hist_proton_ang_low[mode_stack_order[i]]);
    hs_muon_cos->Add(hist_muon_cos[mode_stack_order[i]]);
    hs_proton_cos->Add(hist_proton_cos[mode_stack_order[i]]);
    hs_proton_cos_high->Add(hist_proton_cos_high[mode_stack_order[i]]);
    hs_proton_cos_low->Add(hist_proton_cos_low[mode_stack_order[i]]);

    hs_open_ang->Add(hist_open_ang[mode_stack_order[i]]);
    hs_open_cos->Add(hist_open_cos[mode_stack_order[i]]);
    hs_mom_ratio->Add(hist_mom_ratio[mode_stack_order[i]]);
    hs_mom_vecsum->Add(hist_mom_vecsum[mode_stack_order[i]]);
    hs_mom_scasum->Add(hist_mom_scasum[mode_stack_order[i]]);

    hs_dptt->Add(hist_dptt[mode_stack_order[i]]);
    hs_dpt->Add(hist_dpt[mode_stack_order[i]]);
    hs_pn->Add(hist_pn[mode_stack_order[i]]);
    hs_dalphat->Add(hist_dalphat[mode_stack_order[i]]);
    hs_cosdat->Add(hist_cosdat[mode_stack_order[i]]);
  }

  TH1D *hist_proton_mom_high_tot = (TH1D*)file->Get("hist_proton_mom_high");
  TH1D *hist_proton_mom_low_tot = (TH1D*)file->Get("hist_proton_mom_low");
  TH1D *hist_proton_ang_high_tot = (TH1D*)file->Get("hist_proton_ang_high");
  TH1D *hist_proton_ang_low_tot = (TH1D*)file->Get("hist_proton_ang_low");
  TH1D *hist_proton_cos_high_tot = (TH1D*)file->Get("hist_proton_cos_high");
  TH1D *hist_proton_cos_low_tot = (TH1D*)file->Get("hist_proton_cos_low");
  hist_proton_mom_high_tot->SetFillColor(kRed);
  hist_proton_mom_low_tot->SetFillColor(kBlue);
  hist_proton_ang_high_tot->SetFillColor(kRed);
  hist_proton_ang_low_tot->SetFillColor(kBlue);
  hist_proton_cos_high_tot->SetFillColor(kRed);
  hist_proton_cos_low_tot->SetFillColor(kBlue);  
  
  hs_proton_mom_high_low->Add(hist_proton_mom_low_tot);
  hs_proton_mom_high_low->Add(hist_proton_mom_high_tot);
  hs_proton_ang_high_low->Add(hist_proton_ang_low_tot);
  hs_proton_ang_high_low->Add(hist_proton_ang_high_tot);
  hs_proton_cos_high_low->Add(hist_proton_cos_low_tot);
  hs_proton_cos_high_low->Add(hist_proton_cos_high_tot);
  l_high_low->AddEntry(hist_proton_mom_high_tot, "Higher proton", "f");
  l_high_low->AddEntry(hist_proton_mom_low_tot, "Lower proton", "f");

  TFile *ofile = new TFile("~/stack_0pi2p.root", "recreate");

  ofile->cd();
  hs_muon_mom->Write();
  hs_proton_mom->Write();
  hs_proton_mom_high->Write();
  hs_proton_mom_low->Write();
  hs_muon_ang->Write();
  hs_proton_ang->Write();
  hs_proton_ang_high->Write();
  hs_proton_ang_low->Write();
  hs_muon_cos->Write();
  hs_proton_cos->Write();
  hs_proton_cos_high->Write();
  hs_proton_cos_low->Write();
  hs_open_ang->Write();
  hs_open_cos->Write();
  hs_mom_vecsum->Write();
  hs_mom_scasum->Write();
  hs_dptt->Write();
  hs_dpt->Write();
  hs_pn->Write();
  hs_dalphat->Write();
  hs_cosdat->Write();

  hs_proton_mom_high_low->Write();
  hs_proton_ang_high_low->Write();
  hs_proton_cos_high_low->Write();

  l_muon_mom->Write();
  l_high_low->Write();

  ofile->Close();

}
