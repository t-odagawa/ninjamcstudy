void MergeEfficiencyBreakdown() {

  TString filename = "/hsm/nu/ninja/pra_tmp/mc_tmp_20220620/eff_br/efficiency_breakdown_*.root";
  TChain *tree = new TChain("tree", "tree");
  tree->Add(filename);
  double mom_entry_total[12] = {};
  double ang_entry_total[10] = {};
  double mom_ang_entry_total[12][10] = {{}};
  double mom_cut0[12] = {};
  double mom_cut1[12] = {};
  double mom_cut2[12] = {};
  double mom_cut3[12] = {};
  double mom_cut4[12] = {};
  double mom_cut5[12] = {};
  double mom_cut6[12] = {};
  double ang_cut0[10] = {};
  double ang_cut1[10] = {};
  double ang_cut2[10] = {};
  double ang_cut3[10] = {};
  double ang_cut4[10] = {};
  double ang_cut5[10] = {};
  double ang_cut6[10] = {};
  double mom_ang_cut0[12][10] = {{}};
  double mom_ang_cut1[12][10] = {{}};
  double mom_ang_cut2[12][10] = {{}};
  double mom_ang_cut3[12][10] = {{}};
  double mom_ang_cut4[12][10] = {{}};
  double mom_ang_cut5[12][10] = {{}};
  double mom_ang_cut6[12][10] = {{}};

  tree->SetBranchAddress("mom_entry_total", mom_entry_total);
  tree->SetBranchAddress("ang_entry_total", ang_entry_total);
  tree->SetBranchAddress("mom_ang_entry_total", mom_ang_entry_total);
  tree->SetBranchAddress("mom_cut0", mom_cut0);
  tree->SetBranchAddress("mom_cut1", mom_cut1);
  tree->SetBranchAddress("mom_cut2", mom_cut2);
  tree->SetBranchAddress("mom_cut3", mom_cut3);
  tree->SetBranchAddress("mom_cut4", mom_cut4);
  tree->SetBranchAddress("mom_cut5", mom_cut5);
  tree->SetBranchAddress("mom_cut6", mom_cut5);
  tree->SetBranchAddress("ang_cut0", ang_cut0);
  tree->SetBranchAddress("ang_cut1", ang_cut1);
  tree->SetBranchAddress("ang_cut2", ang_cut2);
  tree->SetBranchAddress("ang_cut3", ang_cut3);
  tree->SetBranchAddress("ang_cut4", ang_cut4);
  tree->SetBranchAddress("ang_cut5", ang_cut5);
  tree->SetBranchAddress("ang_cut6", ang_cut5);
  tree->SetBranchAddress("mom_ang_cut0", mom_ang_cut0);
  tree->SetBranchAddress("mom_ang_cut1", mom_ang_cut1);
  tree->SetBranchAddress("mom_ang_cut2", mom_ang_cut2);
  tree->SetBranchAddress("mom_ang_cut3", mom_ang_cut3);
  tree->SetBranchAddress("mom_ang_cut4", mom_ang_cut4);
  tree->SetBranchAddress("mom_ang_cut5", mom_ang_cut5);
  tree->SetBranchAddress("mom_ang_cut6", mom_ang_cut5);

  std::cout << "Entries : " << tree->GetEntries() << std::endl;

  double mom_n[12] = {};
  double ang_n[10] = {};
  double mom_ang_n[12][10] = {{}};
  double mom_cut0_value[12] = {};
  double mom_cut1_value[12] = {};
  double mom_cut2_value[12] = {};
  double mom_cut3_value[12] = {};
  double mom_cut4_value[12] = {};
  double mom_cut5_value[12] = {};
  double mom_cut6_value[12] = {};
  double ang_cut0_value[12] = {};
  double ang_cut1_value[12] = {};
  double ang_cut2_value[12] = {};
  double ang_cut3_value[12] = {};
  double ang_cut4_value[12] = {};
  double ang_cut5_value[12] = {};
  double ang_cut6_value[12] = {};
  double mom_ang_cut0_value[12][10] = {{}};
  double mom_ang_cut1_value[12][10] = {{}};
  double mom_ang_cut2_value[12][10] = {{}};
  double mom_ang_cut3_value[12][10] = {{}};
  double mom_ang_cut4_value[12][10] = {{}};
  double mom_ang_cut5_value[12][10] = {{}};
  double mom_ang_cut6_value[12][10] = {{}};

  double mom_bin_value[12];
  double ang_bin_value[10];
  double mom_bin_err[12];
  double ang_bin_err[10];
  for ( int i = 0; i < 10; i++ ) {
    mom_bin_value[i] = ( i + 0.5 ) * 100.;
    mom_bin_err[i] = 50.;
  }
  mom_bin_value[10] = 1250.;
  mom_bin_err[10] = 250.;
  mom_bin_value[11] = 1750.;
  mom_bin_err[11] = 250.;

  for ( int i = 0; i < 9; i++ ) {
    ang_bin_value[i] = ( i + 0.5 ) * 5.;
    ang_bin_err[i] = 2.5;
  }
  ang_bin_value[9] = 67.5;
  ang_bin_err[9] = 22.5;
  
  
  for ( int i = 0; i < tree->GetEntries(); i++ ) {
    
    tree->GetEntry(i);

    for ( int imom = 0; imom < 12; imom++ ) {
      mom_n[imom] += mom_entry_total[imom];
      mom_cut0_value[imom] += mom_cut0[imom];
      mom_cut1_value[imom] += mom_cut1[imom];
      mom_cut2_value[imom] += mom_cut2[imom];
      mom_cut3_value[imom] += mom_cut3[imom];
      mom_cut4_value[imom] += mom_cut4[imom];
      mom_cut5_value[imom] += mom_cut5[imom];
      mom_cut6_value[imom] += mom_cut6[imom];
    }
    for ( int iang = 0; iang < 10; iang++ ) {
      ang_n[iang] += ang_entry_total[iang];
      ang_cut0_value[iang] += ang_cut0[iang];
      ang_cut1_value[iang] += ang_cut1[iang];
      ang_cut2_value[iang] += ang_cut2[iang];
      ang_cut3_value[iang] += ang_cut3[iang];
      ang_cut4_value[iang] += ang_cut4[iang];
      ang_cut5_value[iang] += ang_cut5[iang];
      ang_cut6_value[iang] += ang_cut6[iang];
    }
    for ( int imom = 0; imom < 12; imom++ ) {
      for ( int jang = 0; jang < 10; jang++ ) {
	mom_ang_n[imom][jang] += mom_ang_entry_total[imom][jang];
	mom_ang_cut0_value[imom][jang] += mom_ang_cut0[imom][jang];
	mom_ang_cut1_value[imom][jang] += mom_ang_cut1[imom][jang];
	mom_ang_cut2_value[imom][jang] += mom_ang_cut2[imom][jang];
	mom_ang_cut3_value[imom][jang] += mom_ang_cut3[imom][jang];
	mom_ang_cut4_value[imom][jang] += mom_ang_cut4[imom][jang];
	mom_ang_cut5_value[imom][jang] += mom_ang_cut5[imom][jang];	
	mom_ang_cut6_value[imom][jang] += mom_ang_cut6[imom][jang];
      }
    }
  }

  double mom_entry_all = 0.;
  double mom_cut0_value_all = 0.;
  double mom_cut1_value_all = 0.;
  double mom_cut2_value_all = 0.;
  double mom_cut3_value_all = 0.;
  double mom_cut4_value_all = 0.;
  double mom_cut5_value_all = 0.;
  double mom_cut6_value_all = 0.;
  for ( int i = 0; i < 12; i++ ) {
    mom_entry_all += mom_n[i];
    mom_cut0_value_all += mom_cut0_value[i];
    mom_cut1_value_all += mom_cut1_value[i];
    mom_cut2_value_all += mom_cut2_value[i];
    mom_cut3_value_all += mom_cut3_value[i];
    mom_cut4_value_all += mom_cut4_value[i];
    mom_cut5_value_all += mom_cut5_value[i];
    mom_cut6_value_all += mom_cut6_value[i];
  }

  TH1D *mom_eff_history = new TH1D("mom_eff_history", ";;Efficiency", 7, 0, 7);
  mom_eff_history->GetYaxis()->SetRangeUser(0., 1.1);
  mom_eff_history->GetXaxis()->SetBinLabel(1, "No cut");
  mom_eff_history->GetXaxis()->SetBinLabel(2, "Cut 1");
  mom_eff_history->GetXaxis()->SetBinLabel(3, "Cut 2");
  mom_eff_history->GetXaxis()->SetBinLabel(4, "Cut 3");
  mom_eff_history->GetXaxis()->SetBinLabel(5, "Cut 4");
  mom_eff_history->GetXaxis()->SetBinLabel(6, "Cut 5");
  mom_eff_history->GetXaxis()->SetBinLabel(7, "Cut 6");

  mom_eff_history->Fill(0.5, mom_cut0_value_all / mom_cut0_value_all);
  mom_eff_history->Fill(1, mom_cut1_value_all / mom_cut0_value_all);
  mom_eff_history->Fill(2, mom_cut2_value_all / mom_cut0_value_all);
  mom_eff_history->Fill(3, mom_cut3_value_all / mom_cut0_value_all);
  mom_eff_history->Fill(4, mom_cut4_value_all / mom_cut0_value_all);
  mom_eff_history->Fill(5, mom_cut5_value_all / mom_cut0_value_all);
  mom_eff_history->Fill(6, mom_cut6_value_all / mom_cut0_value_all);
  mom_eff_history->SetBinError(1, 0.);
  mom_eff_history->SetBinError(2, std::sqrt(mom_cut1_value_all / mom_cut0_value_all
					    * ( 1 - mom_cut1_value_all / mom_cut0_value_all))
			       / mom_entry_all);
  mom_eff_history->SetBinError(3, std::sqrt(mom_cut2_value_all / mom_cut0_value_all
					    * ( 1 - mom_cut2_value_all / mom_cut0_value_all))
			       / mom_entry_all);
  mom_eff_history->SetBinError(4, std::sqrt(mom_cut3_value_all / mom_cut0_value_all
					    * ( 1 - mom_cut3_value_all / mom_cut0_value_all))
			       / mom_entry_all);
  mom_eff_history->SetBinError(5, std::sqrt(mom_cut4_value_all / mom_cut0_value_all
					    * ( 1 - mom_cut4_value_all / mom_cut0_value_all))
			       / mom_entry_all);
  mom_eff_history->SetBinError(6, std::sqrt(mom_cut5_value_all / mom_cut0_value_all
					    * ( 1 - mom_cut5_value_all / mom_cut0_value_all))
			       / mom_entry_all);
  mom_eff_history->SetBinError(7, std::sqrt(mom_cut6_value_all / mom_cut0_value_all
					    * ( 1 - mom_cut6_value_all / mom_cut0_value_all))
			       / mom_entry_all);

  double mom_bin_edge[] = {0., 100., 200., 300., 400., 500., 600., 700., 800., 900., 1000., 1500., 2000.};
  double ang_bin_edge[] = {0., 5., 10., 15., 20., 25., 30., 35., 40., 45., 90.};
  
  TFile *ofile = new TFile("~/merge_eff_br.root", "recreate");

  ofile->cd();
  mom_eff_history->Write();
  ofile->Close();

}
