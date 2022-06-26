void MergeEfficiency() {

  TString filename = "/hsm/nu/ninja/pra_tmp/mc_tmp_20220620/output/output_mode5_*.root";
  TChain *tree = new TChain("tree", "tree");
  tree->Add(filename);
  double mom_entry_total[10];
  double ang_entry_total[10];
  double mom_ang_entry_total[10][10];
  double mom_total[10];
  double ang_total[10];
  double mom_ang_total[10][10];
  double mom_eff[10];
  double ang_eff[10];
  double mom_ang_eff[10][10];
  
  tree->SetBranchAddress("mom_entry_total", mom_entry_total);
  tree->SetBranchAddress("ang_entry_total", ang_entry_total);
  tree->SetBranchAddress("mom_ang_entry_total", mom_ang_entry_total);
  tree->SetBranchAddress("mom_total", mom_total);
  tree->SetBranchAddress("ang_total", ang_total);
  tree->SetBranchAddress("mom_ang_total", mom_ang_total);
  tree->SetBranchAddress("mom_eff", mom_eff);
  tree->SetBranchAddress("ang_eff", ang_eff);
  tree->SetBranchAddress("mom_ang_eff", mom_ang_eff);


  std::cout << "Entries : " << tree->GetEntries() << std::endl;

  double mom_eff_value[10];
  double ang_eff_value[10];
  double mom_ang_eff_value[10][10];
  double mom_eff_err[10];
  double ang_eff_err[10];
  double mom_bin_value[10];
  double ang_bin_value[10];
  double mom_bin_err[10];
  double ang_bin_err[10];
  for ( int i = 0; i < 10; i++ ) {
    mom_bin_value[i] = ( i + 0.5 ) * 100.;
    ang_bin_value[i] = ( i + 0.5 ) * 5.;
    mom_bin_err[i] = 50.;
    ang_bin_err[i] = 2.5;
    if ( i == 9 ) {
      mom_bin_value[i] = 1200.;
      ang_bin_value[i] = 67.5;
      mom_bin_err[i] = 300.;
      ang_bin_err[i] = 22.5;
    }
  }

  double mom_n[10];
  double ang_n[10];
  double mom_ang_n[10][10];

  // calculate efficiency in each file
  for ( int i = 0; i < tree->GetEntries(); i++ ) {

    tree->GetEntry(i);

    for ( int ibin = 0; ibin < 10; ibin++ ) {
      mom_n[ibin] += mom_entry_total[ibin];
      ang_n[ibin] += ang_entry_total[ibin];
      
      mom_eff_value[ibin] += mom_eff[ibin] / mom_total[ibin];
      ang_eff_value[ibin] += ang_eff[ibin] / ang_total[ibin];

      for ( int jbin = 0; jbin < 10; jbin++ ) {
	mom_ang_n[ibin][jbin] += mom_ang_entry_total[ibin][jbin];
	mom_ang_eff_value[ibin][jbin] += mom_ang_eff[ibin][jbin] / mom_ang_total[ibin][jbin];
      }
    }
  }

  // averaged value of every file
  for ( int ibin = 0; ibin < 10; ibin++ ) {
    mom_eff_value[ibin] /= tree->GetEntries();
    ang_eff_value[ibin] /= tree->GetEntries();
    for ( int jbin = 0; jbin < 10; jbin++ ) {
      mom_ang_eff_value[ibin][jbin] /= tree->GetEntries();
    }
  }

  // error is calculated by total statistics
  for ( int ibin = 0; ibin < 10; ibin++ ) {
    mom_eff_err[ibin] = std::sqrt(mom_eff_value[ibin] * (1-mom_eff_value[ibin]) / mom_n[ibin]);
    ang_eff_err[ibin] = std::sqrt(ang_eff_value[ibin] * (1-ang_eff_value[ibin]) / ang_n[ibin]);
  }

  TFile *ofile = new TFile("~/merge_eff.root", "recreate");

  TGraphErrors *ge_mom_eff = new TGraphErrors(10, 
					      mom_bin_value, mom_eff_value, 
					      mom_bin_err, mom_eff_err);
  ge_mom_eff->SetName("ge_mom_eff");
  TGraphErrors *ge_ang_eff = new TGraphErrors(10, 
					      ang_bin_value, ang_eff_value,
					      ang_bin_err, ang_eff_err);
  ge_ang_eff->SetName("ge_ang_eff");
  double mom_bin_edge[] = {0., 100., 200., 300., 400., 500., 600., 700., 800., 900., 1500.};
  double ang_bin_edge[] = {0., 5., 10., 15., 20., 25., 30., 35., 40., 45., 90.};
  TH2D *hist_mom_ang_eff = new TH2D("hist_mom_ang_eff", "Muon efficiency;Momentum [MeV/c];Angle [deg]",
				    10, mom_bin_edge, 10, ang_bin_edge);
  for ( int i = 0; i < 10; i++ ) {
    for ( int j = 0; j < 10; j++ ) {
      double imom = (mom_bin_edge[i] + mom_bin_edge[i+1]) / 2.;
      double jang = (ang_bin_edge[j] + ang_bin_edge[j+1]) / 2.;
      hist_mom_ang_eff->Fill(imom, jang, mom_ang_eff_value[i][j]);
    }
  }

  ofile->cd();
  ge_mom_eff->Write();
  ge_ang_eff->Write();
  hist_mom_ang_eff->Write();
  ofile->Close();

}
