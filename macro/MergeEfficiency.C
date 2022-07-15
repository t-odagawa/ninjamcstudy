void MergeEfficiency() {

  TString filename = "/hsm/nu/ninja/pra_tmp/mc_tmp_20220620/output/output_mode5_*.root";
  TChain *tree = new TChain("tree", "tree");
  tree->Add(filename);
  double mom_entry_total[12] = {};
  double ang_entry_total[10] = {};
  double mom_ang_entry_total[12][10] = {{}};
  double mom_total[12] = {};
  double ang_total[10] = {};
  double mom_ang_total[12][10] = {{}};
  double mom_eff[12] = {};
  double ang_eff[10] = {};
  double mom_ang_eff[12][10] = {{}};
  
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

  double mom_eff_value[12] = {};
  double ang_eff_value[10] = {};
  double mom_ang_eff_value[12][10] = {{}};
  double mom_total_value[12] = {};
  double ang_total_value[10] = {};
  double mom_ang_total_value[12][10] = {{}};
  double mom_eff_err[12];
  double ang_eff_err[10];
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

  double mom_n[12] = {};
  double ang_n[10] = {};
  double mom_ang_n[12][10]= {{}}; 

  // calculate efficiency in each file
  for ( int i = 0; i < tree->GetEntries(); i++ ) {

    tree->GetEntry(i);

    for ( int imom = 0; imom < 12; imom++ ) {
      mom_n[imom] += mom_entry_total[imom];
      if ( mom_total[imom] > 0. ) {
	mom_eff_value[imom] += mom_eff[imom];
	mom_total_value[imom] += mom_total[imom];
      }
    }
    for ( int iang = 0; iang < 10; iang++ ) {
      ang_n[iang] += ang_entry_total[iang];
      if ( ang_total[iang] > 0. ) {
	ang_eff_value[iang] += ang_eff[iang];
	ang_total_value[iang] += ang_total[iang];
      }
    }
    for ( int imom = 0; imom < 12; imom++ ) {
      for ( int jang = 0; jang < 10; jang++ ) {
	mom_ang_n[imom][jang] += mom_ang_entry_total[imom][jang];
	if ( mom_ang_total[imom][jang] > 0. ) {
	  mom_ang_eff_value[imom][jang] += mom_ang_eff[imom][jang];
	  mom_ang_total_value[imom][jang] += mom_ang_total[imom][jang];
	}
      }
    }
  }

  // averaged value of every file
  for ( int imom = 0; imom < 12; imom++ ) {
    if ( mom_total_value[imom] > 0. )
      mom_eff_value[imom] /= mom_total_value[imom];
    else
      mom_eff_value[imom] = 0.;
  }
  for ( int iang = 0; iang < 10; iang++ ) {
    if ( ang_total_value[iang] > 0. )
      ang_eff_value[iang] /= ang_total_value[iang];
    else
      ang_eff_value[iang] = 0.;
  }
  for ( int imom = 0; imom < 12; imom++ ) {
    for ( int jang = 0; jang < 10; jang++ ) {
      if ( mom_ang_total_value[imom][jang] > 0. )
	mom_ang_eff_value[imom][jang] /= mom_ang_total_value[imom][jang];
      else
	mom_ang_eff_value[imom][jang] = 0.;
    }
  }

  // error is calculated by total statistics
  for ( int imom = 0; imom < 12; imom++ ) {
    mom_eff_err[imom] = std::sqrt(mom_eff_value[imom] * (1-mom_eff_value[imom]) / mom_n[imom]);
  }
  for ( int iang = 0; iang < 10; iang++ ) {
    ang_eff_err[iang] = std::sqrt(ang_eff_value[iang] * (1-ang_eff_value[iang]) / ang_n[iang]);
  }

  TFile *ofile = new TFile("~/merge_eff.root", "recreate");

  TGraphErrors *ge_mom_eff = new TGraphErrors(12, 
					      mom_bin_value, mom_eff_value, 
					      mom_bin_err, mom_eff_err);
  ge_mom_eff->SetName("ge_mom_eff");
  TGraphErrors *ge_ang_eff = new TGraphErrors(10, 
					      ang_bin_value, ang_eff_value,
					      ang_bin_err, ang_eff_err);
  ge_ang_eff->SetName("ge_ang_eff");
  double mom_bin_edge[] = {0., 100., 200., 300., 400., 500., 600., 700., 800., 900., 1000., 1500., 2000.};
  double ang_bin_edge[] = {0., 5., 10., 15., 20., 25., 30., 35., 40., 45., 90.};
  TH2D *hist_mom_ang_eff = new TH2D("hist_mom_ang_eff", "Muon efficiency;Momentum [MeV/c];Angle [deg]",
				    12, mom_bin_edge, 10, ang_bin_edge);
  for ( int i = 0; i < 12; i++ ) {
    for ( int j = 0; j < 10; j++ ) {
      double imom = (mom_bin_edge[i] + mom_bin_edge[i+1]) / 2.;
      double jang = (ang_bin_edge[j] + ang_bin_edge[j+1]) / 2.;
      if ( mom_ang_eff_value[i][j] > 0. )
	hist_mom_ang_eff->Fill(imom, jang, mom_ang_eff_value[i][j]);
    }
  }

  ofile->cd();
  ge_mom_eff->Write();
  ge_ang_eff->Write();
  hist_mom_ang_eff->Write();
  ofile->Close();

}
