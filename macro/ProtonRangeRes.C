
void ProtonRangeRes() {

  TString filename = "/hsm/nu/ninja/pra_tmp/mc_tmp_20220620/mom_root_file/mom_file.root";
  TFile *file = new TFile(filename, "read");
  TTree *tree = (TTree*)file->Get("tree");

  double weight;
  int true_particle_id;
  int recon_particle_id;
  int stop_flag;
  double angle;
  double true_momentum;
  double recon_mom_ecc_range_proton;

  tree->SetBranchAddress("weight", &weight);
  tree->SetBranchAddress("true_particle_id", &true_particle_id);
  tree->SetBranchAddress("recon_particle_id", &recon_particle_id);
  tree->SetBranchAddress("stop_flag", &stop_flag);
  tree->SetBranchAddress("angle", &angle);
  tree->SetBranchAddress("true_momentum", &true_momentum);
  tree->SetBranchAddress("recon_mom_ecc_range_proton", &recon_mom_ecc_range_proton);

  TH1D *hist[10][15];
  Double_t mom_edge[11];
  Double_t ang_edge[16];
  for ( int i = 0; i < 10; i++ ) {
    double mom_min = i * 100.;
    double mom_max = (i + 1) * 100.;
    mom_edge[i] = mom_min;
    if ( i == 9 ) mom_edge[10] = mom_max;
    for ( int j = 0; j < 15; j++ ) {
      double ang_min, ang_max;
      if ( j < 10 ) {
	ang_min = j * 0.1;
	ang_max = (j + 1) * 0.1;
      }
      else if ( j != 14 ) {
	ang_min = 1. + (j - 10) * 0.5;
	ang_max = 1. + (j - 9) * 0.5;
      }
      else {
	ang_min = 3.0;
	ang_max = 5.0;
      }

      ang_edge[j] = ang_min;
      if ( j == 14 ) ang_edge[15] = ang_max;

      hist[i][j] = new TH1D(Form("hist_%d_%d", i, j),
			    Form("(%.1f < tan#theta < %.1f, %.0f < P_{true} < %.0f);(P_{recon} - P_{true})/P_{true};Entries",
				 ang_min, ang_max, mom_min, mom_max),
			    100, -.1, .1);
    }    
  }

  for ( int ientry = 0; ientry < tree->GetEntries(); ientry++ ) {

    tree->GetEntry(ientry);
    
    if ( ientry % 10000 == 0 ) {
      std::cerr << ientry << "\r";
    }

    if ( true_particle_id != 2212 ) continue;
    if ( true_particle_id != recon_particle_id ) continue;
    if ( stop_flag != 2 ) continue;

    int imom = true_momentum / 100;
    if ( imom >= 10 ) imom = 10;
    int jang = angle / 0.1;
    if ( jang > 10 ) {
      jang = 10 + (angle - 1.) / 0.5;
    }
    if ( jang >= 14 ) {
      jang = 14;
    }

    hist[imom][jang]->Fill((recon_mom_ecc_range_proton - true_momentum) / true_momentum, weight);

  }
  std::cerr << std::endl;

  TH2D *hist_res = new TH2D("hist_res", "Relative resolutsions;P_{true} [MeV/c];tan#theta",
			    10, mom_edge, 15, ang_edge);
  TF1 *gaus[10][15];
  for ( int i = 0; i < 10; i++ ) {
    for ( int j = 0; j < 15; j++ ) {
      gaus[i][j] = new TF1(Form("gaus_%d_%d", i, j), "gaus", -0.1, 0.1);
      hist[i][j]->Fit(gaus[i][j], "Q", "", -0.05, 0.05);
      gaus[i][j]->Draw("same");
      hist_res->Fill((mom_edge[i] + mom_edge[i+1]) / 2.,
		     (ang_edge[j] + ang_edge[j+1]) / 2.,
		     gaus[i][j]->GetParameter(2));
      hist_res->SetBinError(i+1, j+1, gaus[i][j]->GetParError(2));
    }
  }

  TFile *ofile = new TFile("~/range_res.root","recreate");
  ofile->cd();
  for ( int i = 0; i < 10; i++ ) {
    for ( int j = 0; j < 15; j++ ) {
      hist[i][j]->Write();
      gaus[i][j]->Write();
    }
  }
  hist_res->Write();
  ofile->Close();

}
