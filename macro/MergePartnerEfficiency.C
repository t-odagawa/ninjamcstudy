void MergePartnerEfficiency() {
  
  TString filename = "/hsm/nu/ninja/pra_tmp/mc_tmp_20220620/output/output_mode8_*.root";
  TChain *tree = new TChain("tree", "tree");
  tree->Add(filename);
  double partner_mom_entry_total[11] = {};
  double partner_ang_entry_total[36] = {};
  double partner_mom_ang_entry_total[11][36] = {{}};
  double partner_mom_total[11] = {};
  double partner_ang_total[36] = {};
  double partner_mom_ang_total[11][36] = {{}};
  double partner_mom_eff[11] = {};
  double partner_ang_eff[36] = {};
  double partner_mom_ang_eff[11][36] = {{}};
  double proton_mom_entry_total[11] = {};
  double proton_ang_entry_total[36] = {};
  double proton_mom_ang_entry_total[11][36] = {{}};
  double proton_mom_total[11] = {};
  double proton_ang_total[36] = {};
  double proton_mom_ang_total[11][36] = {{}};
  double proton_mom_eff[11] = {};
  double proton_ang_eff[36] = {};
  double proton_mom_ang_eff[11][36] = {{}};
  double pion_mom_entry_total[11] = {};
  double pion_ang_entry_total[36] = {};
  double pion_mom_ang_entry_total[11][36] = {{}};
  double pion_mom_total[11] = {};
  double pion_ang_total[36] = {};
  double pion_mom_ang_total[11][36] = {{}};
  double pion_mom_eff[11] = {};
  double pion_ang_eff[36] = {};
  double pion_mom_ang_eff[11][36] = {{}};

  tree->SetBranchAddress("partner_mom_entry_total", partner_mom_entry_total);
  tree->SetBranchAddress("partner_ang_entry_total", partner_ang_entry_total);
  tree->SetBranchAddress("partner_mom_ang_entry_total", partner_mom_ang_entry_total);
  tree->SetBranchAddress("partner_mom_total", partner_mom_total);
  tree->SetBranchAddress("partner_ang_total", partner_ang_total);
  tree->SetBranchAddress("partner_mom_ang_total", partner_mom_ang_total);
  tree->SetBranchAddress("partner_mom_eff", partner_mom_eff);
  tree->SetBranchAddress("partner_ang_eff", partner_ang_eff);
  tree->SetBranchAddress("partner_mom_ang_eff", partner_mom_ang_eff);
  tree->SetBranchAddress("proton_mom_entry_total", proton_mom_entry_total);
  tree->SetBranchAddress("proton_ang_entry_total", proton_ang_entry_total);
  tree->SetBranchAddress("proton_mom_ang_entry_total", proton_mom_ang_entry_total);
  tree->SetBranchAddress("proton_mom_total", proton_mom_total);
  tree->SetBranchAddress("proton_ang_total", proton_ang_total);
  tree->SetBranchAddress("proton_mom_ang_total", proton_mom_ang_total);
  tree->SetBranchAddress("proton_mom_eff", proton_mom_eff);
  tree->SetBranchAddress("proton_ang_eff", proton_ang_eff);
  tree->SetBranchAddress("proton_mom_ang_eff", proton_mom_ang_eff);
  tree->SetBranchAddress("pion_mom_entry_total", pion_mom_entry_total);
  tree->SetBranchAddress("pion_ang_entry_total", pion_ang_entry_total);
  tree->SetBranchAddress("pion_mom_ang_entry_total", pion_mom_ang_entry_total);
  tree->SetBranchAddress("pion_mom_total", pion_mom_total);
  tree->SetBranchAddress("pion_ang_total", pion_ang_total);
  tree->SetBranchAddress("pion_mom_ang_total", pion_mom_ang_total);
  tree->SetBranchAddress("pion_mom_eff", pion_mom_eff);
  tree->SetBranchAddress("pion_ang_eff", pion_ang_eff);
  tree->SetBranchAddress("pion_mom_ang_eff", pion_mom_ang_eff);

  std::cout << "Entries : " << tree->GetEntries() << std::endl;

  double partner_mom_eff_value[11] = {};
  double partner_ang_eff_value[36] = {};
  double partner_mom_ang_eff_value[11][36] = {{}};
  double proton_mom_eff_value[11] = {};
  double proton_ang_eff_value[36] = {};
  double proton_mom_ang_eff_value[11][36] = {{}};
  double pion_mom_eff_value[11] = {};
  double pion_ang_eff_value[36] = {};
  double pion_mom_ang_eff_value[11][36] = {{}};

  double partner_mom_total_value[11] = {};
  double partner_ang_total_value[36] = {};
  double partner_mom_ang_total_value[11][36] = {{}};
  double proton_mom_total_value[11] = {};
  double proton_ang_total_value[36] = {};
  double proton_mom_ang_total_value[11][36] = {{}};
  double pion_mom_total_value[11] = {};
  double pion_ang_total_value[36] = {};
  double pion_mom_ang_total_value[11][36] = {{}};
  
  double partner_mom_eff_err[11] = {};
  double partner_ang_eff_err[36] = {};
  double partner_mom_ang_eff_err[11][36] = {{}};
  double proton_mom_eff_err[11] = {};
  double proton_ang_eff_err[36] = {};
  double proton_mom_ang_eff_err[11][36] = {{}};
  double pion_mom_eff_err[11] = {};
  double pion_ang_eff_err[36] = {};
  double pion_mom_ang_eff_err[11][36] = {{}};
  double mom_bin_value[11];
  double ang_bin_value[36];
  double mom_bin_err[11];
  double ang_bin_err[36];

  for ( int i = 0; i < 10; i++ ) {
    mom_bin_value[i] = ( i + 0.5 ) * 100.;
    mom_bin_err[i] = 50.;
  }
  mom_bin_value[10] = 1250.;
  mom_bin_err[10] = 250.;

  for ( int i = 0; i < 36; i++ ) {
    ang_bin_value[i] = (i + 0.5) * 5.;
    ang_bin_err[i] = 2.5;
  }

  double partner_mom_n[11] = {};
  double partner_ang_n[36] = {};
  double partner_mom_ang_n[11][36] = {{}};
  double proton_mom_n[11] = {};
  double proton_ang_n[36] = {};
  double proton_mom_ang_n[11][36] = {{}};
  double pion_mom_n[11] = {};
  double pion_ang_n[36] = {};
  double pion_mom_ang_n[11][36] = {{}};

  for ( int i = 0; i < tree->GetEntries(); i++ ) {
    tree->GetEntry(i);

    for ( int imom = 0; imom < 11; imom++ ) {
      partner_mom_n[imom] += partner_mom_entry_total[imom];
      proton_mom_n[imom] += proton_mom_entry_total[imom];
      pion_mom_n[imom] += pion_mom_entry_total[imom];
      if ( partner_mom_total[imom] > 0 ) { 
	partner_mom_eff_value[imom] += partner_mom_eff[imom];
	partner_mom_total_value[imom] += partner_mom_total[imom];
      }
      if ( proton_mom_total[imom] > 0 ) {
	proton_mom_eff_value[imom] += proton_mom_eff[imom];
	proton_mom_total_value[imom] += proton_mom_total[imom];
      }
      if ( pion_mom_total[imom] > 0 ) {
	pion_mom_eff_value[imom] += pion_mom_eff[imom];
	pion_mom_total_value[imom] += pion_mom_total[imom];
      }
    }

    for ( int iang = 0; iang < 36; iang++ ) {
      partner_ang_n[iang] += partner_ang_entry_total[iang];
      proton_ang_n[iang] += proton_ang_entry_total[iang];
      pion_ang_n[iang] += pion_ang_entry_total[iang];
      if ( partner_ang_total[iang] > 0 ){
	partner_ang_eff_value[iang] += partner_ang_eff[iang];
	partner_ang_total_value[iang] += partner_ang_total[iang];
      }
      if ( proton_ang_total[iang] > 0 ) {
	proton_ang_eff_value[iang] += proton_ang_eff[iang];
	proton_ang_total_value[iang] += proton_ang_total[iang];
      }
      if ( pion_ang_total[iang] > 0 ) {
	pion_ang_eff_value[iang] += pion_ang_eff[iang];
	pion_ang_total_value[iang] += pion_ang_total[iang];
      }
    }
    
    for ( int imom = 0; imom < 11; imom++ ) {
      for ( int jang = 0; jang < 36; jang++ ) {
	partner_mom_ang_n[imom][jang] += partner_mom_ang_entry_total[imom][jang];
	proton_mom_ang_n[imom][jang] += proton_mom_ang_entry_total[imom][jang];
	pion_mom_ang_n[imom][jang] += pion_mom_ang_entry_total[imom][jang];
	if ( partner_mom_ang_total[imom][jang] > 0 ) {
	  partner_mom_ang_eff_value[imom][jang] += partner_mom_ang_eff[imom][jang];
	  partner_mom_ang_total_value[imom][jang] += partner_mom_ang_total[imom][jang];
	}
	if ( proton_mom_ang_total[imom][jang] > 0 ) {
	  proton_mom_ang_eff_value[imom][jang] += proton_mom_ang_eff[imom][jang];
	  proton_mom_ang_total_value[imom][jang] += proton_mom_ang_total[imom][jang];
	}
	if ( pion_mom_ang_total[imom][jang] > 0 ) {
	  pion_mom_ang_eff_value[imom][jang] += pion_mom_ang_eff[imom][jang];
	  pion_mom_ang_total_value[imom][jang] += pion_mom_ang_total[imom][jang];
	}
      }
    }
  }

  // averaged value of every file
  for ( int imom = 0; imom < 11; imom++ ) {
    if ( partner_mom_total_value[imom] > 0. )
      partner_mom_eff_value[imom] /= partner_mom_total_value[imom];
    else 
      partner_mom_eff_value[imom] = 0.;
    if ( proton_mom_total_value[imom] > 0. )
      proton_mom_eff_value[imom] /= proton_mom_total_value[imom];
    else 
      proton_mom_eff_value[imom] = 0.;
    if ( pion_mom_total_value[imom] > 0. )
      pion_mom_eff_value[imom] /= pion_mom_total_value[imom];
    else 
      pion_mom_eff_value[imom] = 0.;
  }
  for ( int iang = 0; iang < 36; iang++ ) {
    if ( partner_ang_total_value[iang] > 0. )
      partner_ang_eff_value[iang] /= partner_ang_total_value[iang];
    else
      partner_ang_eff_value[iang] = 0.;
    if ( proton_ang_total_value[iang] > 0. )
      proton_ang_eff_value[iang] /= proton_ang_total_value[iang];
    else
      proton_ang_eff_value[iang] = 0.;
    if ( pion_ang_total_value[iang] > 0. )
      pion_ang_eff_value[iang] /= pion_ang_total_value[iang];
    else
      pion_ang_eff_value[iang] = 0.;
  }
  for ( int imom = 0; imom < 11; imom++ ) {
    for (int jang = 0; jang < 36; jang++ ) {
      if ( partner_mom_ang_total_value[imom][jang] > 0. )
	partner_mom_ang_eff_value[imom][jang] /= partner_mom_ang_total_value[imom][jang];
      else
	partner_mom_ang_eff_value[imom][jang] = 0.;
      if ( proton_mom_ang_total_value[imom][jang] > 0. )
	proton_mom_ang_eff_value[imom][jang] /= proton_mom_ang_total_value[imom][jang];
      else 
	proton_mom_ang_eff_value[imom][jang] = 0.;
      if ( pion_mom_ang_total_value[imom][jang] > 0. )
	pion_mom_ang_eff_value[imom][jang] /= pion_mom_ang_total_value[imom][jang];
      else
	pion_mom_ang_eff_value[imom][jang] = 0.;
    }
  }
  
  // error calculation
  for ( int imom = 0; imom < 11; imom++ ) {
    if ( partner_mom_n[imom] > 0 )
      partner_mom_eff_err[imom] = std::sqrt(partner_mom_eff_value[imom] * (1 - partner_mom_eff_value[imom]) / partner_mom_n[imom]);
    else
      partner_mom_eff_err[imom] = 0.;
    if ( proton_mom_n[imom] > 0 )
      proton_mom_eff_err[imom] = std::sqrt(proton_mom_eff_value[imom] * (1 - proton_mom_eff_value[imom]) / proton_mom_n[imom]);
    else
      proton_mom_eff_err[imom] = 0.;
    if ( pion_mom_n[imom] > 0 )
      pion_mom_eff_err[imom] = std::sqrt(pion_mom_eff_value[imom] * (1 - pion_mom_eff_value[imom]) / pion_mom_n[imom]);
    else
      pion_mom_eff_err[imom] = 0.;
  }
  for ( int iang = 0; iang < 36; iang++ ) {
    if ( partner_ang_n[iang] > 0 )
      partner_ang_eff_err[iang] = std::sqrt(partner_ang_eff_value[iang] * (1 - partner_ang_eff_value[iang]) / partner_ang_n[iang]);
    else
      partner_ang_eff_err[iang] = 0.;
    if ( proton_ang_n[iang] > 0 )
      proton_ang_eff_err[iang] = std::sqrt(proton_ang_eff_value[iang] * (1 - proton_ang_eff_value[iang]) / proton_ang_n[iang]);
    else
      proton_ang_eff_err[iang] = 0.;
    if ( pion_ang_n[iang] > 0 )
      pion_ang_eff_err[iang] = std::sqrt(pion_ang_eff_value[iang] * (1 - pion_ang_eff_value[iang]) / pion_ang_n[iang]);
    else
      pion_ang_eff_err[iang] = 0.;
  }
  
  TFile *ofile = new TFile("~/merge_partner_eff.root", "recreate");
  
  TGraphErrors *ge_partner_mom_eff = new TGraphErrors(11,
						      mom_bin_value, partner_mom_eff_value,
						      mom_bin_err, partner_mom_eff_err);
  ge_partner_mom_eff->SetName("ge_partner_mom_eff");
  TGraphErrors *ge_proton_mom_eff = new TGraphErrors(11,
						     mom_bin_value, proton_mom_eff_value,
						     mom_bin_err, proton_mom_eff_err);
  ge_proton_mom_eff->SetName("ge_proton_mom_eff");
  TGraphErrors *ge_pion_mom_eff = new TGraphErrors(11,
						   mom_bin_value, pion_mom_eff_value,
						   mom_bin_err, pion_mom_eff_err);
  ge_pion_mom_eff->SetName("ge_pion_mom_eff");

  TGraphErrors *ge_partner_ang_eff = new TGraphErrors(36,
						      ang_bin_value, partner_ang_eff_value,
						      ang_bin_err, partner_ang_eff_err);
  ge_partner_ang_eff->SetName("ge_partner_ang_eff");
  TGraphErrors *ge_proton_ang_eff = new TGraphErrors(36,
						     ang_bin_value, proton_ang_eff_value,
						     ang_bin_err, proton_ang_eff_err);
  ge_proton_ang_eff->SetName("ge_proton_ang_eff");
  TGraphErrors *ge_pion_ang_eff = new TGraphErrors(36,
						   ang_bin_value, pion_ang_eff_value,
						   ang_bin_err, pion_ang_eff_err);
  ge_pion_ang_eff->SetName("ge_pion_ang_eff");

  double mom_bin_edge[] = {0., 100., 200., 300., 400., 500., 600., 700., 800., 900., 1000., 1500.};
  double ang_bin_edge[] = {0., 5., 10., 15., 20., 25.,
			   30., 35., 40., 45., 50., 55.,
			   60., 65., 70., 75., 80., 85.,
			   90., 95., 100., 105., 110., 115.,
			   120., 125., 130., 135., 140., 145.,
			   150., 155., 160., 165., 170., 175., 
			   180.};
  TH2D *hist_partner_mom_ang_eff = new TH2D("hist_partner_mom_ang_eff", "Partner efficiency;Momentum [MeV/c];Angle [deg]",
					    11, mom_bin_edge, 36, ang_bin_edge);
  TH2D *hist_proton_mom_ang_eff = new TH2D("hist_proton_mom_ang_eff", "Proton efficiency;Momentum [MeV/c];Angle [deg]",
					    11, mom_bin_edge, 36, ang_bin_edge);
  TH2D *hist_pion_mom_ang_eff = new TH2D("hist_pion_mom_ang_eff", "Pion efficiency;Momentum [MeV/c];Angle [deg]",
					    11, mom_bin_edge, 36, ang_bin_edge);
  
  for ( int i = 0; i < 11; i++ ) {
    for ( int j = 0; j < 36; j++ ) {
      double imom = (mom_bin_edge[i] + mom_bin_edge[i+1]) / 2.;
      double jang = (ang_bin_edge[j] + ang_bin_edge[j+1]) / 2.;
      if ( partner_mom_ang_total_value[i][j] > 0 )
	hist_partner_mom_ang_eff->Fill(imom, jang, partner_mom_ang_eff_value[i][j]);
      if ( proton_mom_ang_total_value[i][j] > 0 )
	hist_proton_mom_ang_eff->Fill(imom, jang, proton_mom_ang_eff_value[i][j]);
      if ( pion_mom_ang_total_value[i][j] > 0 )
	hist_pion_mom_ang_eff->Fill(imom, jang, pion_mom_ang_eff_value[i][j]);
    }
  }

  ofile->cd();
  ge_partner_mom_eff->Write();
  ge_proton_mom_eff->Write();
  ge_pion_mom_eff->Write();
  ge_partner_ang_eff->Write();
  ge_proton_ang_eff->Write();
  ge_pion_ang_eff->Write();
  hist_partner_mom_ang_eff->Write();
  hist_proton_mom_ang_eff->Write();
  hist_pion_mom_ang_eff->Write();
  ofile->Close();

}
