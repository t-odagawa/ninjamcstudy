#ifndef NINJA_MC_HISTO_STYLE_HPP
#define NINJA_MC_HISTO_STYLE_HPP

#include <TString.h>

const Int_t num_ninja_mode = 6;
const TString mode_name[num_ninja_mode] = {"CCQE", "2p2h", "CC 1#pi", "CC Multi#pi", "CC Other", "NC"};
const Int_t mode_color[num_ninja_mode] = {424, 624, 394, 395, 401, 408};
const Int_t mode_style[num_ninja_mode] = {1001, 1001, 3006, 3005, 1001, 1001};
const Int_t mode_stack_order[num_ninja_mode] = {5, 4, 3, 2, 0, 1};

Int_t GetNinjaModeId(Int_t mode);

const Int_t nu_ene_bin_size = 20;
const double nu_ene_bins[nu_ene_bin_size] = {.09, 0.1, 0.2, 0.3, 0.4,
					     0.5, 0.6, 0.7, 0.8, 1.0,
					     1.2, 1.5, 2.0, 2.5, 3.0,
					     3.5, 4.0, 5.0, 7.0, 10.0};

const Int_t total_multi_bin_size = 7;
const double total_multi_bins[total_multi_bin_size] = {0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 10.5};
const double total_multi_bin_weight[total_multi_bin_size-1] = {1., 1., 1., 1., 1., 1/5.};

const Int_t hadron_multi_bin_size = 7;
const double hadron_multi_bins[hadron_multi_bin_size] = {-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 9.5};
const double hadron_multi_bin_weight[hadron_multi_bin_size-1] = {1., 1., 1., 1., 1., 1/5.};

const Int_t muon_mom_bin_size = 7;
const double muon_mom_bins[muon_mom_bin_size] = {0., 400., 600., 800., 1000., 1200., 2000.};
const double muon_mom_bin_weight[muon_mom_bin_size-1] = {1/2., 1., 1., 1., 1., 1/4.};
/*
const double muon_mom_bins[muon_mom_bin_size] = {0., 100., 200., 300., 400., 500., 600.,
						 700., 800., 900., 1000., 1100., 1200., 1300.,
						 1400., 1500., 1600., 1700., 1800., 1900., 2000.};
*/
const Int_t hadron_mom_bin_size = 6;
const double hadron_mom_bins[hadron_mom_bin_size] = {0., 200., 400., 600., 800., 1000.};
const double hadron_mom_bin_weight[hadron_mom_bin_size-1] = {1., 1., 1., 1., 1.};
/*
const double hadron_mom_bins[hadron_mom_bin_size] = {0., 100., 200., 300., 400., 500., 600., 700.,
						     800., 900., 1000., 1100., 1200., 1300., 1400., 1500.};
*/

const Int_t muon_cos_bin_size = 21;
const double muon_cos_bins[muon_cos_bin_size] = {0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30,
						 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65,
						 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00};
const Int_t muon_deg_bin_size = 7;
const double muon_deg_bins[muon_deg_bin_size] = {0.0, 10., 20., 30., 40., 50., 90.};
const double muon_deg_bin_weight[muon_deg_bin_size-1] = {1., 1., 1., 1., 1., 1/4.};
/*
const double muon_deg_bins[muon_deg_bin_size] = {0.0, 5.0, 10., 15., 20., 25., 30.,
						 35., 40., 45., 50., 55., 60., 65.,
						 70., 75., 80., 85., 90.};
*/

const Int_t hadron_cos_bin_size = 41;
const double hadron_cos_bins[hadron_cos_bin_size] = {-1.00, -0.95, -0.90, -0.85, -0.80, -0.75, -0.70,
						     -0.65, -0.60, -0.55, -0.50, -0.45, -0.40, -0.35,
						     -0.30, -0.25, -0.20, -0.15, -0.10, -0.05,
						     0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30,
						     0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65,
						     0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00};

const Int_t hadron_deg_bin_size = 10;
const double hadron_deg_bins[hadron_deg_bin_size] = {0.0, 20., 40., 60., 80.,
						     100., 120., 140., 160., 180.};
const double hadron_deg_bin_weight[hadron_deg_bin_size-1] = {1., 1., 1., 1., 1., 1., 1., 1., 1.};
/*
const double hadron_deg_bins[hadron_deg_bin_size] = {0.0, 5.0, 10., 15., 20., 25.,
						     30., 35., 40., 45., 50., 55.,
						     60., 65., 70., 75., 80., 85.,
						     90., 95., 100, 105, 110, 115,
						     120, 125, 130, 135, 140, 145,
						     150, 155, 160, 165, 170, 175, 180};
*/

const Int_t q2_bin_size = 11;
const double q2_bins[q2_bin_size] = {0., 1e6, 2e6, 3e6, 4e6, 5e6,
				     6e6, 7e6, 8e6, 9e6, 10e6};
const Int_t nu_ene_recon_bin_size = 11;
const double nu_ene_recon_bins[nu_ene_recon_bin_size] = {0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000};
const Int_t dpt_bin_size = 21;
const double dpt_bins[dpt_bin_size] = {0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 
				       500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000};
const Int_t dalphat_bin_size = 10;
const double dalphat_bins[dalphat_bin_size] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180};
const Int_t cosdat_bin_size = 11;
const double cosdat_bins[cosdat_bin_size] = {-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1};
const Int_t dphit_bin_size = 19;
const double dphit_bins[dphit_bin_size] = {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180};
const Int_t cosdphit_bin_size = 11;
const double cosdphit_bins[cosdphit_bin_size] = {-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1};
const Int_t dptx_bin_size = 21;
const double dptx_bins[dptx_bin_size] = {-1000, -900, -800, -700, -600, -500, -400, -300, -200, -100, 0, 
					 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000};
const Int_t dpty_bin_size = 21;
const double dpty_bins[dpty_bin_size] = {-1000, -900, -800, -700, -600, -500, -400, -300, -200, -100, 0, 
					 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000};

const Int_t open_deg_bin_size = 37;
const double open_deg_bins[open_deg_bin_size] = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45,
						 50, 55, 60, 65, 70, 75, 80, 85, 90,
						 95, 100, 105, 110, 115, 120, 125, 130, 135,
						 140, 145, 150, 155, 160, 165, 170, 175, 180};
const Int_t open_cos_bin_size = 41;
const double open_cos_bins[open_cos_bin_size] = {-1.00, -0.95, -0.90, -0.85, -0.80, -0.75, -0.70, -0.65, -0.60, -0.55, 
						 -0.50, -0.45, -0.40, -0.35, -0.30, -0.25, -0.20, -0.15, -0.10, -0.05,
						 0.00,
						 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50,
						 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00};
const Int_t mom_ratio_bin_size = 11;
const double mom_ratio_bins[mom_ratio_bin_size] = {0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5};
const Int_t mom_vecsum_bin_size = 21;
const double mom_vecsum_bins[mom_vecsum_bin_size] = {0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,
						     1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000};
const Int_t mom_scasum_bin_size = 21;
const double mom_scasum_bins[mom_scasum_bin_size] = {0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,
						     1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000};
const Int_t dptt_bin_size = 21;
const double dptt_bins[dpt_bin_size] = {-500,-450,-400,-350,-300,-250,-200,-150,-100,-50,0,
					50,100,150,200,250,300,350,400,450,500};
const Int_t dpt_2p_bin_size = 21;
const double dpt_2p_bins[dpt_2p_bin_size] = {0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,
					     1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000};
const Int_t pn_bin_size = 21;
const double pn_bins[pn_bin_size] = {0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,
				     1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000};
const Int_t dalphat_2p_bin_size = 37;
const double dalphat_2p_bins[dalphat_2p_bin_size] = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45,
						     50, 55, 60, 65, 70, 75, 80, 85, 90,
						     95, 100, 105, 110, 115, 120, 125, 130, 135,
						     140, 145, 150, 155, 160, 165, 170, 175, 180};
const Int_t cosdat_2p_bin_size = 41;
const double cosdat_2p_bins[cosdat_2p_bin_size] = {-1, -0.95, -0.9, -0.85, -0.8, -0.75, -0.7, -0.65, -0.6, -0.55,
						   -0.5, -0.45, -0.4, -0.35, -0.3, -0.25, -0.2, -0.15, -0.1, -0.05,
						   0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
						   0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.};

#endif
