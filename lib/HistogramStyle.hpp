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

const Int_t total_multi_bin_size = 11;
const double total_multi_bins[total_multi_bin_size] = {0.5, 1.5, 2.5, 3.5, 4.5,
						       5.5, 6.5, 7.5, 8.5, 9.5, 10.5};

const Int_t hadron_multi_bin_size = 11;
const double hadron_multi_bins[hadron_multi_bin_size] = {-0.5, 0.5, 1.5, 2.5, 3.5,
							 4.5, 5.5, 6.5, 7.5, 8.5, 9.5};

const Int_t muon_mom_bin_size = 21;
const double muon_mom_bins[muon_mom_bin_size] = {0., 100., 200., 300., 400., 500., 600.,
						 700., 800., 900., 1000., 1100., 1200., 1300.,
						 1400., 1500., 1600., 1700., 1800., 1900., 2000.};

const Int_t hadron_mom_bin_size = 16;
const double hadron_mom_bins[hadron_mom_bin_size] = {0., 100., 200., 300., 400., 500., 600., 700.,
						     800., 900., 1000., 1100., 1200., 1300., 1400., 1500.};

const Int_t muon_cos_bin_size = 21;
const double muon_cos_bins[muon_cos_bin_size] = {0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30,
						 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65,
						 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00};
const Int_t muon_deg_bin_size = 19;
const double muon_deg_bins[muon_deg_bin_size] = {0.0, 5.0, 10., 15., 20., 25., 30.,
						 35., 40., 45., 50., 55., 60., 65.,
						 70., 75., 80., 85., 90.};

const Int_t hadron_cos_bin_size = 41;
const double hadron_cos_bins[hadron_cos_bin_size] = {-1.00, -0.95, -0.90, -0.85, -0.80, -0.75, -0.70,
						     -0.65, -0.60, -0.55, -0.50, -0.45, -0.40, -0.35,
						     -0.30, -0.25, -0.20, -0.15, -0.10, -0.05,
						     0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30,
						     0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65,
						     0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00};

const Int_t hadron_deg_bin_size = 37;
const double hadron_deg_bins[hadron_deg_bin_size] = {0.0, 5.0, 10., 15., 20., 25.,
						     30., 35., 40., 45., 50., 55.,
						     60., 65., 70., 75., 80., 85.,
						     90., 95., 100, 105, 110, 115,
						     120, 125, 130, 135, 140, 145,
						     150, 155, 160, 165, 170, 175, 180};

#endif
