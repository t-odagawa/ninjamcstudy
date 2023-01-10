#include <iostream>
#include <vector>
#include <string>

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/filesystem.hpp>

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>

#include <B2Enum.hh>
#include <B2Reader.hh>
#include <B2SpillSummary.hh>
#include <B2EventSummary.hh>
#include <B2VertexSummary.hh>

#include <McsClass.hpp>

#include "HistogramStyle.hpp"
#include "AnalyzeMultiplicity.hpp"

namespace logging = boost::log;
namespace fs = boost::filesystem;

void AnalyzeMultiplicity(std::string b2filename,
			 std::string momchfilename,
			 std::string outputfilename) {

  BOOST_LOG_TRIVIAL(info) << "==========Multiplicity mode==========";
  
  // input B2 file
  B2Reader reader(b2filename);
  
  // input momch file
  if ( !fs::exists(momchfilename) ) {
    throw std::runtime_error("File not found : " + momchfilename);
  }
  std::vector<Momentum_recon::Event_information > ev_vec = Momentum_recon::ReadEventInformationBin(momchfilename);
  
  // output file
  TFile *outputfile = new TFile((TString)outputfilename, "recreate");
  BOOST_LOG_TRIVIAL(info) << "Output filename : " << outputfilename;

  TH1D *hist_total_multi = new TH1D("hist_total_multi", "Total multiplicity;# of tracks;Entries",
				    total_multi_bin_size - 1, total_multi_bins);
  TH2D *hist_total_multi_recon_true = new TH2D("hist_total_multi_recon_true", "Multiplicity;True multiplicity;Reconstructed multiplicity",
					       total_multi_bin_size - 1, total_multi_bins,
					       total_multi_bin_size - 1, total_multi_bins);

  TH1D *hist_water_total_multi = new TH1D("hist_water_total_multi", "Total multiplicity (Water);# of tracks;Entries",
					  total_multi_bin_size - 1, total_multi_bins);
  TH1D *hist_water_proton_multi = new TH1D("hist_water_proton_multi", "Proton multiplicity (Water);# of protons;Entries",
					   hadron_multi_bin_size - 1, hadron_multi_bins);
  TH1D *hist_water_pion_multi = new TH1D("hist_water_pion_multi", "Pion multiplicity (Water);# of pions;Entries",
					 hadron_multi_bin_size - 1, hadron_multi_bins);
  TH1D *hist_iron_total_multi = new TH1D("hist_iron_total_multi", "Total multiplicity (Iron);# of tracks;Entries",
					 total_multi_bin_size - 1, total_multi_bins);


  // Muon is correctly id-ed
  TH1D *hist_water_mode_multi[num_ninja_mode];
  TH1D *hist_iron_mode_multi[num_ninja_mode];
  TH1D *hist_water_mode_proton_multi[num_ninja_mode];
  TH1D *hist_water_mode_pion_multi[num_ninja_mode];
  TH1D *hist_water_mode_other_multi[num_ninja_mode];
  for ( int i = 0; i < num_ninja_mode; i++ ) {
    hist_water_mode_multi[i] = new TH1D(Form("hist_water_mode_multi_%d", i), "", total_multi_bin_size - 1, total_multi_bins);
    hist_iron_mode_multi[i] = new TH1D(Form("hist_iron_mode_multi_%d", i), "", total_multi_bin_size - 1, total_multi_bins);
    hist_water_mode_proton_multi[i] = new TH1D(Form("hist_water_mode_proton_multi_%d", i),"",
					       hadron_multi_bin_size - 1, hadron_multi_bins);
    hist_water_mode_pion_multi[i] = new TH1D(Form("hist_water_mode_pion_multi_%d", i),"",
					     hadron_multi_bin_size - 1, hadron_multi_bins);
    hist_water_mode_other_multi[i] = new TH1D(Form("hist_water_mode_other_multi_%d", i),"",
					      hadron_multi_bin_size - 1, hadron_multi_bins);
    hist_water_mode_multi[i]->SetFillColor(mode_color[i]);
    hist_iron_mode_multi[i]->SetFillColor(mode_color[i]);
    hist_water_mode_proton_multi[i]->SetFillColor(mode_color[i]);
    hist_water_mode_pion_multi[i]->SetFillColor(mode_color[i]);
    hist_water_mode_other_multi[i]->SetFillColor(mode_color[i]);
    hist_water_mode_multi[i]->SetFillStyle(mode_style[i]);
    hist_iron_mode_multi[i]->SetFillStyle(mode_style[i]);
    hist_water_mode_proton_multi[i]->SetFillStyle(mode_style[i]);
    hist_water_mode_pion_multi[i]->SetFillStyle(mode_style[i]);
    hist_water_mode_other_multi[i]->SetFillStyle(mode_style[i]);
  }

  // Packing
  TH1D *hist_single_proton_multi = new TH1D("hist_single_proton_multi", "",
					    hadron_multi_bin_size - 1, hadron_multi_bins);
  TH1D *hist_single_pion_multi = new TH1D("hist_single_pion_multi", "", hadron_multi_bin_size - 1, hadron_multi_bins);

  // Muon is mis-id-ed
  TH1D *hist_water_multi_muon_misid = new TH1D("hist_water_multi_muon_misid",
					       "", total_multi_bin_size - 1, total_multi_bins);
  TH1D *hist_water_proton_multi_muon_misid = new TH1D("hist_water_proton_multi_muon_misid",
						      "", hadron_multi_bin_size - 1, hadron_multi_bins);
  TH1D *hist_water_pion_multi_muon_misid = new TH1D("hist_water_pion_multi_muon_misid",
						    "", hadron_multi_bin_size - 1, hadron_multi_bins);
  TH1D *hist_water_other_multi_muon_misid = new TH1D("hist_water_other_multi_muon_misid",
						     "", hadron_multi_bin_size - 1, hadron_multi_bins);

  // Muon is correctly id-ed but partner is not
  TH1D *hist_water_proton_misid_multi = new TH1D("hist_water_proton_misid_multi",
						 "Proton multiplicity (Water);# of tracks;Entries",
						 hadron_multi_bin_size - 1, hadron_multi_bins);
  TH1D *hist_water_pion_misid_multi = new TH1D("hist_water_pion_misid_multi",
					       "Pion multiplicity (Water);# of tracks;Entries",
					       hadron_multi_bin_size - 1, hadron_multi_bins);
  TH1D *hist_water_other_misid_multi = new TH1D("hist_water_other_misid_multi",
						"Non-ided hadron multiplicity (Water);# of tracks;Entries",
						hadron_multi_bin_size - 1, hadron_multi_bins);
  // charge check
  TH1D *hist_recon_charge = new TH1D("hist_recon_charge", "Charge ID;Reconstructed charge sign;Entries", 3, -1.5, 1.5);


  // For flux systematic uncertainty study
  // mis-pid もそれぞれいるかも？
  TH2D *hist_flux = new TH2D("hist_flux", ";# of tracks;E_{#nu} [GeV]",
			     total_multi_bin_size - 1, total_multi_bins, nu_ene_bin_size - 1, nu_ene_bins);
  TH2D *hist_flux_p = new TH2D("hist_flux_p", ";# of protons;E_{#nu} [GeV]",
			       hadron_multi_bin_size - 1, hadron_multi_bins, nu_ene_bin_size - 1, nu_ene_bins);
  TH2D *hist_flux_pi = new TH2D("hist_flux_pi", ";# of pions;E_{#nu} [GeV]",
				hadron_multi_bin_size - 1, hadron_multi_bins, nu_ene_bin_size - 1, nu_ene_bins);
  

  for ( auto ev : ev_vec ) {

    reader.ReadSpill(ev.groupid);
    auto &spill_summary = reader.GetSpillSummary();
    auto it_event = spill_summary.BeginTrueEvent();
    const auto *event = it_event.Next();

    auto &vertex = event->GetPrimaryVertex();
    int mode_id = GetNinjaModeId(vertex.GetInteractionType());

    if ( !ev.chains.empty() ) {
      hist_total_multi->Fill(ev.chains.size(), ev.weight);
      hist_total_multi_recon_true->Fill(ev.true_chains.size(), ev.chains.size(), ev.weight);
      if ( ev.vertex_material == B2Material::kWater ) {
	hist_water_total_multi->Fill(ev.chains.size(), ev.weight);
	bool muon_correct_id_flag = false;
	int num_proton_water = 0;
	int num_true_proton_water = 0;
	int num_pion_water = 0;
	int num_true_pion_water = 0;
	for ( auto chain : ev.chains ) {
	  if ( chain.particle_flag % 10000 == 13 ) {
	    if ( chain.particle_flag / 10000 == 13 ) muon_correct_id_flag = true;
	    hist_recon_charge->Fill(chain.charge_sign, ev.weight);
	  }
	  if ( chain.particle_flag % 10000 == 2212 ) {
	    if ( chain.particle_flag / 10000 == 2212 ) num_true_proton_water++;
	    num_proton_water++;
	  }
	  else if ( chain.particle_flag % 10000 == 211 ) {
	    if ( chain.particle_flag / 10000 == 211 ) num_true_pion_water++;
	    num_pion_water++;	
	  }
	}

	if ( ev.chains.size() != num_proton_water + 1 ) continue;       	

	if ( ev.chains.size() == 1 ) {
	  hist_single_proton_multi->Fill(0., ev.weight);
	  hist_single_pion_multi->Fill(0., ev.weight);
	}

	if ( muon_correct_id_flag ) {
	  hist_water_mode_multi[mode_id]->Fill(ev.chains.size(), ev.weight);
	  hist_flux->Fill(ev.chains.size(), ev.nu_energy / 1000., ev.weight);

	  if ( num_proton_water == num_true_proton_water &&
	       num_proton_water + num_pion_water + 1 == ev.chains.size() ) {
	    hist_water_proton_multi->Fill((double)num_proton_water, ev.weight);
	    hist_water_mode_proton_multi[mode_id]->Fill((double)num_proton_water, ev.weight);
	    hist_flux_p->Fill((double)num_proton_water, ev.nu_energy / 1000., ev.weight);
	  }
	  else 
	    hist_water_proton_misid_multi->Fill((double)num_proton_water, ev.weight);
	  
	  if ( num_pion_water == num_true_pion_water &&
	       num_proton_water + num_pion_water + 1 == ev.chains.size() ) {
	    hist_water_pion_multi->Fill((double)num_pion_water, ev.weight);
	    hist_water_mode_pion_multi[mode_id]->Fill((double)num_pion_water, ev.weight);
	    hist_flux_pi->Fill((double)num_pion_water, ev.nu_energy / 1000., ev.weight);
	  }
	  else
	    hist_water_pion_misid_multi->Fill((double)num_pion_water, ev.weight);

	  if ( num_proton_water + num_pion_water + 1 == ev.chains.size() ) {
	    hist_water_mode_other_multi[mode_id]->Fill(ev.chains.size() - num_proton_water
						       - num_pion_water - 1, ev.weight);
	  }
	  else {
	    hist_water_other_misid_multi->Fill(ev.chains.size() - num_proton_water
					       - num_pion_water - 1, ev.weight);
	  }

	}
	else {
	  std::cout << "Muon mis id : " << ev.groupid << std::endl;
	  hist_water_multi_muon_misid->Fill(ev.chains.size(), ev.weight);
	  hist_water_proton_multi_muon_misid->Fill((double)num_proton_water, ev.weight);
	  hist_water_pion_multi_muon_misid->Fill((double)num_pion_water, ev.weight);
	  hist_water_other_multi_muon_misid->Fill(ev.chains.size() - num_proton_water - num_pion_water - 1, ev.weight);
	}
      }
      else if ( ev.vertex_material == B2Material::kIron ) {
	hist_iron_total_multi->Fill(ev.chains.size(), ev.weight);
	hist_iron_mode_multi[mode_id]->Fill(ev.chains.size(), ev.weight);
      }	
    }
    
  }

  outputfile->cd();
  hist_total_multi->Write();
  hist_total_multi_recon_true->Write();
  hist_water_total_multi->Write();
  hist_water_proton_multi->Write();
  hist_water_pion_multi->Write();
  hist_iron_total_multi->Write();
  for ( int i = 0; i < num_ninja_mode; i++ ) {
    hist_water_mode_multi[i]->Write();
    hist_iron_mode_multi[i]->Write();
    hist_water_mode_proton_multi[i]->Write();
    hist_water_mode_pion_multi[i]->Write();
    hist_water_mode_other_multi[i]->Write();
  }
  hist_water_multi_muon_misid->Write();
  hist_water_proton_multi_muon_misid->Write();
  hist_water_pion_multi_muon_misid->Write();
  hist_water_other_multi_muon_misid->Write();
  hist_water_proton_misid_multi->Write();
  hist_water_pion_misid_multi->Write();
  hist_water_other_misid_multi->Write();

  hist_single_proton_multi->Write();
  hist_single_pion_multi->Write();

  hist_recon_charge->Write();

  hist_flux->Write();
  hist_flux_p->Write();
  hist_flux_pi->Write();

  outputfile->Close();
  
}
