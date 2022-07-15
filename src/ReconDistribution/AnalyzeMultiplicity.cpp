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
				    10, 0.5, 10.5);
  TH2D *hist_total_multi_recon_true = new TH2D("hist_total_multi_recon_true", "Multiplicity;True multiplicity;Reconstructed multiplicity",
					       10, 0.5, 10.5, 10, 0.5, 10.5);

  TH1D *hist_water_total_multi = new TH1D("hist_water_total_multi", "Total multiplicity (Water);# of tracks;Entries",
					  10, 0.5, 10.5);
  TH1D *hist_iron_total_multi = new TH1D("hist_iron_total_multi", "Total multiplicity (Iron);# of tracks;Entries",
					 10, 0.5, 10.5);

  TH1D *hist_water_mode_multi[num_ninja_mode];
  TH1D *hist_iron_mode_multi[num_ninja_mode];
  TH1D *hist_water_mode_proton_multi[num_ninja_mode];
  TH1D *hist_water_mode_pion_multi[num_ninja_mode];
  for ( int i = 0; i < num_ninja_mode; i++ ) {
    hist_water_mode_multi[i] = new TH1D(Form("hist_water_mode_multi_%d", i), "", 10, 0.5, 10.5);
    hist_iron_mode_multi[i] = new TH1D(Form("hist_iron_mode_multi_%d", i), "", 10, 0.5, 10.5);
    hist_water_mode_proton_multi[i] = new TH1D(Form("hist_water_mode_proton_multi_%d", i),"", 10, 0.5, 10.5);
    hist_water_mode_pion_multi[i] = new TH1D(Form("hist_water_mode_pion_multi_%d", i),"", 10, 0.5, 10.5);
    hist_water_mode_multi[i]->SetFillColor(mode_color[i]);
    hist_iron_mode_multi[i]->SetFillColor(mode_color[i]);
    hist_water_mode_proton_multi[i]->SetFillColor(mode_color[i]);
    hist_water_mode_pion_multi[i]->SetFillColor(mode_color[i]);
    hist_water_mode_multi[i]->SetFillStyle(mode_style[i]);
    hist_iron_mode_multi[i]->SetFillStyle(mode_style[i]);
    hist_water_mode_proton_multi[i]->SetFillStyle(mode_style[i]);
    hist_water_mode_pion_multi[i]->SetFillStyle(mode_style[i]);
  }

  TH1D *hist_water_proton_misid_multi = new TH1D("hist_water_proton_misid_multi",
						 "Proton multiplicity (Water);# of tracks;Entries",
						 10, 0.5, 10.5);
  TH1D *hist_water_pion_misid_multi = new TH1D("hist_water_pion_misid_multi",
					       "Pion multiplicity (Water);# of tracks;Entries",
					       10, 0.5, 10.5);
  // charge check
  TH1D *hist_recon_charge = new TH1D("hist_recon_charge", "Charge ID;Reconstructed charge sign;Entries", 3, -1.5, 1.5);

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
	hist_water_mode_multi[mode_id]->Fill(ev.chains.size(), ev.weight);
	int num_proton_water = 0;
	int num_pion_water = 0;
	for ( auto chain : ev.chains ) {
	  if ( chain.particle_flag % 10000 == 13 ) hist_recon_charge->Fill(chain.charge_sign, ev.weight);
	  if ( chain.particle_flag % 10000 == 2212 ) num_proton_water++;
	  else if ( chain.particle_flag % 10000 == 211 ) num_pion_water++;	
	}
	hist_water_mode_proton_multi[mode_id]->Fill(num_proton_water, ev.weight);
	hist_water_mode_pion_multi[mode_id]->Fill(num_pion_water, ev.weight);
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
  hist_iron_total_multi->Write();
  for ( int i = 0; i < num_ninja_mode; i++ ) {
    hist_water_mode_multi[i]->Write();
    hist_iron_mode_multi[i]->Write();
    hist_water_mode_proton_multi[i]->Write();
    hist_water_mode_pion_multi[i]->Write();
  }
  hist_water_proton_misid_multi->Write();
  hist_water_pion_misid_multi->Write();
  hist_recon_charge->Write();
  outputfile->Close();
  
}
