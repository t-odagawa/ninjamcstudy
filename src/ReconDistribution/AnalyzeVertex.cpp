#include <iostream>
#include <vector>
#include <cmath>
#include <string>

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/filesystem.hpp>

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>

#include <B2Reader.hh>
#include <B2SpillSummary.hh>
#include <B2EventSummary.hh>
#include <B2VertexSummary.hh>

#include <McsClass.hpp>

#include "HistogramStyle.hpp"
#include "AnalyzeVertex.hpp"

namespace logging = boost::log;
namespace fs = boost::filesystem;

void AnalyzeVertex(std::string b2filename,
		   std::string momchfilename,
		   std::string outputfilename) {

  BOOST_LOG_TRIVIAL(info) << "==========Vertex mode==========";
  
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
  
  TH1D *hist_recon_vertex_x = new TH1D("hist_recon_vertex_x", "Reconstructed horizontal position;x [#mum];Entries", 25, 0., 250.e3);
  TH1D *hist_recon_vertex_y = new TH1D("hist_recon_vertex_y", "Reconstructed vertical position;y [#mum];Entries", 25, 0., 250.e3);
  TH1D *hist_recon_vertex_z = new TH1D("hist_recon_vertex_z", "Reconstructed depth;z [#mum];Entries", 25, -250.e3, 0.);
  TH2D *hist_recon_vertex_xy = new TH2D("hist_recon_vertex_xy", "Reconstructed film position; x [#mum];y [#mum]",
					25, 0., 250.e3, 25, 0., 250.e3);
  TH2D *hist_recon_vertex_xz = new TH2D("hist_recon_vertex_xz", "Reconstructed film position; z [#mum];x [#mum]",
					25, -250.e3, 0., 25, 0., 250.e3);
  TH2D *hist_recon_vertex_yz = new TH2D("hist_recon_vertex_yz", "Reconstructed film position; z [#mum];y [#mum]",
					25, -250.e3, 0., 25, 0., 250.e3);
  // Z distributions

  // Muon is correctly id-ed
  TH1D *hist_recon_vertex_pl[num_ninja_mode];
  for ( int i = 0; i < num_ninja_mode; i++ ) {
    hist_recon_vertex_pl[i] = new TH1D(Form("hist_vertex_pl_%d", i), "",
				       133, 0.5, 133.5);
    hist_recon_vertex_pl[i]->SetFillColor(mode_color[i]);
    hist_recon_vertex_pl[i]->SetFillStyle(mode_style[i]);
  }
  // Muon is mis-id-ed
  TH1D *hist_recon_vertex_pl_muon_misid = new TH1D("hist_recon_vertex_pl_muon_misid", "",
						   133, 0.5, 133.5);

  TH1D *hist_resolution_x = new TH1D("hist_resolution_x", "Horizontal positional resolution;#Deltax [#mum];Entries", 100, -500, 500);
  TH1D *hist_resolution_y = new TH1D("hist_resolution_y", "Vertical positional resolution;#Deltay [#mum];Entries", 100, -500, 500);
  TH1D *hist_resolution_z = new TH1D("hist_resolution_z", "Depth resolution;#Deltaz [#mum];Entries", 100, -500, 500);
  
  TH1D *hist_target_material = new TH1D("hist_target_material", "Target material;Material;Entries", 6, -0.5, 5.5);
  hist_target_material->GetXaxis()->SetBinLabel(1, "Water");
  hist_target_material->GetXaxis()->SetBinLabel(2, "Base");
  hist_target_material->GetXaxis()->SetBinLabel(3, "Iron");
  hist_target_material->GetXaxis()->SetBinLabel(4, "Unknown 1");
  hist_target_material->GetXaxis()->SetBinLabel(5, "Unknown 2");
  hist_target_material->GetXaxis()->SetBinLabel(6, "Emulsion");

  // chamber reconstructed success
  TH1D *hist_vertex_pl_diff = new TH1D("hist_vertex_pl_diff", "Vertex plate difference;#Deltapl;Entries", 260, -130.5, 129.5);
  TH2D *hist_vertex_pl_recon_true = new TH2D("hist_vertex_pl_recon_true", "Vertex plate;True plate;Reconstructed plate",
					     130, 0.5, 130.5, 130, 0.5, 130.5);
  
  // chamber reconstructed insuccess  
  TH2D *hist_ecc_recon_true = new TH2D("hist_ecc_recon_true", "Reconstructed ECC;True ECC;Reconstructed ECC",
				       9, 0.5, 9.5, 9, 0.5, 9.5);
  TH2D *hist_recon_vertex_xy_pene = new TH2D("hist_recon_vertex_xy_pene", "Reconstructed film position; x [#mum];y [#mum]",
					25, 0., 250.e3, 25, 0., 250.e3);  
  TH2D *hist_recon_vertex_xy_side = new TH2D("hist_recon_vertex_xy_side", "Reconstructed film position; x [#mum];y [#mum]",
					25, 0., 250.e3, 25, 0., 250.e3);

  for ( auto ev : ev_vec ) {

    reader.ReadSpill(ev.groupid);
    auto &spill_summary = reader.GetSpillSummary();
    auto it_event = spill_summary.BeginTrueEvent();
    const auto *event = it_event.Next();

    auto &vertex = event->GetPrimaryVertex();
    int mode_id = GetNinjaModeId(vertex.GetInteractionType());

    if ( ev.chains.empty() ) continue;

    int true_ecc_id = ev.ecc_id / 10;
    int ecc_id = ev.ecc_id % 10;
    int true_vertex_pl = ev.vertex_pl / 1000;
    int vertex_pl = ev.vertex_pl % 1000;

    bool muon_correct_flag = false;
    for ( auto chain : ev.chains ) {
      if ( chain.particle_flag % 10000 == 13 &&
	   chain.particle_flag / 10000 == 13 )
	muon_correct_flag = true;
    }

    if ( ecc_id == 0 ) continue;

    if ( ev.vertex_material >= 0 ) { // stop identified only
      hist_recon_vertex_x->Fill(ev.recon_vertex_position[0], ev.weight);
      hist_recon_vertex_y->Fill(ev.recon_vertex_position[1], ev.weight);
      hist_recon_vertex_z->Fill(ev.recon_vertex_position[2], ev.weight);
      hist_recon_vertex_xy->Fill(ev.recon_vertex_position[0], ev.recon_vertex_position[1], ev.weight);
      hist_recon_vertex_xz->Fill(ev.recon_vertex_position[2], ev.recon_vertex_position[0], ev.weight);
      hist_recon_vertex_yz->Fill(ev.recon_vertex_position[2], ev.recon_vertex_position[1], ev.weight);

      hist_resolution_x->Fill(ev.recon_vertex_position[0] - ev.true_vertex_position[0],
			      ev.weight);
      hist_resolution_y->Fill(ev.recon_vertex_position[1] - ev.true_vertex_position[1],
			      ev.weight);
      hist_resolution_z->Fill(ev.recon_vertex_position[2] - ev.true_vertex_position[2],
			      ev.weight);
      
      if ( muon_correct_flag )
	hist_recon_vertex_pl[mode_id]->Fill(vertex_pl, ev.weight);
      else 
	hist_recon_vertex_pl_muon_misid->Fill(vertex_pl, ev.weight);

      hist_target_material->Fill(ev.vertex_material, ev.weight);
    }

    if ( true_ecc_id == ecc_id ) {
      hist_vertex_pl_diff->Fill(vertex_pl - true_vertex_pl, ev.weight);
      hist_vertex_pl_recon_true->Fill(true_vertex_pl, vertex_pl, ev.weight);
    }
    else {
      hist_ecc_recon_true->Fill(true_ecc_id, ecc_id, ev.weight);
      if ( ev.vertex_pl > 131 )
	hist_recon_vertex_xy_pene->Fill(ev.recon_vertex_position[0], ev.recon_vertex_position[1], ev.weight);
      else
	hist_recon_vertex_xy_side->Fill(ev.recon_vertex_position[0], ev.recon_vertex_position[1], ev.weight);      
    }


  }

  outputfile->cd();
  hist_recon_vertex_x->Write();
  hist_recon_vertex_y->Write();
  hist_recon_vertex_z->Write();
  hist_recon_vertex_xy->Write();
  hist_recon_vertex_xz->Write();
  hist_recon_vertex_yz->Write();
  hist_resolution_x->Write();
  hist_resolution_y->Write();
  hist_resolution_z->Write();
  for ( int i = 0; i < num_ninja_mode; i++ )
    hist_recon_vertex_pl[i]->Write();
  hist_recon_vertex_pl_muon_misid->Write();
  hist_target_material->Write();
  hist_vertex_pl_diff->Write();
  hist_vertex_pl_recon_true->Write();
  hist_ecc_recon_true->Write();
  hist_recon_vertex_xy_pene->Write();
  hist_recon_vertex_xy_side->Write();
  outputfile->Close();
  
}
