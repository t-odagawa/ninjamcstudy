#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

#include <B2Reader.hh>
#include <B2SpillSummary.hh>
#include <B2EventSummary.hh>
#include <B2VertexSummary.hh>
#include <B2TrackSummary.hh>
#include <B2Pdg.hh>

#include <TFile.h>
#include <TH1D.h>
#include <TVector3.h>

namespace logging = boost::log;

int main (int argc, char* argv[]) {

  logging::core::get()->set_filter
    (
     logging::trivial::severity >= logging::trivial::info
     );
  
  if ( argc != 3 ) {
    BOOST_LOG_TRIVIAL(error) << "Usage : " << argv[0]
			     << " <input B2 file path> <outpu data file path>";
    std::exit(1);
  }

  BOOST_LOG_TRIVIAL(info) << "==========Start==========";

  try {

    B2Reader reader(argv[1]);
    
    TFile *ofile = new TFile(argv[2], "recreate");
    
    TH1D *kink_momentum_hist_ = new TH1D("kink_momentum_hist_","kink_momentum_hist_",50,0,5000);
    TH1D *kink_ax_hist_ = new TH1D("kink_ax_hist_","kink_ax_hist_",100,-5,5);
    TH1D *kink_ay_hist_ = new TH1D("kink_ay_hist_","kink_ay_hist_",100,-5,5);
    
    while ( reader.ReadNextSpill() > 0 ) {
      
      auto &spill_summary = reader.GetSpillSummary();
      
      auto it_event = spill_summary.BeginTrueEvent();
      const auto *event = it_event.Next();
      Double_t normalization = event->GetNormalization();

      auto &primary_vertex_summary = event->GetPrimaryVertex();
      
      const Int_t num_outgoing_tracks = primary_vertex_summary.GetNumOutgoingTracks();
      if ( num_outgoing_tracks < 1 ) continue;

      const Double_t total_cross_section = primary_vertex_summary.GetTotalCrossSection();
      const Double_t weight = normalization * total_cross_section * 1.e-38 * 6.02e23;

      auto it_track = primary_vertex_summary.BeginTrack();
      while ( const auto *track = it_track.Next() ) {       
	if ( !B2Pdg::IsMuonPlusOrMinus(track->GetParticlePdg()) ) continue;
	Double_t absolute_momentum = track->GetInitialAbsoluteMomentum().GetValue();
	TVector3 direction = track->GetInitialDirection().GetValue();
	
	kink_momentum_hist_->Fill(absolute_momentum, weight);
	kink_ax_hist_->Fill(direction.X() / direction.Z(), weight);
	kink_ay_hist_->Fill(direction.Y() / direction.Z(), weight);
      }


  }

  ofile->cd();
  kink_momentum_hist_->Write();
  kink_ax_hist_->Write();
  kink_ay_hist_->Write();
  ofile->Close();

  } catch (const std::runtime_error &error) {
    BOOST_LOG_TRIVIAL(fatal) << "Runtime error : " << error.what();
    std::exit(1);
  } catch (const std::invalid_argument &error) {
    BOOST_LOG_TRIVIAL(fatal) << "Invalid argument error : " << error.what();
    std::exit(1);
  }
  
  BOOST_LOG_TRIVIAL(info) << "==========Finish==========";
  std::exit(0);
  
}
