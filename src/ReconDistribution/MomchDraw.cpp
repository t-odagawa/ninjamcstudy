// system includes
#include <sstream>

// boost includes
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp> 
#include <boost/log/expressions.hpp>

#include "AnalyzeMultiplicity.hpp"
#include "AnalyzeVertex.hpp"
#include "AnalyzeMomentum.hpp"
#include "AnalyzeAngle.hpp"
#include "AnalyzeEfficiency.hpp"
#include "AnalyzePid.hpp"
#include "Analyze0pi1p.hpp"
#include "Analyze0pi2p.hpp"
#include "AnalyzePartnerEfficiency.hpp"

namespace logging = boost::log;

int main (int argc, char *argv[]) {

  logging::core::get()->set_filter
    (
     logging::trivial::severity >= logging::trivial::info
     //logging::trivial::severity >= logging::trivial::debug
     );

  if ( argc != 5 ) {
    BOOST_LOG_TRIVIAL(error) << "Usage : " << argv[0]
			     << " <file prefix> <fileid> <ecc> <mode>";
    std::exit(1);
  }

  BOOST_LOG_TRIVIAL(info) << "==========Start==========";

  try {

    std::string prefix = argv[1];

    int fileid = std::atoi(argv[2]);
    int eccid = std::atoi(argv[3]);
    int modeid = std::atoi(argv[4]);

    // input B2 file   
    std::stringstream b2filename_ss;
    b2filename_ss << prefix << "/ninja_mc_h2o_"
		  << fileid
		  << ".root";
    std::string b2filename = b2filename_ss.str();

    // input momch file
    std::stringstream momchfilename_ss;
    momchfilename_ss << prefix << "/momch/momch_ecc"
		     << eccid << "_" << fileid;
    if ( modeid == 1 || modeid == 5 ) momchfilename_ss << "_addbm.momch";
    else if ( modeid == 0 || modeid == 2 || modeid == 3 || modeid == 4 ||
	      modeid == 6 || modeid == 7 || modeid == 8 ) momchfilename_ss << "_pid.momch";
    std::string momchfilename = momchfilename_ss.str();

    // output file
    std::stringstream outputfilename_ss;
    outputfilename_ss << prefix << "/output/output_mode"
		      << modeid << "_" << fileid << ".root";
    std::string outputfilename = outputfilename_ss.str();

    switch (modeid) {
    case 0 :
      AnalyzeMultiplicity(b2filename, momchfilename, outputfilename);
      break;
    case 1 :
      AnalyzeVertex(b2filename, momchfilename, outputfilename);
      break;
    case 2 :
      AnalyzeMomentum(b2filename, momchfilename, outputfilename);
      break;
    case 3 :
      AnalyzeAngle(b2filename, momchfilename, outputfilename);
      break;
    case 4 :
      AnalyzePid(b2filename, momchfilename, outputfilename);
      break;
    case 5 :
      AnalyzeEfficiency(b2filename, momchfilename, outputfilename);
      break;
    case 6 :
      Analyze0pi1p(b2filename, momchfilename, outputfilename);
      break;
    case 7 :
      Analyze0pi2p(b2filename, momchfilename, outputfilename);
      break;
    case 8 :
      AnalyzePartnerEfficiency(b2filename, momchfilename, outputfilename);
      break;
    default :
      throw std::runtime_error("Select one analyze mode : 0:Multiplicity, 1:Vertex, 2:Momentum, 3:Angle, 4:PID, 5:Efficiency, 6:0pi1p, 7:0pi2p, \n8:Partner efficiency");
    }

  } catch (const std::runtime_error &error) {
    BOOST_LOG_TRIVIAL(fatal) << "Runtime error : "<< error.what();
    std::exit(1);
  } catch (const std::invalid_argument &error) {
    BOOST_LOG_TRIVIAL(fatal) << "Invalid argument error : " << error.what();
    std::exit(1);
  }

  BOOST_LOG_TRIVIAL(info) << "==========Finish==========";
  std::exit(0);

}
