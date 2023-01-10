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
			     << " <b2filename> <momchfilename> <ofilenname> <mode id>";
    std::exit(1);
  }

  BOOST_LOG_TRIVIAL(info) << "==========Start==========";

  try {

    // input B2 file
    std::string b2filename = std::string(argv[1]);
    // input momch file
    std::string momchfilename = std::string(argv[2]);
    // output file
    std::string outputfilename = std::string(argv[3]);
    int modeid = std::atoi(argv[4]);

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
