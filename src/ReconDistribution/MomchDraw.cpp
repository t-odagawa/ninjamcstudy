// system includes
#include <sstream>

// boost includes
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp> 
#include <boost/log/expressions.hpp>

#include "AnalyzeMultiplicity.hpp"
#include "AnalyzeVertex.hpp"
#include "AnalyzeMomentum.hpp"
#include "AnalyzePid.hpp"

namespace logging = boost::log;

int main (int argc, char *argv[]) {

  logging::core::get()->set_filter
    (
     // logging::trivial::severity >= logging::trivial::info
     logging::trivial::severity >= logging::trivial::debug
     );

  if ( argc != 4 ) {
    BOOST_LOG_TRIVIAL(error) << "Usage : " << argv[0]
			     << " <fileid> <ecc> <mode>";
    std::exit(1);
  }

  BOOST_LOG_TRIVIAL(info) << "==========Start==========";

  try {

    int fileid = std::atoi(argv[1]);
    int eccid = std::atoi(argv[2]);
    int modeid = std::atoi(argv[3]);

    // input B2 file
    std::stringstream b2filename_ss;
    b2filename_ss << "/hsm/nu/ninja/pra_tmp/mc_tmp_20220505/ninja_mc_h2o_"
		  << fileid
		  << ".root";
    std::string b2filename = b2filename_ss.str();

    // input momch file
    std::stringstream momchfilename_ss;
    momchfilename_ss << "/hsm/nu/ninja/pra_tmp/mc_tmp_20220505/momch/momch_ecc"
		     << eccid << "_" << fileid;
    if ( modeid == 0 || modeid == 1 ) momchfilename_ss << ".momch";
    else if ( modeid == 2 ) momchfilename_ss << "_addmcs.momch";
    else if ( modeid == 4 ) momchfilename_ss << "_pid.momch";
    std::string momchfilename = momchfilename_ss.str();

    // output file
    std::stringstream outputfilename_ss;
    outputfilename_ss << "/hsm/nu/ninja/pra_tmp/mc_tmp_20220505/output/output_mode"
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
      // AnalyzeAngle(b2filename, momchfilename, outputfilename);
      break;
    case 4 :
      AnalyzePid(b2filename, momchfilename, outputfilename);
      break;
    default :
      throw std::runtime_error("Select one analyze mode : 0:Multiplicity, 1:Vertex, 2:Momentum, 3:Angle, 4:PID");
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
