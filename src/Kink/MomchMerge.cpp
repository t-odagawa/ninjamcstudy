#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <McsClass.hpp>

#include "MomchMerge.hpp"

int main (int argc, char* argv[]) {

  if ( argc != 3 ) {
    std::cerr << "Usage : " << argv[0]
	      << " <input file name> <output file name>" << std::endl;
    std::exit(1);
  }
  std::string input_dir_name = argv[1];
  std::string output_file_name = argv[2];
  std::vector<Momentum_recon::Event_information> ev_vec;

  for ( int i = 1; i < 1000; i++ ) {
    ReadMergeMomch(input_dir_name, i, ev_vec);
  }

  Momentum_recon::WriteEventInformationBin(output_file_name, ev_vec);

  std::cout << "Finish" << std::endl;
  std::exit(0);

}

void ReadMergeMomch(std::string input_dir_name, int i, std::vector<Momentum_recon::Event_information> &ev_vec) {

  std::stringstream ss_filename;
  ss_filename << input_dir_name << i << ".momch";
  // std::cout << ss_filename.str() << std::endl;

  std::vector<Momentum_recon::Event_information> ev_vec_i = Momentum_recon::ReadEventInformationBin(ss_filename.str());
  if ( ev_vec_i.empty() ) return;

  for ( auto &ev : ev_vec_i ) {
    ev.groupid += i * 10000;
  }

  ev_vec.insert(ev_vec.end(), ev_vec_i.begin(), ev_vec_i.end());

  return;
}
