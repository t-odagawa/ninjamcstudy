#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

#include <B2Const.hh>
#include <B2VertexSummary.hh>

#include <TRandom3.h>
#include <TVector3.h>

#include "MyFunction.hpp"

Int_t GetInteractionEcc(const B2VertexSummary &primary_vertex_summary) {
  if ( primary_vertex_summary.GetDetector() != B2Detector::kNinja )
    return -1;
  TVector3 vertex_position = primary_vertex_summary.GetAbsolutePosition().GetValue();
  vertex_position.SetX(vertex_position.X() - NINJA_POS_X - NINJA_ECC_POS_X);
  vertex_position.SetX(vertex_position.Y() - NINJA_POS_Y - NINJA_ECC_POS_Y);

  Int_t top_id, side_id;

  if ( std::fabs(vertex_position.X() + NINJA_ECC_GAP_X) < 0.5 * NINJA_DESIC_WIDTH )
    top_id = 0;
  else if ( std::fabs(vertex_position.X()) < 0.5 * NINJA_DESIC_WIDTH )
    top_id = 1;
  else if ( std::fabs(vertex_position.X() - NINJA_ECC_GAP_X) < 0.5 * NINJA_DESIC_WIDTH )
    top_id = 2;
  else return -1;

  if ( std::fabs(vertex_position.Y() - NINJA_ECC_GAP_Y) < 0.5 * NINJA_DESIC_HEIGHT )
    side_id = 0;
  else if ( std::fabs(vertex_position.Y()) < 0.5 * NINJA_DESIC_HEIGHT )
    side_id = 1;
  else if ( std::fabs(vertex_position.Y() + NINJA_ECC_GAP_Y) < 0.5 * NINJA_DESIC_HEIGHT )
    side_id = 2;
  else return -1;

  return 3 * side_id + top_id;

}

Int_t GetVertexPlate(Double_t z_pos /*um*/) { // z_pos is in ECC coordinate

  z_pos /= 1.e3; // um -> mm
  if ( z_pos > - NINJA_EMULSION_LAYER_THICK
       - 14 * NINJA_FILM_THICK
       - 11 * NINJA_IRON_LAYER_THICK
       - NINJA_SS_AC_THICK
       - NINJA_ENV_THICK ) { // iron ECC
    BOOST_LOG_TRIVIAL(debug) << "Iron ECC interaction";
    z_pos = z_pos
      + NINJA_EMULSION_LAYER_THICK 
      + 3 * NINJA_FILM_THICK
      + NINJA_SS_AC_THICK; // iron most downstream position -> origin
    Int_t film_id = (Int_t)(-z_pos / (NINJA_FILM_THICK + NINJA_IRON_LAYER_THICK));
    z_pos = z_pos
      + film_id * (NINJA_FILM_THICK + NINJA_IRON_LAYER_THICK);
    if (- NINJA_IRON_LAYER_THICK < z_pos &&
	z_pos <= 0.) {
      BOOST_LOG_TRIVIAL(debug) << "Iron interaction in iron ECC";
      return film_id + 3;
    }
    else {
      BOOST_LOG_TRIVIAL(warning) << "Non iron interaction in iron ECC?";
      return -1;
    }
  }
  else if ( z_pos > - NINJA_EMULSION_LAYER_THICK
	    - 132 * NINJA_FILM_THICK
	    - 58 * NINJA_WATER_LAYER_THICK
	    - (59 * 2 + 1) * NINJA_ENV_THICK
	    - 70 * NINJA_IRON_LAYER_THICK
	    - NINJA_SS_AC_THICK) { // water ECC
    BOOST_LOG_TRIVIAL(debug) << "Water ECC interaction";
    z_pos = z_pos
      + NINJA_EMULSION_LAYER_THICK
      + 15 * NINJA_FILM_THICK
      + 11 * NINJA_IRON_LAYER_THICK
      + NINJA_SS_AC_THICK
      + 2 * NINJA_ENV_THICK; // iron most downstream in water ECC -> origin
    Int_t unit_id = (Int_t)(-z_pos / (2 * NINJA_FILM_THICK + NINJA_IRON_LAYER_THICK
				      + NINJA_WATER_LAYER_THICK + 2 * NINJA_ENV_THICK));
    z_pos = z_pos
      + unit_id * (2 * NINJA_FILM_THICK + NINJA_IRON_LAYER_THICK
		   + NINJA_WATER_LAYER_THICK + 2 * NINJA_ENV_THICK); // iron most downstream in one unit -> origin
    if (- NINJA_IRON_LAYER_THICK < z_pos &&
	z_pos <= 0.) { // iron interaction
      BOOST_LOG_TRIVIAL(debug) << "Iron interaction in water ECC";
      return 2 * (unit_id + 8) - 1;
    }
    else if (- NINJA_IRON_LAYER_THICK - NINJA_FILM_THICK - NINJA_WATER_LAYER_THICK - NINJA_ENV_THICK < z_pos &&
	     z_pos <= - NINJA_IRON_LAYER_THICK - NINJA_FILM_THICK - NINJA_ENV_THICK) { // water interaction
      BOOST_LOG_TRIVIAL(debug) << "Water interaction in water ECC";
      return 2 * (unit_id + 8);
    }
  }
  else {
    BOOST_LOG_TRIVIAL(warning) << "Non recognized material interaction in water ECC?";
    return -1;
  }

}

void CalcPosInEccCoordinate(TVector3 &position, Int_t ecc_id, Bool_t absolute_flag) {
  // center of ECC5 dessicator
  position.SetX(position.X() - NINJA_POS_X - NINJA_ECC_POS_X);
  position.SetY(position.Y() - NINJA_POS_Y - NINJA_ECC_POS_Y);
  position.SetZ(position.Z() - NINJA_POS_Z - NINJA_ECC_POS_Z);

  // film coordinate
  position.SetX(position.X()
		+ 0.5 * NINJA_ECC_FILM_XY);
  position.SetY(position.Y()
		+ 0.5 * NINJA_DESIC_HEIGHT
		- NINJA_DESIC_THICK
		- NINJA_ENV_THICK);
  position.SetZ(position.Z()
		- 0.5 * NINJA_DESIC_DEPTH
		+ NINJA_DESIC_THICK
		+ NINJA_ENV_THICK
		+ NINJA_EMULSION_LAYER_THICK
		+ NINJA_BASE_LAYER_THICK);

  // move to each ECC
  position.SetX(position.X() 
		+ NINJA_ECC_GAP_X * (1 - ecc_id % 3));
  position.SetY(position.Y()
		+ NINJA_ECC_GAP_Y * (ecc_id / 3 - 1));


  // mm -> um
  position.SetX(position.X() * 1.e3);
  position.SetY(position.Y() * 1.e3);
  position.SetZ(position.Z() * 1.e3);
}

void SmearPosition(TVector3 &position /*um*/) {
  position.SetX(gRandom->Gaus(position.X(), 0.2));
  position.SetY(gRandom->Gaus(position.Y(), 0.2));
  position.SetZ(gRandom->Gaus(position.Z(), 2.));
}

Double_t GetMinimumDistance(TVector3 parent_pos, TVector3 daughter_pos, TVector3 parent_dir, TVector3 daughter_dir,
			    std::array<Double_t,2> z_range, std::array<Double_t,2> &extrapolate_z,
			    std::vector<Double_t> &recon_vertex) {
  if ( recon_vertex.size() != 3 )
    throw std::invalid_argument("size of recon_vertex should be three");

  std::array<Double_t,2> extrapolate_distance;
  TVector3 position_difference = daughter_pos - parent_pos;
  // Almost parallel
  if ( TMath::ACos((parent_dir * daughter_dir) / (parent_dir.Mag() * daughter_dir.Mag())) < 1.e-4) {
    extrapolate_distance.at(0) = (parent_pos.Z() + daughter_pos.Z()) / 2. - parent_pos.Z();
    extrapolate_distance.at(1) = (parent_pos.Z() + daughter_pos.Z()) / 2. - daughter_pos.Z();
  }
  else {
    Double_t delta = parent_dir.Mag2() * daughter_dir.Mag2() - (parent_dir * daughter_dir) * (parent_dir * daughter_dir);
    extrapolate_distance.at(0) = ( 1 * (position_difference * parent_dir) * daughter_dir.Mag2()
				   - (parent_dir * daughter_dir) * (position_difference * daughter_dir) ) / delta;
    extrapolate_distance.at(1) = (-1 * (position_difference * daughter_dir) * parent_dir.Mag2()
				   + (parent_dir * daughter_dir) * (position_difference * parent_dir) ) / delta;
  }
  // z_range.at(0) : small, z_range.at(1) : large
  if ( z_range.at(0) > z_range.at(1) ) {
    std::swap(z_range.at(0), z_range.at(1));
  }

  if ( parent_pos.Z() + extrapolate_distance.at(0) < z_range.at(0) ||
       daughter_pos.Z() + extrapolate_distance.at(1) * daughter_dir.Z() < z_range.at(0)) {
    extrapolate_distance.at(0) = z_range.at(0) - parent_pos.Z();
    extrapolate_distance.at(1) = z_range.at(0) - daughter_pos.Z();
  }
  else if ( parent_pos.Z() + extrapolate_distance.at(0) > z_range.at(1) ||
	    daughter_pos.Z() + extrapolate_distance.at(1) * daughter_dir.Z() > z_range.at(1)) {
    extrapolate_distance.at(0) = z_range.at(1) - parent_pos.Z();
    extrapolate_distance.at(1) = z_range.at(1) - daughter_pos.Z();
  }

  extrapolate_z.at(0) = extrapolate_distance.at(0);
  extrapolate_z.at(1) = extrapolate_distance.at(1);

  TVector3 calculate_parent_position = parent_pos + extrapolate_distance.at(0) * parent_dir;
  TVector3 calculate_daughter_position = daughter_pos + extrapolate_distance.at(1) * daughter_dir;

  TVector3 distance_vec = calculate_parent_position - calculate_daughter_position;
  recon_vertex.at(0) = (calculate_daughter_position + distance_vec).X();
  recon_vertex.at(1) = (calculate_daughter_position + distance_vec).Y();
  recon_vertex.at(2) = (calculate_daughter_position + distance_vec).Z();
  return distance_vec.Mag();

}
