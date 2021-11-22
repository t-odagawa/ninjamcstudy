#ifndef MINIMUM_DISTANCE_HPP
#define MINIMUM_DISTANCE_HPP

Int_t GetInteractionEcc(const B2VertexSummary &primary_vertex_summary);

Int_t GetVertexPlate(Double_t z_pos);

void CalcPosInEccCoordinate(TVector3 &vertex_position, Int_t ecc_id);

TVector3 SmearPosition(TVector3 true_position);

TVector3 SmearTangent(TVector3 true_tangent);

Double_t RadialSmearFunction(Double_t tangent);

Double_t GetMinimumDistance(TVector3 parent_pos, TVector3 daughter_pos,
			    TVector3 parent_die, TVector3 daughter_dir,
			    std::array<Double_t, 2 > z_range,
			    std::array<Double_t, 2 > &exrapolate_z,
			    std::vector<Double_t> &recon_vertex);

#endif
