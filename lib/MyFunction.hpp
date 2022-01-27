#ifndef NINJA_MC_STUDY_MYFUNC_HPP
#define NINJA_MC_STUDY_MYFUNC_HPP

#include <vector>
#include <array>

#include <B2VertexSummary.hh>

#include <TVector3.h>

Int_t GetInteractionEcc(const B2VertexSummary &primary_vertex_summary);

Int_t GetVertexPlate(Double_t z_pos);

void CalcPosInEccCoordinate(TVector3 &position, Int_t ecc_id, Bool_t absolute_flag);

void SmearPosition(TVector3 &position);

Double_t GetMinimumDistance(TVector3 parent_pos, TVector3 daughter_pos, TVector3 parent_dir, TVector3 daughter_dir,
			    std::array<Double_t ,2> z_range, std::array<Double_t, 2> &extrapolate_z,
			    std::vector<Double_t> &recon_vertex);

#endif
