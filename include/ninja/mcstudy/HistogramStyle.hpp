#ifndef NINJA_MC_HISTO_STYLE_HPP
#define NINJA_MC_HISTO_STYLE_HPP

#include <TString.h>

const Int_t num_ninja_mode = 6;
const TString mode_name[num_ninja_mode] = {"CCQE", "2p2h", "CC 1#pi", "CC Multi#pi", "CC Other", "NC"};
const Int_t mode_color[num_ninja_mode] = {424, 624, 394, 395, 401, 408};
const Int_t mode_style[num_ninja_mode] = {1001, 1001, 3006, 3005, 1001, 1001};
const Int_t mode_stack_order[num_ninja_mode] = {5, 4, 3, 2, 0, 1};

Int_t GetNinjaModeId(Int_t mode);

#endif
