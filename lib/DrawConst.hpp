#ifndef NINJA_MC_DRAW_CONST_HPP
#define NINJA_MC_DRAW_CONST_HPP

#include <cmath>

const double neutrino_beam_ax = -2.325e-2;
const double neutrino_beam_ay = -8.075e-2;
const double neutrino_beam_thetax = std::atan(neutrino_beam_ax);
const double neutrino_beam_thetay = std::atan(neutrino_beam_ay); // rad

const double muon_mass = 105.658; // MeV/c2
const double pion_mass = 139.571;
const double proton_mass = 938.272;

#endif
