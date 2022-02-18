#include <B2Enum.hh>

#include "HistogramStyle.hpp"

Int_t GetNinjaModeId(Int_t mode) {
  switch (mode) {
    case B2InteractionMode::MODE_CCQE : 
      return 0;
  case B2InteractionMode::MODE_2P2H : 
    return 1;
  case B2InteractionMode::MODE_CC_1PROTON_1PI_PLUS :
  case B2InteractionMode::MODE_CC_1PROTON_1PI_ZERO :
  case B2InteractionMode::MODE_CC_1NEUTRON_1PI_PLUS :
  case B2InteractionMode::MODE_CC_COHERENT_PI_PLUS :
    return 2;
  case B2InteractionMode::MODE_CC_MULTI_PI : 
    return 3;
  case B2InteractionMode::MODE_CC_1PROTON_1GAMMA : 
  case B2InteractionMode::MODE_CC_ETA :
  case B2InteractionMode::MODE_CC_LAMBDA :
  case B2InteractionMode::MODE_CC_DIS :
  case B2InteractionMode::MODE_CC_DIFFRACTIVE_PI :
    return 4;
  case B2InteractionMode::MODE_NC_1NEUTRON_1PI_ZERO :
  case B2InteractionMode::MODE_NC_1PROTON_1PI_ZERO :
  case B2InteractionMode::MODE_NC_1PROTON_1PI_MINUS :
  case B2InteractionMode::MODE_NC_1NEUTRON_1PI_PLUS :
  case B2InteractionMode::MODE_NC_COHERENT_PI_ZERO :
  case B2InteractionMode::MODE_NC_1NEUTRON_1GAMMA :
  case B2InteractionMode::MODE_NC_1PROTON_1GAMMA :
  case B2InteractionMode::MODE_NC_MULTI_PI :
  case B2InteractionMode::MODE_NC_1NEUTRON_ETA :
  case B2InteractionMode::MODE_NC_1PROTON_ETA :
  case B2InteractionMode::MODE_NC_LAMBDA_KAON_ZERO :
  case B2InteractionMode::MODE_NC_LAMBDA_KAON_PLUS :
  case B2InteractionMode::MODE_NC_MESONS :
  case B2InteractionMode::MODE_NC_DIFFRACTIVE_PI :
  case B2InteractionMode::MODE_NC_ELASTIC_1PROTON :
  case B2InteractionMode::MODE_NC_ELASTIC_1NEUTRON :
    return 5;
  default :
    return -1;
  }
}
