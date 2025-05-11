#include "ProjAxes.h"

EEC::ProjAxes::ProjAxes(const std::vector<double>& R) noexcept :
    R(R) {}

#ifdef CMSSW_GIT_HASH

EEC::ProjAxes::ProjAxes(const edm::ParameterSet& iConfig) noexcept :
    R(iConfig.getParameter<std::vector<double>>("R")) {}

void EEC::ProjAxes::fillPSetDescription(edm::ParameterSetDescription& desc){
    desc.add<std::vector<double>>("R");
}

#endif
