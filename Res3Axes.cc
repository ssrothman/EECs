#include "Res3Axes.h"

EEC::Res3Axes::Res3Axes(
        const std::vector<double>& R,
        const std::vector<double>& r,
        const std::vector<double>& c) noexcept :
    R(R), r(r), c(c) {}

#ifdef CMSSW_GIT_HASH

EEC::Res3Axes::Res3Axes(const edm::ParameterSet& iConfig) noexcept :
    R(iConfig.getParameter<std::vector<double>>("R")),
    r(iConfig.getParameter<std::vector<double>>("r")),
    c(iConfig.getParameter<std::vector<double>>("c")) {}

void EEC::Res3Axes::fillPSetDescription(edm::ParameterSetDescription& desc){
    desc.add<std::vector<double>>("R");
    desc.add<std::vector<double>>("r");
    desc.add<std::vector<double>>("c");
}
    
#endif

