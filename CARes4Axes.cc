#include "CARes4Axes.h"

EEC::CARes4Axes::CARes4Axes(
        const std::vector<double>& R,
        const std::vector<double>& r_chain,
        const std::vector<double>& c_chain,
        const std::vector<double>& r_symmetric,
        const std::vector<double>& c_symmetric) noexcept :
    R(R),
    r_chain(r_chain), c_chain(c_chain),
    r_symmetric(r_symmetric), c_symmetric(c_symmetric) {}

#ifdef CMSSW_GIT_HASH

EEC::CARes4Axes::CARes4Axes(const edm::ParameterSet& iConfig) noexcept :
    R(iConfig.getParameter<std::vector<double>>("R")),
    r_chain(iConfig.getParameter<std::vector<double>>("r_chain")),
    c_chain(iConfig.getParameter<std::vector<double>>("c_chain")),
    r_symmetric(iConfig.getParameter<std::vector<double>>("r_symmetric")),
    c_symmetric(iConfig.getParameter<std::vector<double>>("c_symmetric")) {}

void EEC::CARes4Axes::fillPSetDescription(edm::ParameterSetDescription& desc){
    desc.add<std::vector<double>>("R");
    desc.add<std::vector<double>>("r_chain");
    desc.add<std::vector<double>>("c_chain");
    desc.add<std::vector<double>>("r_symmetric");
    desc.add<std::vector<double>>("c_symmetric");
}
    
#endif
