#include "Res4Axes.h"

EEC::Res4Axes::Res4Axes(
        const std::vector<double>& R,
        const std::vector<double>& r_dipole,
        const std::vector<double>& c_dipole,
        const std::vector<double>& r_tee,
        const std::vector<double>& c_tee,
        const std::vector<double>& r_triangle,
        const std::vector<double>& c_triangle) noexcept :
    R(R),
    r_dipole(r_dipole), c_dipole(c_dipole),
    r_tee(r_tee), c_tee(c_tee),
    r_triangle(r_triangle), c_triangle(c_triangle) {}

#ifdef CMSSW_GIT_HASH

EEC::Res4Axes::Res4Axes(const edm::ParameterSet& iConfig) noexcept :
    R(iConfig.getParameter<std::vector<double>>("R")),
    r_dipole(iConfig.getParameter<std::vector<double>>("r_dipole")),
    c_dipole(iConfig.getParameter<std::vector<double>>("c_dipole")),
    r_tee(iConfig.getParameter<std::vector<double>>("r_tee")),
    c_tee(iConfig.getParameter<std::vector<double>>("c_tee")),
    r_triangle(iConfig.getParameter<std::vector<double>>("r_triangle")),
    c_triangle(iConfig.getParameter<std::vector<double>>("c_triangle")) {}

void EEC::Res4Axes::fillPSetDescription(edm::ParameterSetDescription& desc){
    desc.add<std::vector<double>>("R");
    desc.add<std::vector<double>>("r_dipole");
    desc.add<std::vector<double>>("c_dipole");
    desc.add<std::vector<double>>("r_tee");
    desc.add<std::vector<double>>("c_tee");
    desc.add<std::vector<double>>("r_triangle");
    desc.add<std::vector<double>>("c_triangle");
}
    
#endif
