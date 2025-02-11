#ifndef SROTHMAN_EEC_V2_RES4_AXES_H
#define SROTHMAN_EEC_V2_RES4_AXES_H

#include "usings.h"
#include <vector>

#ifdef CMSSW_GIT_HASH
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#endif

namespace EEC {
    class Res4Axes {
    public:
        Res4Axes(const std::vector<double>& R,
                 const std::vector<double>& r_dipole,
                 const std::vector<double>& c_dipole,
                 const std::vector<double>& r_tee,
                 const std::vector<double>& c_tee,
                 const std::vector<double>& r_triangle,
                 const std::vector<double>& c_triangle) noexcept;

#ifdef CMSSW_GIT_HASH
        Res4Axes(const edm::ParameterSet& iConfig) noexcept;
        static void fillPSetDescription(edm::ParameterSetDescription& iConfig);
#endif
        
        const axis R;
        const axis r_dipole, c_dipole;
        const axis r_tee, c_tee;
        const axis r_triangle, c_triangle;
    };
};

#endif
