#ifndef SROTHMAN_EEC_V2_RES3_AXES_H
#define SROTHMAN_EEC_V2_RES3_AXES_H

#include "usings.h"
#include <vector>

#ifdef CMSSW_GIT_HASH
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#endif

namespace EEC {
    class Res3Axes {
    public:
        Res3Axes(const std::vector<double>& R,
                 const std::vector<double>& r,
                 const std::vector<double>& c) noexcept;

#ifdef CMSSW_GIT_HASH
        Res3Axes(const edm::ParameterSet& iConfig) noexcept;
        static void fillPSetDescription(edm::ParameterSetDescription& iConfig);
#endif

        const axis R, r, c;
    };
};

#endif
