#ifndef SROTHMN_EEC_PROJ_AXES_H
#define SROTHMN_EEC_PROJ_AXES_H

#include "usings.h"
#include <vector>

#ifdef CMSSW_GIT_HASH
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#endif

namespace EEC {
    class ProjAxes {
    public:
        ProjAxes(const std::vector<double>& R) noexcept;
        
#ifdef CMSSW_GIT_HASH
        ProjAxes(const edm::ParameterSet& iConfig) noexcept;
        static void fillPSetDescription(edm::ParameterSetDescription& iConfig);
#endif

        const axis R;
    };
};

#endif
