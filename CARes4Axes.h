#ifndef SROTHMAN_EEC_V2_CA_RES4_AXES_H
#define SROTHMAN_EEC_V2_CA_RES4_AXES_H

#include "usings.h"
#include <vector>

#ifdef CMSSW_GIT_HASH
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#endif

namespace EEC {
    class CARes4Axes {
    public:
        CARes4Axes(const std::vector<double>& R,
                 const std::vector<double>& r_chain,
                 const std::vector<double>& c_chain,
                 const std::vector<double>& r_symmetric,
                 const std::vector<double>& c_symmetric) noexcept;

#ifdef CMSSW_GIT_HASH
        CARes4Axes(const edm::ParameterSet& iConfig) noexcept;
        static void fillPSetDescription(edm::ParameterSetDescription& iConfig);
#endif
        
        const axis R;
        const axis r_chain, c_chain;
        const axis r_symmetric, c_symmetric;
    };
};

#endif
