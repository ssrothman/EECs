#ifndef SROTHMAN_EEC_PROJ_CALCULATOR_H
#define SROTHMAN_EEC_PROJ_CALCULATOR_H

#include "ProjAxes.h"
#include "EECjet.h"

#include "ProjResult.h"
#include "ProjTransferResult.h"

#ifdef CMSSW_GIT_HASH
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#endif

namespace EEC{
    class ProjCalculator {
    public:
        ProjCalculator(const std::vector<double>& R,
                       const normType nt):
            axes(R),
            nt(nt) {}

#ifdef CMSSW_GIT_HASH
        ProjCalculator(const edm::ParameterSet& iConfig);
        static void fillPSetDescription(edm::ParameterSetDescription& desc);
#endif

        void compute(const simon::jet& J,
                     ProjResult_Array& result) const noexcept;

        void compute(const simon::jet& J,
                     ProjResult_Vector& result) const noexcept;

        void compute(const simon::jet& J,
                     ProjResult_Unbinned& result) const noexcept;

        void compute_matched(const simon::jet& J,
                             const std::vector<bool>& matched,
                             ProjResult_Array& result,
                             ProjResult_Array& unmatched) const noexcept;

        void compute_matched(const simon::jet& J,
                             const std::vector<bool>& matched,
                             ProjResult_Vector& result,
                             ProjResult_Vector& unmatched) const noexcept;

        void compute_matched(const simon::jet& J,
                             const std::vector<bool>& matched,
                             ProjResult_Unbinned& result,
                             ProjResult_Unbinned& unmatched) const noexcept;

        const ProjAxes& get_axes() const noexcept {
            return axes;
        }
    private:
        ProjAxes axes;
        normType nt;
    };
};

#endif
