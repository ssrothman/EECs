#ifndef SROTHMAN_EECS_ONLYHEADER_H
#define SROTHMAN_EECS_ONLYHEADER_H

#include "fastStructs.h"
#include "flags.h"

void runFastEEC(fastEEC::result_t<double>& ans,

        const jet& J,
        const fastEEC::axisptr& ax,
        const fastEEC::normType nt,
        const unsigned maxOrder,
        const unsigned flags,
        
        const fastEEC::axisptr coarseRLax = nullptr,
        const fastEEC::axisptr xiax = nullptr,
        const fastEEC::axisptr phiax = nullptr,

        const fastEEC::axisptr rax_dipole = nullptr,
        const fastEEC::axisptr ctax_dipole = nullptr,

        const fastEEC::axisptr rax_tee = nullptr,
        const fastEEC::axisptr ctax_tee = nullptr,

        const fastEEC::axisptr rax_triangle = nullptr,
        const fastEEC::axisptr ctax_triangle = nullptr,

        const double shapetol = 0,

        [[maybe_unused]] const std::vector<bool>* const PU = nullptr,
        const jet* const J_Reco = nullptr,
        const arma::mat* ptrans = nullptr);

#endif
