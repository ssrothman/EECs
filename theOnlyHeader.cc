#include "run.h"

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

        const fastEEC::axisptr rax_minR = nullptr,
        const fastEEC::axisptr phiax_minR = nullptr,

        const double shapetol = 0,

[[maybe_unused]] const std::vector<bool>* const PU = nullptr,
        const jet* const J_Reco = nullptr,
        const Eigen::MatrixXd* ptrans = nullptr){
            
            fastEEC::runSuperSpecific<double>(
                    ans,
                    J, ax, nt, 
                    maxOrder, flags,
                    coarseRLax, xiax, phiax,
                    rax_dipole, ctax_dipole,
                    rax_tee, ctax_tee,
                    rax_triangle, ctax_triangle,
                    rax_minR, phiax_minR,
                    shapetol,
                    PU, J_Reco, ptrans
            );
}
