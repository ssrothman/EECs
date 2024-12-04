#ifndef EECS_FAST_RUN_H
#define EECS_FAST_RUN_H

#include <vector>
#include "fastStructs.h"
#include "fastN.h"
#include "clear2.h"
#include "util.h"
#include "SRothman/SimonTools/src/recursive_reduce.h"
#include "flags.h"

namespace fastEEC{

    template <typename T,
             bool doPU, bool doTransfer,
             bool doRes3, bool doRes4,
             unsigned maxOrder>
    void run(result_t<T>& ans,

             const jet& J,
             const axisptr& ax,
             const normType nt,

             const axisptr coarseRLax = nullptr,
             const axisptr xiax = nullptr,
             const axisptr phiax = nullptr,

             const axisptr rax_dipole = nullptr,
             const axisptr ctax_dipole = nullptr,
             
             const axisptr rax_tee = nullptr,
             const axisptr ctax_tee = nullptr,

             const axisptr rax_triangle = nullptr,
             const axisptr ctax_triangle = nullptr,

             const axisptr rax_minR = nullptr,
             const axisptr phiax_minR = nullptr,
             const axisptr phidiffax_minR = nullptr,
             const axisptr thetaax_minR = nullptr,

             const T shapetol = 0,

             [[maybe_unused]] const std::vector<bool>* const PU = nullptr,
             const jet* const J_Reco = nullptr,
             const Eigen::MatrixXd* ptrans = nullptr){

        struct jetDetails_t<T> jetDetails(J, ax, nt);

        struct res3axes_t res3axes;
        if constexpr (doRes3){
            res3axes.RL = coarseRLax;
            res3axes.xi = xiax;
            res3axes.phi = phiax;
        }
        struct res4shapesAxes_t res4shapesAxes;
        if constexpr (doRes4){
            res4shapesAxes.RL = coarseRLax;

            res4shapesAxes.r_dipole = rax_dipole;
            res4shapesAxes.ct_dipole = ctax_dipole;

            res4shapesAxes.r_tee = rax_tee;
            res4shapesAxes.ct_tee = ctax_tee;

            res4shapesAxes.r_triangle = rax_triangle;
            res4shapesAxes.ct_triangle = ctax_triangle;

            res4shapesAxes.shapetol = shapetol;

            res4shapesAxes.r_minR = rax_minR;
            res4shapesAxes.phi_minR = phiax_minR;
            res4shapesAxes.phidiff_minR = phidiffax_minR;
            res4shapesAxes.theta_minR = thetaax_minR;
        }

        struct transferInputs<T> tin;
        if constexpr (doTransfer){
            tin.setup(J_Reco, ptrans, ax, nt);

            Eigen::VectorXd ptrec = J_Reco->ptvec();
            Eigen::VectorXd ptgen = J.ptvec();
            Eigen::VectorXd ptfwd = (*tin.ptrans).transpose() * ptgen;

            //std::cout << *tin.ptrans << std::endl;
            for (unsigned iReco=0; iReco<J_Reco->particles.size(); ++iReco){
                for(unsigned iGen=0; iGen<J.particles.size(); ++iGen){
                    if (ptfwd[iReco] > 0){
                        (*tin.ptrans)(iGen,iReco) = (*tin.ptrans)(iGen,iReco) * ptrec(iReco) / ptfwd[iReco];
                    }
                }
            }

            switch(nt){
                case normType::RAWPT:
                    *tin.ptrans *= J.rawpt / J_Reco->rawpt;
                    break;
                case normType::SUMPT:
                    *tin.ptrans *= J.sumpt / J_Reco->sumpt;
                    break;
                case normType::CORRPT:
                    *tin.ptrans *= J.pt / J_Reco->pt;
                    break;
                default:
                    assert(false);
            }
        }

        clear<T, doPU, doTransfer>(ans, ax, res3axes, res4shapesAxes);

        unsigned nPart = jetDetails.Es.size();

        prev_t<T, 1> prev;

        doN<T, doPU, doTransfer,
            doRes3, doRes4, 
            maxOrder, 1, 0>(
                ans,
                jetDetails,
                nPart,

                res3axes,
                res4shapesAxes,

                tin,
                PU,
                
                prev
        );
    }

    template <typename T, unsigned flags>
    void runSpecific(result_t<T>& ans,
            
            const jet& J,
            const axisptr& ax,
            const normType nt,
            const unsigned maxOrder,

            const axisptr coarseRLax = nullptr,
            const axisptr xiax = nullptr,
            const axisptr phiax = nullptr,

            const axisptr rax_dipole = nullptr,
            const axisptr ctax_dipole = nullptr,
             
            const axisptr rax_tee = nullptr,
            const axisptr ctax_tee = nullptr,

            const axisptr rax_triangle = nullptr,
            const axisptr ctax_triangle = nullptr,

             const axisptr rax_minR = nullptr,
             const axisptr phiax_minR = nullptr,
             const axisptr phidiffax_minR = nullptr,
             const axisptr thetaax_minR = nullptr,

            const T shapetol = 0,

            [[maybe_unused]] const std::vector<bool>* const PU = nullptr,
            const jet* const J_Reco = nullptr,
            const Eigen::MatrixXd* ptrans = nullptr){
        switch(maxOrder){
            case 2:
                run<T, bool(flags & DOPU), 
                       bool(flags & DOTRANSFER), 
                       bool(flags & DORES3),
                       bool(flags & DORES4), 2>(
                    ans,
                    J, ax, nt,
                    coarseRLax, xiax, phiax,
                    rax_dipole, ctax_dipole,
                    rax_tee, ctax_tee,
                    rax_triangle, ctax_triangle,
                    rax_minR, phiax_minR, phidiffax_minR, thetaax_minR,
                    shapetol,
                    PU, J_Reco, ptrans
                );
                break;
            case 3:
                run<T, bool(flags & DOPU), 
                       bool(flags & DOTRANSFER), 
                       bool(flags & DORES3),
                       bool(flags & DORES4), 3>(
                    ans,
                    J, ax, nt,
                    coarseRLax, xiax, phiax,
                    rax_dipole, ctax_dipole,
                    rax_tee, ctax_tee,
                    rax_triangle, ctax_triangle,
                    rax_minR, phiax_minR, phidiffax_minR, thetaax_minR,
                    shapetol,
                    PU, J_Reco, ptrans
                );
                break;
            case 4:
                run<T, bool(flags & DOPU), 
                       bool(flags & DOTRANSFER), 
                       bool(flags & DORES3),
                       bool(flags & DORES4), 4>(
                    ans,
                    J, ax, nt,
                    coarseRLax, xiax, phiax,
                    rax_dipole, ctax_dipole,
                    rax_tee, ctax_tee,
                    rax_triangle, ctax_triangle,
                    rax_minR, phiax_minR, phidiffax_minR, thetaax_minR,
                    shapetol,
                    PU, J_Reco, ptrans
                );
                break;
            case 5:
                run<T, bool(flags & DOPU), 
                       bool(flags & DOTRANSFER), 
                       bool(flags & DORES3),
                       bool(flags & DORES4), 5>(
                    ans,
                    J, ax, nt,
                    coarseRLax, xiax, phiax,
                    rax_dipole, ctax_dipole,
                    rax_tee, ctax_tee,
                    rax_triangle, ctax_triangle,
                    rax_minR, phiax_minR, phidiffax_minR, thetaax_minR,
                    shapetol,
                    PU, J_Reco, ptrans
                );
                break;
            case 6:
                run<T, bool(flags & DOPU), 
                       bool(flags & DOTRANSFER), 
                       bool(flags & DORES3),
                       bool(flags & DORES4), 6>(
                    ans,
                    J, ax, nt,
                    coarseRLax, xiax, phiax,
                    rax_dipole, ctax_dipole,
                    rax_tee, ctax_tee,
                    rax_triangle, ctax_triangle,
                    rax_minR, phiax_minR, phidiffax_minR, thetaax_minR,
                    shapetol,
                    PU, J_Reco, ptrans
                );
                break;
            default:
                assert(false);
        }
    }

    template <typename T>
    void runSuperSpecific(result_t<T>& ans,

            const jet& J,
            const axisptr& ax,
            const normType nt,
            const unsigned maxOrder,
            const unsigned flags,

            const axisptr coarseRLax = nullptr,
            const axisptr xiax = nullptr,
            const axisptr phiax = nullptr,

            const axisptr rax_dipole = nullptr,
            const axisptr ctax_dipole = nullptr,
             
            const axisptr rax_tee = nullptr,
            const axisptr ctax_tee = nullptr,

            const axisptr rax_triangle = nullptr,
            const axisptr ctax_triangle = nullptr,

             const axisptr rax_minR = nullptr,
             const axisptr phiax_minR = nullptr,
             const axisptr phidiffax_minR = nullptr,
             const axisptr thetaax_minR = nullptr,

            const T shapetol = 0,

             [[maybe_unused]] const std::vector<bool>* const PU = nullptr,
             const jet* const J_Reco = nullptr,
             const Eigen::MatrixXd* ptrans = nullptr){


        switch(flags){
            case 0:
                runSpecific<T, 0>(
                    ans,
                    J, ax, nt, maxOrder,
                    coarseRLax, xiax, phiax,
                    rax_dipole, ctax_dipole,
                    rax_tee, ctax_tee,
                    rax_triangle, ctax_triangle,
                    rax_minR, phiax_minR, phidiffax_minR, thetaax_minR,
                    shapetol,
                    PU, J_Reco, ptrans
                );
                break;
            case 1:
                runSpecific<T, 1>(
                    ans,
                    J, ax, nt, maxOrder,
                    coarseRLax, xiax, phiax,
                    rax_dipole, ctax_dipole,
                    rax_tee, ctax_tee,
                    rax_triangle, ctax_triangle,
                    rax_minR, phiax_minR, phidiffax_minR, thetaax_minR,
                    shapetol,
                    PU, J_Reco, ptrans
                );
                break;
            case 2:
                runSpecific<T, 2>(
                    ans,
                    J, ax, nt, maxOrder,
                    coarseRLax, xiax, phiax,
                    rax_dipole, ctax_dipole,
                    rax_tee, ctax_tee,
                    rax_triangle, ctax_triangle,
                    rax_minR, phiax_minR, phidiffax_minR, thetaax_minR,
                    shapetol,
                    PU, J_Reco, ptrans
                );
                break;
            case 3:
                runSpecific<T, 3>(
                    ans,
                    J, ax, nt, maxOrder,
                    coarseRLax, xiax, phiax,
                    rax_dipole, ctax_dipole,
                    rax_tee, ctax_tee,
                    rax_triangle, ctax_triangle,
                    rax_minR, phiax_minR, phidiffax_minR, thetaax_minR,
                    shapetol,
                    PU, J_Reco, ptrans
                );
                break;
            case 4:
                runSpecific<T, 4>(
                    ans,
                    J, ax, nt, maxOrder,
                    coarseRLax, xiax, phiax,
                    rax_dipole, ctax_dipole,
                    rax_tee, ctax_tee,
                    rax_triangle, ctax_triangle,
                    rax_minR, phiax_minR, phidiffax_minR, thetaax_minR,
                    shapetol,
                    PU, J_Reco, ptrans
                );
                break;
            case 5:
                runSpecific<T, 5>(
                    ans,
                    J, ax, nt, maxOrder,
                    coarseRLax, xiax, phiax,
                    rax_dipole, ctax_dipole,
                    rax_tee, ctax_tee,
                    rax_triangle, ctax_triangle,
                    rax_minR, phiax_minR, phidiffax_minR, thetaax_minR,
                    shapetol,
                    PU, J_Reco, ptrans
                );
                break;
            case 6:
                runSpecific<T, 6>(
                    ans,
                    J, ax, nt, maxOrder,
                    coarseRLax, xiax, phiax,
                    rax_dipole, ctax_dipole,
                    rax_tee, ctax_tee,
                    rax_triangle, ctax_triangle,
                    rax_minR, phiax_minR, phidiffax_minR, thetaax_minR,
                    shapetol,
                    PU, J_Reco, ptrans
                );
                break;
            case 7:
                runSpecific<T, 7>(
                    ans,
                    J, ax, nt, maxOrder,
                    coarseRLax, xiax, phiax,
                    rax_dipole, ctax_dipole,
                    rax_tee, ctax_tee,
                    rax_triangle, ctax_triangle,
                    rax_minR, phiax_minR, phidiffax_minR, thetaax_minR,
                    shapetol,
                    PU, J_Reco, ptrans
                );
                break;
            case 8:
                runSpecific<T, 8>(
                    ans,
                    J, ax, nt, maxOrder,
                    coarseRLax, xiax, phiax,
                    rax_dipole, ctax_dipole,
                    rax_tee, ctax_tee,
                    rax_triangle, ctax_triangle,
                    rax_minR, phiax_minR, phidiffax_minR, thetaax_minR,
                    shapetol,
                    PU, J_Reco, ptrans
                );
                break;
            case 9:
                runSpecific<T, 9>(
                    ans,
                    J, ax, nt, maxOrder,
                    coarseRLax, xiax, phiax,
                    rax_dipole, ctax_dipole,
                    rax_tee, ctax_tee,
                    rax_triangle, ctax_triangle,
                    rax_minR, phiax_minR, phidiffax_minR, thetaax_minR,
                    shapetol,
                    PU, J_Reco, ptrans
                );
                break;
            case 10:
                runSpecific<T, 10>(
                    ans,
                    J, ax, nt, maxOrder,
                    coarseRLax, xiax, phiax,
                    rax_dipole, ctax_dipole,
                    rax_tee, ctax_tee,
                    rax_triangle, ctax_triangle,
                    rax_minR, phiax_minR, phidiffax_minR, thetaax_minR,
                    shapetol,
                    PU, J_Reco, ptrans
                );
                break;
            case 11:
                runSpecific<T, 11>(
                    ans,
                    J, ax, nt, maxOrder,
                    coarseRLax, xiax, phiax,
                    rax_dipole, ctax_dipole,
                    rax_tee, ctax_tee,
                    rax_triangle, ctax_triangle,
                    rax_minR, phiax_minR, phidiffax_minR, thetaax_minR,
                    shapetol,
                    PU, J_Reco, ptrans
                );
                break;
            case 12:
                runSpecific<T, 12>(
                    ans,
                    J, ax, nt, maxOrder,
                    coarseRLax, xiax, phiax,
                    rax_dipole, ctax_dipole,
                    rax_tee, ctax_tee,
                    rax_triangle, ctax_triangle,
                    rax_minR, phiax_minR, phidiffax_minR, thetaax_minR,
                    shapetol,
                    PU, J_Reco, ptrans
                );
                break;
            case 13:
                runSpecific<T, 13>(
                    ans,
                    J, ax, nt, maxOrder,
                    coarseRLax, xiax, phiax,
                    rax_dipole, ctax_dipole,
                    rax_tee, ctax_tee,
                    rax_triangle, ctax_triangle,
                    rax_minR, phiax_minR, phidiffax_minR, thetaax_minR,
                    shapetol,
                    PU, J_Reco, ptrans
                );
                break;
            case 14:
                runSpecific<T, 14>(
                    ans,
                    J, ax, nt, maxOrder,
                    coarseRLax, xiax, phiax,
                    rax_dipole, ctax_dipole,
                    rax_tee, ctax_tee,
                    rax_triangle, ctax_triangle,
                    rax_minR, phiax_minR, phidiffax_minR, thetaax_minR,
                    shapetol,
                    PU, J_Reco, ptrans
                );
                break;
            case 15:
                runSpecific<T, 15>(
                    ans,
                    J, ax, nt, maxOrder,
                    coarseRLax, xiax, phiax,
                    rax_dipole, ctax_dipole,
                    rax_tee, ctax_tee,
                    rax_triangle, ctax_triangle,
                    rax_minR, phiax_minR, phidiffax_minR, thetaax_minR,
                    shapetol,
                    PU, J_Reco, ptrans
                );
                break;
            case 16:
                runSpecific<T, 16>(
                    ans,
                    J, ax, nt, maxOrder,
                    coarseRLax, xiax, phiax,
                    rax_dipole, ctax_dipole,
                    rax_tee, ctax_tee,
                    rax_triangle, ctax_triangle,
                    rax_minR, phiax_minR, phidiffax_minR, thetaax_minR,
                    shapetol,
                    PU, J_Reco, ptrans
                );
                break;
            case 17:
                runSpecific<T, 17>(
                    ans,
                    J, ax, nt, maxOrder,
                    coarseRLax, xiax, phiax,
                    rax_dipole, ctax_dipole,
                    rax_tee, ctax_tee,
                    rax_triangle, ctax_triangle,
                    rax_minR, phiax_minR, phidiffax_minR, thetaax_minR,
                    shapetol,
                    PU, J_Reco, ptrans
                );
                break;
            case 18:
                runSpecific<T, 18>(
                    ans,
                    J, ax, nt, maxOrder,
                    coarseRLax, xiax, phiax,
                    rax_dipole, ctax_dipole,
                    rax_tee, ctax_tee,
                    rax_triangle, ctax_triangle,
                    rax_minR, phiax_minR, phidiffax_minR, thetaax_minR,
                    shapetol,
                    PU, J_Reco, ptrans
                );
                break;
            case 19:
                runSpecific<T, 19>(
                    ans,
                    J, ax, nt, maxOrder,
                    coarseRLax, xiax, phiax,
                    rax_dipole, ctax_dipole,
                    rax_tee, ctax_tee,
                    rax_triangle, ctax_triangle,
                    rax_minR, phiax_minR, phidiffax_minR, thetaax_minR,
                    shapetol,
                    PU, J_Reco, ptrans
                );
                break;
            case 20:
                runSpecific<T, 20>(
                    ans,
                    J, ax, nt, maxOrder,
                    coarseRLax, xiax, phiax,
                    rax_dipole, ctax_dipole,
                    rax_tee, ctax_tee,
                    rax_triangle, ctax_triangle,
                    rax_minR, phiax_minR, phidiffax_minR, thetaax_minR,
                    shapetol,
                    PU, J_Reco, ptrans
                );
                break;
            case 21:
                runSpecific<T, 21>(
                    ans,
                    J, ax, nt, maxOrder,
                    coarseRLax, xiax, phiax,
                    rax_dipole, ctax_dipole,
                    rax_tee, ctax_tee,
                    rax_triangle, ctax_triangle,
                    rax_minR, phiax_minR, phidiffax_minR, thetaax_minR,
                    shapetol,
                    PU, J_Reco, ptrans
                );
                break;
            case 22:
                runSpecific<T, 22>(
                    ans,
                    J, ax, nt, maxOrder,
                    coarseRLax, xiax, phiax,
                    rax_dipole, ctax_dipole,
                    rax_tee, ctax_tee,
                    rax_triangle, ctax_triangle,
                    rax_minR, phiax_minR, phidiffax_minR, thetaax_minR,
                    shapetol,
                    PU, J_Reco, ptrans
                );
                break;
            case 23:
                runSpecific<T, 23>(
                    ans,
                    J, ax, nt, maxOrder,
                    coarseRLax, xiax, phiax,
                    rax_dipole, ctax_dipole,
                    rax_tee, ctax_tee,
                    rax_triangle, ctax_triangle,
                    rax_minR, phiax_minR, phidiffax_minR, thetaax_minR,
                    shapetol,
                    PU, J_Reco, ptrans
                );
                break;
            case 24:
                runSpecific<T, 24>(
                    ans,
                    J, ax, nt, maxOrder,
                    coarseRLax, xiax, phiax,
                    rax_dipole, ctax_dipole,
                    rax_tee, ctax_tee,
                    rax_triangle, ctax_triangle,
                    rax_minR, phiax_minR, phidiffax_minR, thetaax_minR,
                    shapetol,
                    PU, J_Reco, ptrans
                );
                break;
            case 25:
                runSpecific<T, 25>(
                    ans,
                    J, ax, nt, maxOrder,
                    coarseRLax, xiax, phiax,
                    rax_dipole, ctax_dipole,
                    rax_tee, ctax_tee,
                    rax_triangle, ctax_triangle,
                    rax_minR, phiax_minR, phidiffax_minR, thetaax_minR,
                    shapetol,
                    PU, J_Reco, ptrans
                );
                break;
            case 26:
                runSpecific<T, 26>(
                    ans,
                    J, ax, nt, maxOrder,
                    coarseRLax, xiax, phiax,
                    rax_dipole, ctax_dipole,
                    rax_tee, ctax_tee,
                    rax_triangle, ctax_triangle,
                    rax_minR, phiax_minR, phidiffax_minR, thetaax_minR,
                    shapetol,
                    PU, J_Reco, ptrans
                );
                break;
            case 27:
                runSpecific<T, 27>(
                    ans,
                    J, ax, nt, maxOrder,
                    coarseRLax, xiax, phiax,
                    rax_dipole, ctax_dipole,
                    rax_tee, ctax_tee,
                    rax_triangle, ctax_triangle,
                    rax_minR, phiax_minR, phidiffax_minR, thetaax_minR,
                    shapetol,
                    PU, J_Reco, ptrans
                );
                break;
            case 28:
                runSpecific<T, 28>(
                    ans,
                    J, ax, nt, maxOrder,
                    coarseRLax, xiax, phiax,
                    rax_dipole, ctax_dipole,
                    rax_tee, ctax_tee,
                    rax_triangle, ctax_triangle,
                    rax_minR, phiax_minR, phidiffax_minR, thetaax_minR,
                    shapetol,
                    PU, J_Reco, ptrans
                );
                break;
            case 29:
                runSpecific<T, 29>(
                    ans,
                    J, ax, nt, maxOrder,
                    coarseRLax, xiax, phiax,
                    rax_dipole, ctax_dipole,
                    rax_tee, ctax_tee,
                    rax_triangle, ctax_triangle,
                    rax_minR, phiax_minR, phidiffax_minR, thetaax_minR,
                    shapetol,
                    PU, J_Reco, ptrans
                );
                break;
            case 30:
                runSpecific<T, 30>(
                    ans,
                    J, ax, nt, maxOrder,
                    coarseRLax, xiax, phiax,
                    rax_dipole, ctax_dipole,
                    rax_tee, ctax_tee,
                    rax_triangle, ctax_triangle,
                    rax_minR, phiax_minR, phidiffax_minR, thetaax_minR,
                    shapetol,
                    PU, J_Reco, ptrans
                );
                break;
            case 31:
                runSpecific<T, 31>(
                    ans,
                    J, ax, nt, maxOrder,
                    coarseRLax, xiax, phiax,
                    rax_dipole, ctax_dipole,
                    rax_tee, ctax_tee,
                    rax_triangle, ctax_triangle,
                    rax_minR, phiax_minR, phidiffax_minR, thetaax_minR,
                    shapetol,
                    PU, J_Reco, ptrans
                );
                break;
            default:
                assert(false);
        };
    }
}

#endif
