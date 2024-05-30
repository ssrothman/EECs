#ifndef EECS_FAST_H

#include <boost/histogram.hpp>
#include <boost/multi_array.hpp>

#include <cassert>

#include "SRothman/SimonTools/src/jets.h"

#include "fastStructs.h"
#include "faststart.h"

namespace fastEEC{
    template <typename T, bool doPU, bool doTransfer, bool doRes3, bool doRes4, bool doRes4Fixed>
    void fastEEC(result<T>& ans,

                const jet& J, const axisptr& ax, 
                const int order, const normType nt,

                axisptr coarseRLax = nullptr,
                axisptr xiax = nullptr,
                axisptr phiax = nullptr,

                axisptr rax_dipole = nullptr,
                axisptr ctax_dipole = nullptr,

                axisptr rax_tee = nullptr,
                axisptr ctax_tee = nullptr,

                axisptr rax_triangle = nullptr,
                axisptr ctax_triangle = nullptr,
                T shapetol = 0,

                const std::vector<bool>* const PU = nullptr,
                const jet * const J_Reco = nullptr,
                const arma::mat* ptrans = nullptr){

        assert(order >= 2 && order <= 6);

        unsigned NDR = histogram::axis::traits::extent(*ax);

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
        }
        struct res4fixedAxes_t res4fixedAxes;
        if constexpr (doRes4Fixed){
            res4fixedAxes.RL = coarseRLax;
            res4fixedAxes.shapetol = shapetol;
        }

        struct transferInputs<T> tin;
        if constexpr (doTransfer){
            struct jetDetails_t<T> recoJetDetails(*J_Reco, ax, nt);
            tin.recoJet = recoJetDetails;

            tin.adj = adjacency(J, *J_Reco);
            tin.ptrans = *ptrans;
        }

        if (order == 2){
            start<T, doPU, doTransfer, 2, doRes3, doRes4, doRes4Fixed>(
                    ans, jetDetails, NDR, res3axes, res4shapesAxes, res4fixedAxes, tin, PU
            );
        } else if(order == 3){
            start<T, doPU, doTransfer, 3, doRes3, doRes4, doRes4Fixed>(
                    ans, jetDetails, NDR, res3axes, res4shapesAxes, res4fixedAxes, tin, PU
            );
        } else if(order == 4){
            start<T, doPU, doTransfer, 4, doRes3, doRes4, doRes4Fixed>(
                    ans, jetDetails, NDR, res3axes, res4shapesAxes, res4fixedAxes, tin, PU
            );
        } else if(order == 5){
            start<T, doPU, doTransfer, 5, doRes3, doRes4, doRes4Fixed>(
                    ans, jetDetails, NDR, res3axes, res4shapesAxes, res4fixedAxes, tin, PU
            );
        } else if(order == 6){
            start<T, doPU, doTransfer, 6, doRes3, doRes4, doRes4Fixed>(
                    ans, jetDetails, NDR, res3axes, res4shapesAxes, res4fixedAxes, tin, PU
            );
        }
    }
};

#endif
