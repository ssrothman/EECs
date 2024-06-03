#ifndef EECS_FAST_RUN_H
#define EECS_FAST_RUN_H

#include <vector>
#include "fastStructs.h"
#include "fastN.h"
#include "clear2.h"
#include "util.h"
#include "SRothman/SimonTools/src/recursive_reduce.h"

namespace fastEEC{
    template <typename T,
             bool doPU, bool doTransfer,
             bool doRes3, bool doRes4, bool doRes4Fixed,
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

             const T shapetol = 0,

             [[maybe_unused]] const std::vector<bool>* const PU = nullptr,
             const jet* const J_Reco = nullptr,
             const arma::mat* ptrans = nullptr){

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

        clear(ans, ax, res3axes, res4shapesAxes, res4fixedAxes);

        unsigned nPart = jetDetails.Es.size();

        prev_t<T, 1> prev;

        doN<T, doPU, doTransfer,
            doRes3, doRes4, doRes4Fixed,
            maxOrder, 1>(
                ans,
                jetDetails,
                nPart,

                res3axes,
                res4shapesAxes,
                res4fixedAxes,

                tin,
                nullptr,
                
                prev
        );
        printf("RESULTS\n");
        for(unsigned order=0; order < 5; ++order){
            printf("Order %d\n", order);
            printf("\tsumwt = %f\n", recursive_reduce(*ans.wts[order], 0.));
        }
        printf("RES3\n");
        printf("\tsumwt = %f\n", recursive_reduce(*ans.resolved3, 0.));
        fflush(stdout);
        printf("RES4\n");
        printf("\tsumwt = %f\n", recursive_reduce(*ans.resolved4_shapes, 0.));
        printf("\t\tshape0 = %f\n", recursive_reduce((*ans.resolved4_shapes)[0], 0.));
        printf("\t\tshape1 = %f\n", recursive_reduce((*ans.resolved4_shapes)[1], 0.));
        printf("\t\tshape2 = %f\n", recursive_reduce((*ans.resolved4_shapes)[2], 0.));
        printf("\t\tshape3 = %f\n", recursive_reduce((*ans.resolved4_shapes)[3], 0.));
    }
}

#endif
