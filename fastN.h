#ifndef EECS_FASTN_H
#define EECS_FASTN_H

#include "util.h"
#include "fastStructs.h"
#include "symFact.h"
#include "prev.h"
#include "res3.h"
#include "res4.h"

namespace fastEEC{
           template <typename T,                               //result type
             bool doPU, bool doTransfer,                   //flags
             bool doRes3, bool doRes4, bool doRes4Fixed,  //more flags
             unsigned maxOrder,                           //max EEC order
             unsigned order                               //current EEC order
    >
    void doN(result_t<T> ans,
             const jetDetails_t<T>& jetDetails,
             unsigned nPart,

             const res3axes_t& res3ax,
             const res4shapesAxes_t& res4ax,
             const res4fixedAxes_t& res4fixedax,

             const transferInputs<T>& tin,
             const vector<bool>* const PU,

             const prev_t<T, order>& prev){

        //printf("(ORDER = %u)\n", order);
        //printf("(prev.is = ");
        //printVec(prev.is);
        //printf(")\n");

        unsigned istart = 0;
        if constexpr (order > 1){
            istart = prev.is[order-2];
        }
        //printf("got istart = %d\n", istart);
        //fflush(stdout);

        for(unsigned inew = istart; inew < nPart; ++inew){
            for(unsigned q=0; q<order-1; ++q){
                //printf("\t");
                //fflush(stdout);
            }
            //printf("inew = %d\n", inew);
            //printf("\n");
            //fflush(stdout);
            prev_t<T, order+1> next;

            if constexpr(order == 1){
                next.partial = jetDetails.Es[inew];
            } else {
                next.partial = prev.partial * jetDetails.Es[inew];
            }
            //printf("got partial\n");
            //fflush(stdout);

            if constexpr (order == 1){
                next.maxDR = 0;
                next.maxDRbin = 0;
            } else {
                std::array<T, order> dRlist;
                std::array<unsigned, order> dRbin_list;
                for(unsigned iold=0; iold<order-1; ++iold){
                    dRlist[iold] = jetDetails.floatDRs[prev.is[iold]][inew];
                    dRbin_list[iold] = jetDetails.dRbins[prev.is[iold]][inew];
                }
                dRlist[order-1] = prev.maxDR;
                dRbin_list[order-1] = prev.maxDRbin;
                //printf("got drlist\n");
                //fflush(stdout);

                auto maxel = max_element(dRlist.begin(), dRlist.end());
                next.maxDR = *maxel;
                next.maxDRbin = dRbin_list[std::distance(dRlist.begin(), maxel)];
            }
            //printf("got maxel\n");
            //fflush(stdout);

            for (unsigned i=0; i<order-1; ++i){
                next.is[i] = prev.is[i];
            }
            next.is[order-1] = inew;
            //printf("got next.is\n");
            //fflush(stdout);
            
            T symfac = getSymfac(next);
            T weight = symfac * next.partial;
            //printf("got weight\n");
            //fflush(stdout);

            if constexpr (doPU){
                next.isPU = prev.isPU || PU->at(inew);
            }

            if constexpr(order > 1){
                //accumulate
                (*ans.wts[order-2])[next.maxDRbin] += weight;

                if constexpr (doPU){
                    if(next.isPU){
                        (*ans.wts_PU[order-2])[next.maxDRbin] += weight;
                    }//end if isPU
                }//end if doPU
                //printf("filled\n");
                //fflush(stdout);
                printf("(");
                printVec(next.is);
                printf(")\n");
                printf("\tsymfac = %g\n", symfac);
                printf("\tweight = %g\n", weight);
            }
            
            //TODO: resolved
            if constexpr (doRes3 && order == 3){
                runRes3<T>(jetDetails, res3ax, next);
                printf("\tres3: %u %u %u\n", next.RL_res3_idx, next.xi_res3_idx, next.phi_res3_idx);
                fflush(stdout);
                (*ans.resolved3)[next.RL_res3_idx][next.xi_res3_idx][next.phi_res3_idx] += weight;
                if constexpr (doPU){
                    if(next.isPU){
                        (*ans.resolved3_PU)[next.RL_res3_idx][next.xi_res3_idx][next.phi_res3_idx] += weight;
                    }
                }
            }

            if constexpr (doRes4 && order == 4){
                runRes4<T>(jetDetails, res4ax, next);
                printf("\tres4: %u %u %u %u\n", next.shape_res4_idx, next.RL_res4_idx, next.r_res4_idx, next.ct_res4_idx);
                fflush(stdout);
                (*ans.resolved4_shapes)[next.shape_res4_idx][next.RL_res4_idx][next.r_res4_idx][next.ct_res4_idx] += weight;
                if constexpr(doPU){
                    if(next.isPU){
                        (*ans.resolved4_shapes_PU)[next.shape_res4_idx][next.RL_res4_idx][next.r_res4_idx][next.ct_res4_idx] += weight;
                    }
                }
            }
            
            if constexpr (doRes4Fixed && order == 4){

            }

            //TODO: transfer
            
            if constexpr (order < maxOrder){
                //printf("punting...\n");
                doN<T, 
                    doPU, doTransfer, 
                    doRes3, doRes4, doRes4Fixed, 
                    maxOrder, order+1
                >(
                        ans, 
                        jetDetails, nPart, 
                        res3ax, res4ax, res4fixedax, 
                        tin, PU, 
                        next
                );
            }//end if order < maxOrder
        }//end loop
    }//end function
}

#endif
