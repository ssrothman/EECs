#ifndef EECS_FASTN_INNER_H
#define EECS_FASTN_INNER_H

#include "util.h"
#include "fastStructs.h"
#include "symFact.h"
#include "prev.h"
#include "res3.h"
#include "res4.h"
#include <assert.h>

namespace fastEEC{
    template <typename T,                                 //result type
             bool doPU, bool doTransfer,                  //flags
             bool doRes3, bool doRes4, bool doRes4Fixed,  //more flags
             unsigned maxOrder,                           //max EEC order
             unsigned order                               //current EEC order
    >
    void doN(result_t<T>& ans,
             const jetDetails_t<T>& jetDetails,
             const unsigned nPart,

             const res3axes_t& res3ax,
             const res4shapesAxes_t& res4ax,
             const res4fixedAxes_t& res4fixedax,

             const transferInputs<T>& tin,
             const vector<bool>* const PU,

             const prev_t<T, order>& prev);
    

    template <typename T,                                 //result type
             bool doPU, bool doTransfer,                  //flags
             bool doRes3, bool doRes4, bool doRes4Fixed,  //more flags
             unsigned maxOrder,                           //max EEC order
             unsigned order                               //current EEC order
    >
    void doNinner(result_t<T>& ans,
                  const unsigned inew,
                  
                  const jetDetails_t<T>& jetDetails,
                  const unsigned nPart,

                  const res3axes_t& res3ax,
                  const res4shapesAxes_t& res4ax,
                  const res4fixedAxes_t& res4fixedax,

                  const transferInputs<T>& tin,
                  const vector<bool>* const PU,

                  const prev_t<T, order>& prev,
                  const T symfac,
                  const unsigned symfacIndex){

        prev_t<T, order+1> next;

        if constexpr(order == 1){
            next.partial = jetDetails.Es[inew];
        } else {
            next.partial = prev.partial * jetDetails.Es[inew];
        }

        next.symfacIndex = symfacIndex;

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

            auto maxel = max_element(dRlist.begin(), dRlist.end());
            next.maxDR = *maxel;
            next.maxDRbin = dRbin_list[std::distance(dRlist.begin(), maxel)];
        }

        for (unsigned i=0; i<order-1; ++i){
            next.is[i] = prev.is[i];
        }
        next.is[order-1] = inew;

        //printVec(next.is);
        //printf("\n");
        //fflush(stdout);

        if constexpr (doPU){
            //printf("about to do PU\n");
            next.isPU = prev.isPU || PU->at(inew);
            //printf("\tset isPU PU\n");
            //fflush(stdout);
        }

        if constexpr(order > 1){
            T weight = symfac * next.partial;

            for (unsigned i=0; i<order-2; ++i){
                //printf("passing forward %u\n", i);
                //printf("\t= %g\n", prev.wts[i]);
                next.wts[i] = prev.wts[i];
                next.dRbins[i] = prev.dRbins[i];
            }
            //printf("setting this weight %u\n", order-2);
            //printf("\t= %g\n", weight);
            next.wts[order-2] = weight;
            next.dRbins[order-2] = next.maxDRbin;
            //printf("\tgot weights\n");
            //fflush(stdout);

            //accumulate
            (*ans.wts[order-2])[next.maxDRbin] += weight;
            //printf("\taccumulated weight\n");
            //fflush(stdout);

            if constexpr (doPU){
                if(next.isPU){
                    (*ans.wts_PU[order-2])[next.maxDRbin] += weight;
                    //printf("\taccumulated PU weight\n");
                    //fflush(stdout);
                }//end if isPU
            }//end if doPU
            /*
            printf("(");
            printVec(next.is);
            printf(")\n");
            printf("\tsymfac = %g\n", symfac);
            printf("\tweight = %g\n", weight);
            printf("\tnextwts: ");
            printVec(next.wts);
            printf("\n");
            printf("\tnextDRs: ");
            printVec(next.dRbins);
            printf("\n");
            */
            if constexpr (doRes3 && order == 3){
                runRes3<T>(jetDetails, res3ax, next);
                //printf("\tres3: %u %u %u\n", next.RL_res3_idx, next.xi_res3_idx, next.phi_res3_idx);
                //fflush(stdout);
                (*ans.resolved3)[next.RL_res3_idx][next.xi_res3_idx][next.phi_res3_idx] += weight;
                if constexpr (doPU){
                    if(next.isPU){
                        (*ans.resolved3_PU)[next.RL_res3_idx][next.xi_res3_idx][next.phi_res3_idx] += weight;
                    }
                }
                //printf("\thandled res3\n");
                //fflush(stdout);
            }

            if constexpr (doRes4 && order == 4){
                runRes4<T>(jetDetails, res4ax, next);
                //printf("\tres4: %u %u %u %u\n", next.shape_res4_idx, next.RL_res4_idx, next.r_res4_idx, next.ct_res4_idx);
                //fflush(stdout);
                for(unsigned q=0; q<3; ++q){
                    ans.resolved4_shapes->fill(weight,
                            next.RL_res4_idx[q],
                            next.shape_res4_idx[q],
                            next.r_res4_idx[q],
                            next.ct_res4_idx[q]);
                }

                if constexpr(doPU){
                    if(next.isPU){
                        for(unsigned q=0; q<1; ++q){
                            ans.resolved4_shapes_PU->fill(weight,
                                    next.RL_res4_idx[q],
                                    next.shape_res4_idx[q],
                                    next.r_res4_idx[q],
                                    next.ct_res4_idx[q]);
                        }
                    }
                }
                //printf("\thandled res4\n");
                //fflush(stdout);
            }
            
            //TODO: res4fixed
            if constexpr (doRes4Fixed && order == 4){

            }
        }
        
        //TODO: transfer
        
        if constexpr (order < maxOrder){
            //printf("\tpunting...\n");
            //fflush(stdout);
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
    }
};

#endif
