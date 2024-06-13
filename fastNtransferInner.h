#ifndef FASTN_TRANSFER_INNER_H
#define FASTN_TRANSFER_INNER_H

#include "util.h"
#include "fastStructs.h"
#include "symFact.h"
#include "prev.h"
#include "res3.h"
#include "res4.h"
#include <assert.h>

namespace fastEEC{
    template <typename T,
             bool doPU, 
             bool doRes3, bool doRes4,
             unsigned maxOrder,
             unsigned order,
             unsigned symfacIndex,
             bool accumulate_prev
    >
    void transferN(result_t<T>& ans,
                   const jetDetails_t<T>& jetDetails,
                   const unsigned nPart,

                   const res3axes_t& res3ax,
                   const res4shapesAxes_t& res4ax,

                   const transferInputs<T>& tin,

                   const prev_t<T, maxOrder+1>& prevGen,

                   const prev_t<T, order>& prevTrans);

    template <typename T,
             bool doPU,
             bool doRes3, bool doRes4,
             unsigned maxOrder,
             unsigned order,
             unsigned symfacIndex,
             bool accumulate_new
    >
    void transferNinner(result_t<T>& ans,
                        const jetDetails_t<T>& jetDetails,
                        const unsigned nPart,

                        const res3axes_t& res3ax,
                        const res4shapesAxes_t& res4ax,

                        const transferInputs<T>& tin,

                        const prev_t<T, maxOrder+1>& prevGen,

                        const prev_t<T, order>& prevTrans,

                        const unsigned jnew){
        const unsigned inew = prevGen.is[order-1];

        prev_t<T, order+1> nextTrans;

        for (unsigned o=0; o<order-1; ++o){
            nextTrans.is[o] = prevTrans.is[o];
        }
        nextTrans.is[order-1] = jnew;

        //printf("TOP OF TRANSFER INNER for (");
        //printVec(prevGen.is);
        //printf(")- > (");
        //printVec(nextTrans.is);
        //printf(")\n");
        //fflush(stdout);
        if constexpr(order == 1){
            nextTrans.partial = (*tin.ptrans)(inew, jnew);
        } else {
            nextTrans.partial = prevTrans.partial * (*tin.ptrans)(inew, jnew);
        }
        //printf("\tpartial = %f\n", nextTrans.partial);
        //fflush(stdout);
        /*printf("Got partial\n");
        fflush(stdout);*/

        if constexpr(order == 1){
            nextTrans.maxDR = 0;
            nextTrans.maxDRbin = 0;
        } else {
            std::array<T, order> dRlist;
            std::array<unsigned, order> dRbin_list;
            for(unsigned jold=0; jold<order-1; ++jold){
                dRlist[jold] = tin.recoJet->floatDRs[prevTrans.is[jold]][jnew];
                dRbin_list[jold] = tin.recoJet->dRbins[prevTrans.is[jold]][jnew];
            }
            dRlist[order-1] = prevTrans.maxDR;
            dRbin_list[order-1] = prevTrans.maxDRbin;

            auto maxel = max_element(dRlist.begin(), dRlist.end());
            nextTrans.maxDR = *maxel;
            nextTrans.maxDRbin = dRbin_list[std::distance(dRlist.begin(), maxel)];
        }
        //printf("\tmaxDR = %f\n", nextTrans.maxDR);
        //printf("\tmaxDRbin = %u\n", nextTrans.maxDRbin);
        //fflush(stdout);
        /*printf("GOT MAX DR\n");
        fflush(stdout);*/


        if constexpr(accumulate_new && order > 1){
            T weight = nextTrans.partial * prevGen.wts[order-2];

            /*printf("DR transfer is %u -> %u\n", prevGen.dRbins[order-2], 
                    nextTrans.maxDRbin);*/
            /*printf("\taccumulating proj %u->%u += %g\n", 
                    prevGen.dRbins[order-2], 
                    nextTrans.maxDRbin,
                    weight);*/
            //fflush(stdout);
            (*ans.transfer_wts[order-2])[prevGen.dRbins[order-2]][nextTrans.maxDRbin] += weight;

            if constexpr(doRes3 && order == 3){
                runRes3<T>(*tin.recoJet, res3ax, nextTrans);

                /*printf("res3 transfer is (%u, %u, %u) -> (%u, %u, %u)\n",
                        prevGen.RL_res3_idx, 
                        prevGen.xi_res3_idx,
                        prevGen.phi_res3_idx,
                        nextTrans.RL_res3_idx,
                        nextTrans.xi_res3_idx,
                        nextTrans.phi_res3_idx);*/
                /*printf("\tacc res3 (%u, %u, %u) -> (%u, %u, %u) += %g\n",
                        prevGen.RL_res3_idx, 
                        prevGen.xi_res3_idx,
                        prevGen.phi_res3_idx,
                        nextTrans.RL_res3_idx,
                        nextTrans.xi_res3_idx,
                        nextTrans.phi_res3_idx,
                        weight);*/
                (*ans.transfer_res3)[prevGen.RL_res3_idx][prevGen.xi_res3_idx][prevGen.phi_res3_idx][nextTrans.RL_res3_idx][nextTrans.xi_res3_idx][nextTrans.phi_res3_idx] += weight;
            }

            if constexpr(doRes4 && order == 4){
                runRes4<T>(*tin.recoJet, res4ax, nextTrans);

                for (unsigned q=0; q<3; ++q){
                    /*if (prevGen.shape_res4_idx[q] == 0){
                        printf("\tacc res4 noshape -> shape %u\n", 
                                nextTrans.shape_res4_idx[q]);
                    } else if (prevGen.shape_res4_idx[q] == 1){
                        printf("\tacc res4 dipole -> shape %u\n", 
                                nextTrans.shape_res4_idx[q]);
                        printf("\t\twith (%u, %u, %u) -> (%u, %u, %u)\n",
                                prevGen.RL_res4_idx[q],
                                prevGen.r_res4_idx[q],
                                prevGen.ct_res4_idx[q],
                                nextTrans.RL_res4_idx[q],
                                nextTrans.r_res4_idx[q],
                                nextTrans.ct_res4_idx[q]);
                    } else {
                        printf("\tacc res4 tee -> shape %u\n",
                                nextTrans.shape_res4_idx[q]);
                        printf("\t\twith (%u, %u, %u) -> (%u, %u, %u)\n",
                                prevGen.RL_res4_idx[q],
                                prevGen.r_res4_idx[q],
                                prevGen.ct_res4_idx[q],
                                nextTrans.RL_res4_idx[q],
                                nextTrans.r_res4_idx[q],
                                nextTrans.ct_res4_idx[q]);
                    }*/

                    ans.transfer_res4_shapes->fill(weight,
                            prevGen.RL_res4_idx[q],
                            prevGen.shape_res4_idx[q],
                            prevGen.r_res4_idx[q],
                            prevGen.ct_res4_idx[q],
                            nextTrans.RL_res4_idx[q],
                            nextTrans.shape_res4_idx[q],
                            nextTrans.r_res4_idx[q],
                            nextTrans.ct_res4_idx[q]);
                }
            }
        }
        //printf("\tdone.\n");

        if constexpr(order < maxOrder){
            transferN<T, 
                      doPU, doRes3, doRes4,
                      maxOrder, order+1,
         symfacIndex,
                      accumulate_new>(
                                ans,
                                jetDetails, nPart,
                                res3ax, res4ax,
                                tin, prevGen, 
                                nextTrans);
        }
    }
}

#endif
