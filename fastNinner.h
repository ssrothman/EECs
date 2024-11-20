#ifndef EECS_FASTN_INNER_H
#define EECS_FASTN_INNER_H

#include "util.h"
#include "fastStructs.h"
#include "symFact.h"
#include "prev.h"
#include "res3.h"
#include "res4.h"
#include <assert.h>
#include "fastNtransfer.h"

namespace fastEEC{
    template <typename T,                    //result type
             bool doPU, bool doTransfer,     //flags
             bool doRes3, bool doRes4,       //more flags
             unsigned maxOrder,              //max EEC order
             unsigned order,                 //current EEC order
             unsigned symfacIndex_prev       //symfac index 
    >
    void doN(result_t<T>& ans,
             const jetDetails_t<T>& jetDetails,
             const unsigned nPart,

             const res3axes_t& res3ax,
             const res4shapesAxes_t& res4ax,

             const transferInputs<T>& tin,
             const vector<bool>* const PU,

             const prev_t<T, order>& prev) noexcept;
    

    template <typename T,                         //result type
             bool doPU, bool doTransfer,          //flags
             bool doRes3, bool doRes4,            //more flags
             unsigned maxOrder,                   //max EEC order
             unsigned order,                      //current EEC order
             unsigned symfacIndex                 //symfac index (passthrough)
    >
    void doNinner(result_t<T>& ans,
                  const unsigned inew,
                  
                  const jetDetails_t<T>& jetDetails,
                  const unsigned nPart,

                  const res3axes_t& res3ax,
                  const res4shapesAxes_t& res4ax,

                  const transferInputs<T>& tin,
                  const vector<bool>* const PU,

                  const prev_t<T, order>& prev,
                  const T symfac) noexcept {

        prev_t<T, order+1> next;

        if constexpr(order == 1){
            next.partial = jetDetails.Es[inew];
        } else {
            next.partial = prev.partial * jetDetails.Es[inew];
        }

        if constexpr (order == 1){
            next.maxDR = 0;
            next.maxDRbin = 0;
        } else {
            std::array<T, order> dRlist;
            std::array<unsigned, order> dRbin_list;
            for(unsigned iiold =0; iiold<order-1; ++iiold){
                const unsigned& iold = prev.is[iiold];
                dRlist[iiold] = jetDetails.floatDRs[inew][iold];
                dRbin_list[iiold] = jetDetails.dRbins[inew][iold];
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

        if constexpr (doPU){
            next.isPU = prev.isPU || PU->at(inew);
        }

        if constexpr(order > 1){
            T weight = symfac * next.partial;

            for (unsigned i=0; i<order-2; ++i){
                next.wts[i] = prev.wts[i];
                next.dRbins[i] = prev.dRbins[i];
            }
            next.wts[order-2] = weight;
            next.dRbins[order-2] = next.maxDRbin;

            //accumulate
            (*ans.wts[order-2])[next.maxDRbin] += weight;

            if constexpr (doPU){
                if(next.isPU){
                    (*ans.wts_PU[order-2])[next.maxDRbin] += weight;
                }//end if isPU
            }//end if doPU

            if constexpr (doRes3 && order == 3){
                runRes3<T>(jetDetails, res3ax, next);
                (*ans.resolved3)[next.RL_res3_idx][next.xi_res3_idx][next.phi_res3_idx] += weight;
                if constexpr (doPU){
                    if(next.isPU){
                        (*ans.resolved3_PU)[next.RL_res3_idx][next.xi_res3_idx][next.phi_res3_idx] += weight;
                    }
                }
            } else if constexpr(doRes3 && order > 3){
                next.RL_res3_idx = prev.RL_res3_idx;
                next.xi_res3_idx = prev.xi_res3_idx;
                next.phi_res3_idx = prev.phi_res3_idx;
            }

            if constexpr (doRes4 && order == 4){
                runRes4<T>(jetDetails, res4ax, next);
                for(unsigned q=0; q<3; ++q){
                    ans.resolved4_shapes->fill(weight,
                            next.RL_res4_idx[q],
                            next.shape_res4_idx[q],
                            next.r_res4_idx[q],
                            next.ct_res4_idx[q]);
                }
                for (unsigned q=0; q<4; ++q){
                    ans.resolved4_shapes->fillTri(weight,
                            next.istri_res4[q],
                            next.RL_res4tri_idx[q],
                            next.r_res4tri_idx[q],
                            next.ct_res4tri_idx[q]);
                }

                if constexpr(doPU){
                    if(next.isPU){
                        for(unsigned q=0; q<3; ++q){
                            ans.resolved4_shapes_PU->fill(weight,
                                    next.RL_res4_idx[q],
                                    next.shape_res4_idx[q],
                                    next.r_res4_idx[q],
                                    next.ct_res4_idx[q]);
                        }
                        for (unsigned q=0; q<4; ++q){
                            ans.resolved4_shapes_PU->fillTri(weight,
                                    next.istri_res4[q],
                                    next.RL_res4tri_idx[q],
                                    next.r_res4tri_idx[q],
                                    next.ct_res4tri_idx[q]);
                        }
                    }
                }
            } else if constexpr(doRes4 && order > 4){
                for(unsigned q=0; q<3; ++q){
                    next.RL_res4_idx[q] = prev.RL_res4_idx[q];
                    next.shape_res4_idx[q] = prev.shape_res4_idx[q];
                    next.r_res4_idx[q] = prev.r_res4_idx[q];
                    next.ct_res4_idx[q] = prev.ct_res4_idx[q];

                    next.istri_res4[q] = prev.istri_res4[q];
                    next.RL_res4tri_idx[q] = prev.RL_res4tri_idx[q];
                    next.r_res4tri_idx[q] = prev.r_res4tri_idx[q];
                    next.ct_res4tri_idx[q] = prev.ct_res4tri_idx[q];
                }
            }
        }
        
        if constexpr (order < maxOrder){
            doN<T, 
                doPU, doTransfer, 
                doRes3, doRes4, 
                maxOrder, order+1,
                symfacIndex
            >(
                    ans, 
                    jetDetails, nPart, 
                    res3ax, res4ax, 
                    tin, PU, 
                    next
            );
        } else if constexpr (doTransfer){
            prev_t<T, 1> prevTrans;
            transferN<T, 
                doPU, doRes3, doRes4, 
                maxOrder, 1,
                symfacIndex,
                0>(
                        ans, 
                        jetDetails, nPart, 
                        res3ax, res4ax, 
                        tin, 
                        next, 
                        prevTrans
            );
        }
    }
};

#endif
