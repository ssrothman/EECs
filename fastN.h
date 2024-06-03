#ifndef EECS_FASTN_H
#define EECS_FASTN_H

#include "util.h"
#include "fastStructs.h"
#include "symFact.h"
#include "prev.h"
#include "res3.h"
#include "res4.h"

#include "fastNinner.h"

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

             const prev_t<T, order>& prev){

        /*
         * I think we can do something cute with symfacs
         *
         * The symfac in the first iteration is something cute related to the symfac from the above level of nesting
         * The symfac in all subsequent iterations is identical, and can be taken as a constant (based on the prev_symfac)
         *
         * I think the thing to do will be to factorize the inside of the loop into a different function
         * to which you pass the symfac
         *
         * So then you can call it once with the first-iteration symfac,
         * and then  N-1 times with the all-other-iterations symfac
         *
         * The symfac lookup can be in a constexpr lookup table
         * which I think we should specify as an array (constant-time lookups)
         * where the index is a bitmask for each of the comparisons
         *
         * The new bitmask is just the prev one << 1 + 1
         * and then the all-other-iterations one is just that -1
         *
         * easy easy
         */
        //printf("top of doN\n");
        //fflush(stdout);

        unsigned istart = 0;
        if constexpr (order > 1){
            istart = prev.is[order-2];
        }

        T symfac;
        unsigned symfacIndex = 0;
        if constexpr (order < 2){
            symfac = 1;
        } else {
            symfacIndex = (prev.symfacIndex << 1);
            symfac = symfacLookup<T, order>()[symfacIndex];
        }
        //printVec(prev.is);
        //printf("%u\n", istart);
        //printf("\tidx: %u\n", symfacIndex);
        //printf("\tsymfac: %f\n", symfac);
        //fflush(stdout);

        doNinner<T, doPU, doTransfer,
            doRes3, doRes4, doRes4Fixed,
            maxOrder, order>(
                    ans, istart,
                    jetDetails, nPart, 
                    res3ax, res4ax, res4fixedax, 
                    tin, PU, 
                    prev, symfac, symfacIndex);

        if constexpr (order < 2){
            symfac = 1;
        } else {
            symfacIndex += 1;
            symfac = symfacLookup<T, order>()[symfacIndex];
        }

        for(unsigned inew = istart+1; inew < nPart; ++inew){
            //printVec(prev.is);
            //printf("%u\n", inew);
            //printf("\tidx: %u\n", symfacIndex);
            //printf("\tsymfac: %f\n", symfac);
            //fflush(stdout);

            doNinner<T, doPU, doTransfer, 
                doRes3, doRes4, doRes4Fixed, 
                maxOrder, order>(
                        ans, inew, 
                        jetDetails, nPart, 
                        res3ax, res4ax, res4fixedax, 
                        tin, PU, 
                        prev, symfac, symfacIndex);
        }//end loop
    }//end function
}

#endif
