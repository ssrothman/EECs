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
             bool doRes3, bool doRes4,                    //more flags
             unsigned maxOrder,                           //max EEC order
             unsigned order,                              //current EEC order
             unsigned symfacIndex_prev                    //symfac index
    >
    void doN(result_t<T>& ans,
             const jetDetails_t<T>& jetDetails,
             const unsigned nPart,

             const res3axes_t& res3ax,
             const res4shapesAxes_t& res4ax,

             const transferInputs<T>& tin,
             const vector<bool>* const PU,

             const prev_t<T, order>& prev) noexcept {

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
        unsigned istart = 0;
        if constexpr (order > 1){
            istart = prev.is[order-2];
        }

        constexpr unsigned symfacIndex_new_first = (order < 2) ? 0 : (symfacIndex_prev << 1);
        constexpr T symfac_first = (order < 2) ? 1 : symfacLookup<T, order>()[symfacIndex_new_first];

        doNinner<T, doPU, doTransfer,
            doRes3, doRes4, 
            maxOrder, order,
            symfacIndex_new_first>(
                    ans, istart,
                    jetDetails, nPart, 
                    res3ax, res4ax, 
                    tin, PU, 
                    prev, symfac_first);

        constexpr unsigned symfacIndex_new_all = (order < 2) ? symfacIndex_new_first : symfacIndex_new_first + 1;
        constexpr T symfac_all = (order < 2) ? 1 : symfacLookup<T, order>()[symfacIndex_new_all];

        for(unsigned inew = istart+1; inew < nPart; ++inew){
            doNinner<T, doPU, doTransfer, 
                doRes3, doRes4, 
                maxOrder, order,
                symfacIndex_new_all>(
                        ans, inew, 
                        jetDetails, nPart, 
                        res3ax, res4ax, 
                        tin, PU, 
                        prev, symfac_all);
        }//end loop
    }//end function
}

#endif
