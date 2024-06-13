#ifndef EECS_FASTNTRANSFER_H
#define EECS_FASTNTRANSFER_H

#include "util.h"
#include "fastStructs.h"
#include "symFact.h"
#include "prev.h"
#include "res3.h"
#include "res4.h"
#include "fastNtransferInner.h"
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

                   const prev_t<T, order>& prevTrans){


        //printf("TOP OF TRANSFERN\n");
        //fflush(stdout);

        const std::vector<unsigned>& jnews = tin.adj->at(prevGen.is[order-1]);
        if(jnews.empty()){//if there are no neighbors, we can break early
            return;
        }

        //printf("TRANSFERN: %d\n", order);
        //fflush(stdout);

        constexpr unsigned ander = (1 << (maxOrder - order + 1)) - 1;
        constexpr unsigned comparator = 1 << (maxOrder - order);
        constexpr bool accumulate_new = accumulate_prev | 
                ((order == 1) ? symfacIndex == 0 : 
                               ((symfacIndex&ander)==comparator));

        //printf("accumulate_new = %d\n", accumulate_new);
        //fflush(stdout);

        for (const unsigned& jnew : jnews){
            //printf("jnew = %d\n", jnew);
            //fflush(stdout);
            transferNinner<T, 
                           doPU, 
                           doRes3, doRes4, 
                           maxOrder, order, 
                           symfacIndex, accumulate_new>(
                ans, jetDetails, nPart,
                res3ax, res4ax,
                tin,
                prevGen, prevTrans, jnew
            );
        }
    }
}

#endif
