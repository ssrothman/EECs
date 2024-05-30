#ifndef EECS_FASTSTART_H
#define EECS_FASTSTART_H

#include "clear.h"
#include "fast1.h"

namespace fastEEC{
    template <typename T, bool doPU, bool doTransfer, unsigned maxOrder, bool doRes3, bool doRes4, bool doRes4Fixed>
    void start(result<T>& ans,

               const jetDetails_t& J,

               const unsigned NDR,
               const res3axes_t& res3ax,
               const res4shapesAxes_t& res4ax,
               const res4fixedAxes_t& res4fixedax,

               const transferInputs<T>& tin,

               const vector<bool> *const PU = nullptr){
        static_assert(maxOrder >= 2 && maxOrder <= 6);

        clear<T, doPU, doTransfer, maxOrder, doRes3, doRes4, doRes4Fixed>(ans, NDR, res3ax, res4ax, res4fixedax);

        unsigned nPart = J.Es.size();

        do1<T, doPU, doTransfer, maxOrder, doRes3, doRes4, doRes4Fixed>(
                ans, J, nPart, res3ax, res4ax, res4fixedax, tin, PU
        );
    }
};

#endif
