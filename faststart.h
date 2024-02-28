#ifndef EECS_FASTSTART_H
#define EECS_FASTSTART_H

#include "clear.h"
#include "fast1.h"

namespace fastEEC{
    template <typename T, bool doPU, bool doTransfer, unsigned maxOrder>
    void start(const umat& dRs, 
               const vector<T>& Es,
               const unsigned NDR,
               const resolvedInputs<T>& rin,
               result<T>& ans,
               const vector<bool> *const PU = nullptr,
               const transferInputs<T>* const tin = nullptr) {
        static_assert(maxOrder >= 2 && maxOrder <= 6);

        clear<T, doPU, doTransfer, maxOrder>(ans, NDR, rin);

        printf("running fast EEC\n");
        printf("\tdoPU = %d\n", doPU);
        printf("\tdoTransfer = %d\n", doTransfer);
        printf("\tmaxOrder = %d\n", maxOrder);

        unsigned nPart = Es.size();

        do1<T, doPU, doTransfer, maxOrder>(
            dRs, Es, nPart, rin, ans, PU, tin
        );
    }
};

#endif
