#ifndef EECS_FASTSTART_H
#define EECS_FASTSTART_H

#include "clear.h"
#include "fast1.h"

namespace fastEEC{
    template <typename T, bool doPU, bool doTransfer, unsigned maxOrder, bool doRes3, bool doRes4, bool doRes4Fixed>
    void start(const umat& dRs, 
               const vector<T>& Es,
               const unsigned NDR,
               const resolvedInputs<T>& rin,
               result<T>& ans,
               const vector<bool> *const PU = nullptr,
               const transferInputs<T>* const tin = nullptr) {
        static_assert(maxOrder >= 2 && maxOrder <= 6);

        printf("top of start\n");
        fflush(stdout);

        clear<T, doPU, doTransfer, maxOrder, doRes3, doRes4, doRes4Fixed>(ans, NDR, rin);
        printf("after clear\n");
        fflush(stdout);

        unsigned nPart = Es.size();

        printf("about to do1\n");
        fflush(stdout);
        do1<T, doPU, doTransfer, maxOrder, doRes3, doRes4, doRes4Fixed>(
            dRs, Es, nPart, rin, ans, PU, tin
        );
        printf("done\n");
        fflush(stdout);
    }
};

#endif
