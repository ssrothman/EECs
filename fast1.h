#ifndef EECS_FAST1_H
#define EECS_FAST1_H

#include "fast2.h"

namespace fastEEC{
    template <typename T, bool doPU, bool doTransfer, unsigned maxOrder>
    void do1(const umat& dRs, 
             const vector<T>& Es,
             const unsigned nPart,
             const resolvedInputs<T>& rin,
             result<T>& ans,
             const vector<bool>* const PU = nullptr,
             const transferInputs<T>* const tin = nullptr) {

        for(unsigned i0=0; i0 < nPart; ++i0){
            T partial0 = Es[i0];

            bool isPU=false;
            if constexpr(doPU){
                isPU = PU->at(i0);
            }

            if constexpr (doTransfer){
                const uvec& adj0 = tin->adj.at(i0);
                if (adj0.empty()){
                    do2<T, doPU, doTransfer, maxOrder, true>(
                        dRs, Es, nPart, rin, ans,
                        i0, partial0, isPU,
                        0, 0,
                        PU, tin
                    );
                } else {
                    unsigned j0 = adj0[0];
                    T partialtrans0 = tin->ptrans[i0][j0];
                    do2<T, doPU, doTransfer, maxOrder, true>(
                        dRs, Es, nPart, rin, ans,
                        i0, partial0, isPU,
                        j0, partialtrans0,
                        PU, tin
                    );
                    for(unsigned j=1; j<adj0.size(); ++j){
                        j0 = adj0[j];
                        partialtrans0 = tin->ptrans[i0][j0];
                        do2<T, doPU, doTransfer, maxOrder, false>(
                            dRs, Es, nPart, rin, ans,
                            i0, partial0, isPU,
                            j0, partialtrans0,
                            PU, tin
                        );
                    }
                }
            } else {
                do2<T, doPU, doTransfer, maxOrder, true>(
                    dRs, Es, nPart, rin, ans,
                    i0, partial0, isPU,
                    0, 0,
                    PU, tin
                );
            }
        }
    }
};

#endif
