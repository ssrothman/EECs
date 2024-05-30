#ifndef EECS_FAST1_H
#define EECS_FAST1_H

#include "fast2.h"

namespace fastEEC{
    template <typename T, bool doPU, bool doTransfer, unsigned maxOrder, bool doRes3, bool doRes4, bool doRes4Fixed>
    void do1(result<T>& ans,
             const jetDetails_t<T>& J,
             unsigned nPart,

             const res3axes_t& res3ax,
             const res4shapesAxes_t& res4ax,
             const res4fixedAxes_t& res4fixedax,

             const transferInputs<T>& tin,
             const vector<bool>* const PU){

        for(unsigned i0=0; i0 < nPart; ++i0){
            T partial0 = J.Es[i0];

            bool isPU=false;
            if constexpr(doPU){
                isPU = PU->at(i0);
            }

            if constexpr (doTransfer){
                const uvec& adj0 = tin.adj.at(i0);
                if (adj0.empty()){
                    do2<T, doPU, doTransfer, maxOrder, true, doRes3, doRes4, doRes4Fixed>(
                            ans, J, nPart,
                            res3ax, res4ax, res4fixedax,
                            tin, PU,

                            i0, partial0, isPU,

                            0, 0 //transfer information
                    );
                } else {
                    unsigned j0 = adj0[0];
                    T partialtrans0 = tin.ptrans[i0][j0];
                    do2<T, doPU, doTransfer, maxOrder, true, doRes3, doRes4, doRes4Fixed>(
                        ans, J, nPart,
                        res3ax, res4ax, res4fixedax,
                        tin, PU,

                        i0, partial0, isPU,

                        j0, partialtrans0,
                    );
                    for(unsigned j=1; j<adj0.size(); ++j){
                        j0 = adj0[j];
                        partialtrans0 = tin.ptrans[i0][j0];
                        do2<T, doPU, doTransfer, maxOrder, false, doRes3, doRes4, doRes4Fixed>(
                            ans, J, nPart,
                            res3ax, res4ax, res4fixedax,
                            tin, PU,

                            i0, partial0, isPU,

                            j0, partialtrans0,
                        );
                    }
                }
            } else {
                do2<T, doPU, doTransfer, maxOrder, true, doRes3, doRes4, doRes4Fixed>(
                    ans, J, nPart,
                    res3ax, res4ax, res4fixedax,
                    tin, PU,

                    i0, partial0, isPU,

                    0, 0,
                );
            }
        }
    }
};

#endif
