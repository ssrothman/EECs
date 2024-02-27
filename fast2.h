#ifndef EECS_FAST2_H
#define EECS_FAST2_H

#include "fast3.h"

namespace fastEEC{
    template <typename T, bool doPU, bool doTransfer, unsigned maxOrder, bool nontransfer>
    void do2(const umat& dRs,
             const vector<T>& Es,
             const unsigned nPart,
             result<T>& ans,

             const unsigned i0,
             const T partial0,
             bool isPU0, 

             const unsigned j0,
             const T partialtrans0,

             const vector<bool>* const PU = nullptr,
             const transferInputs<T>* const tin = nullptr) {

        T weight2;
        bool isPU1=isPU0;
        for(unsigned i1=i0; i1<nPart; ++i1){
            T partial1 = partial0 * Es[i1];
            unsigned DR1 = dRs[i0][i1];

            T symfac = (i0==i1) ? 1 : 2;
            weight2 = symfac * partial1;

            if constexpr(nontransfer){
                //accumulate
                ans.wts2[DR1] += weight2;
                if constexpr(doPU){
                    isPU1 = isPU0 || PU->at(i1);
                    if(isPU1){
                        ans.wts2_PU[DR1] += weight2;
                    }
                }
            }

            if constexpr(doTransfer){
                const uvec& adj1 = tin->adj.at(i1);
                if(adj1.empty() || partialtrans0 == 0){
                    if constexpr (maxOrder >=3 && nontransfer){
                        do3<T, doPU, doTransfer, maxOrder, nontransfer>(
                            dRs, Es, nPart, ans,
                            i0, i1, partial1, DR1, isPU1,
                            0, 0, 0, 0,
                            PU, tin
                        );
                    }
                } else {
                    unsigned j1 = adj1[0];
                    T partialtrans1 = partialtrans0 * tin->ptrans[i1][j1];
                    unsigned DR1_Reco = tin->dRs[j0][j1];
                    ans.transfer2[DR1][DR1_Reco] += partialtrans1 * weight2;

                    if constexpr (maxOrder >=3){
                        do3<T, doPU, doTransfer, maxOrder, nontransfer>(
                            dRs, Es, nPart, ans,
                            i0, i1, partial1, DR1, isPU1,
                            j0, j1, partialtrans1, DR1_Reco,
                            PU, tin
                        );
                    }

                    for(unsigned j=1; j<adj1.size(); ++j){
                        j1 = adj1[j];
                        partialtrans1 = partialtrans0 * tin->ptrans[i1][j1];
                        unsigned DR1_Reco = tin->dRs[j0][j1];
                        ans.transfer2[DR1][DR1_Reco] += partialtrans1 * weight2;

                        if constexpr (maxOrder >= 3){
                            do3<T, doPU, doTransfer, maxOrder, false>(
                                dRs, Es, nPart, ans,
                                i0, i1, partial1, DR1, isPU1,
                                j0, j1, partialtrans1, DR1_Reco,
                                PU, tin
                            );
                        }
                    }
                }
            } else {
                if constexpr(maxOrder >= 3){
                    do3<T, doPU, doTransfer, maxOrder, true>(
                        dRs, Es, nPart, ans,
                        i0, i1, partial1, DR1, isPU1,
                        0, 0, 0, 0,
                        PU, tin
                    );
                }
            }
        }
        if(isPU1){
            isPU1 = true;
        }
    }
}

#endif
