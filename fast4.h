#ifndef EECS_FAST4_H
#define EECS_FAST4_H

#include "fast5.h"

namespace fastEEC{
    template <typename T, bool doPU, bool doTransfer, unsigned maxOrder, bool nontransfer>
    void do4(const umat& dRs,
             const vector<T>& Es,
             const unsigned nPart,
             result<T>& ans,

             const unsigned i0,
             const unsigned i1,
             const unsigned i2,
             const T partial2,
             const unsigned DR2,
             const unsigned i0max, 
             const unsigned i1max,
             bool isPU2,

             const unsigned j0,
             const unsigned j1,
             const unsigned j2,
             const T partialtrans2,
             const unsigned DR2_Reco,
             const unsigned j0max,
             const unsigned j1max,

             const vector<bool>* const PU = nullptr,
             const transferInputs<T>* const tin = nullptr) {

        T weight4;
        bool isPU3=isPU2;
        unsigned i0max_new=0, i1max_new=0;
        unsigned j0max_new=0, j1max_new=0;
        for(unsigned i3=i2; i3<nPart; ++i3){
            T partial3 = partial2 * Es[i3];

            uvec dRlist = {DR2, dRs[i0][i3], 
                                dRs[i1][i3], 
                                dRs[i2][i3]};
            auto maxel = max_element(dRlist.begin(), dRlist.end());
            unsigned DR3 = *maxel;
            switch(std::distance(dRlist.begin(), maxel)){
                case 0:
                    i0max_new = i0max;
                    i1max_new = i1max;
                    break;
                case 1:
                    i0max_new = 0;
                    i1max_new = 3;
                    break;
                case 2:
                    i0max_new = 1;
                    i1max_new = 3;
                    break;
                case 3:
                    i0max_new = 2;
                    i1max_new = 3;
                    break;
            };

            T symfac;
            if (i0==i1){
                if(i1==i2){
                    symfac = (i2==i3) ? 1 : 4; //(4) vs (3, 1)
                } else {
                    symfac = (i2==i3) ? 6 : 12; //(2, 2) vs (2, 1, 1)
                }
            } else {
                if(i1==i2){
                    symfac = (i2==i3) ? 4 : 12; //(1, 3) vs (1, 2, 1)
                } else {
                    symfac = (i2==i3) ? 12 : 24; //(1, 1, 2) vs (1, 1, 1, 1)
                }
            }

            weight4 = symfac * partial3;

            if constexpr(nontransfer){
                //accumulate
                ans.wts4[DR3] += weight4;
                if constexpr(doPU){
                    isPU3 = isPU2 || PU->at(i3);
                    if(isPU3){
                        ans.wts4_PU[DR3] += weight4;
                    }
                }
            }

            if constexpr(doTransfer){
                const uvec& adj3 = tin->adj.at(i3);
                if(adj3.empty() || partialtrans2==0){
                    if constexpr (maxOrder >=5 && nontransfer){
                        do5<T, doPU, doTransfer, maxOrder, nontransfer>(
                            dRs, Es, nPart, ans,
                            i0, i1, i2, i3, partial3, 
                            DR3, i0max, i1max,
                            isPU3,
                            0, 0, 0, 0, 0, 
                            0, 0, 0,
                            PU, tin
                        );
                    }
                } else {
                    unsigned j3 = adj3[0];
                    T partialtrans3 = partialtrans2 * tin->ptrans[i3][j3];

                    uvec dRlist_reco = {DR2_Reco, tin->dRs[j0][j3],
                                                  tin->dRs[j1][j3], 
                                                  tin->dRs[j2][j3]};
                    auto maxel_reco = max_element(dRlist_reco.begin(), dRlist_reco.end());
                    unsigned DR3_Reco = *maxel_reco;
                    switch(std::distance(dRlist_reco.begin(), maxel_reco)){
                        case 0:
                            j0max_new = j0max;
                            j1max_new = j1max;
                            break;
                        case 1:
                            j0max_new = 0;
                            j1max_new = 3;
                            break;
                        case 2:
                            j0max_new = 1;
                            j1max_new = 3;
                            break;
                        case 3:
                            j0max_new = 2;
                            j1max_new = 3;
                            break;
                    };
                    ans.transfer4[DR3][DR3_Reco] += partialtrans3 * weight4;

                    if constexpr (maxOrder >=5){
                        do5<T, doPU, doTransfer, maxOrder, nontransfer>(
                            dRs, Es, nPart, ans,
                            i0, i1, i2, i3, partial3, 
                            DR3, i0max_new, i1max_new,
                            isPU3,
                            j0, j1, j2, j3, partialtrans3, 
                            DR3_Reco, j0max_new, j1max_new,
                            PU, tin
                        );
                    }

                    for(unsigned j=1; j<adj3.size(); ++j){
                        j3 = adj3[j];
                        partialtrans3 = partialtrans2 * tin->ptrans[i3][j3];

                        dRlist_reco[1] = tin->dRs[j0][j3];
                        dRlist_reco[2] = tin->dRs[j1][j3];
                        dRlist_reco[3] = tin->dRs[j2][j3];
                        maxel_reco = max_element(dRlist_reco.begin(), dRlist_reco.end());
                        DR3_Reco = *maxel_reco;
                        switch(std::distance(dRlist_reco.begin(), maxel_reco)){
                            case 0:
                                j0max_new = j0max;
                                j1max_new = j1max;
                                break;
                            case 1:
                                j0max_new = 0;
                                j1max_new = 3;
                                break;
                            case 2:
                                j0max_new = 1;
                                j1max_new = 3;
                                break;
                            case 3:
                                j0max_new = 2;
                                j1max_new = 3;
                                break;
                        };

                        ans.transfer4[DR3][DR3_Reco] += partialtrans3 * weight4;

                        if constexpr (maxOrder >= 5){
                            do5<T, doPU, doTransfer, maxOrder, false>(
                                dRs, Es, nPart, ans,
                                i0, i1, i2, i3, partial3, 
                                DR3, i0max_new, i1max_new,
                                isPU3,
                                j0, j1, j2, j3, partialtrans3, 
                                DR3_Reco, j0max_new, j1max_new,
                                PU, tin
                            );
                        }
                    }
                }
            } else {
                if constexpr(maxOrder >= 5){
                    do5<T, doPU, doTransfer, maxOrder, true>(
                        dRs, Es, nPart, ans,
                        i0, i1, i2, i3, partial3, 
                        DR3, i0max_new, i1max_new,
                        isPU3,
                        0, 0, 0, 0, 0, 
                        0, 0, 0,
                        PU, tin
                    );
                }
            }
        }
        if(isPU3){
            isPU3 = true;
        } else if(i0max_new){
            i0max_new = 0;
        } else if(i1max_new){
            i1max_new = 0;
        } else if(j0max_new){
            j0max_new = 0;
        } else if(j1max_new){
            j1max_new = 0;
        }
    }
};

#endif
