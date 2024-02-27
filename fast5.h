#ifndef EECS_FAST5_H
#define EECS_FAST5_H

#include "fast6.h"

namespace fastEEC{
    template <typename T, bool doPU, bool doTransfer, unsigned maxOrder, bool nontransfer>
    void do5(const umat& dRs,
             const vector<T>& Es,
             const unsigned nPart,
             result<T>& ans,

             const unsigned i0,
             const unsigned i1,
             const unsigned i2,
             const unsigned i3,
             const T partial3,
             const unsigned DR3,
             const unsigned i0max,
             const unsigned i1max,
             bool isPU3,

             const unsigned j0,
             const unsigned j1,
             const unsigned j2,
             const unsigned j3,
             const T partialtrans3,
             const unsigned DR3_Reco,
             const unsigned j0max, 
             const unsigned j1max,

             const vector<bool>* const PU = nullptr,
             const transferInputs<T>* const tin = nullptr) {

        T weight5;
        bool isPU4=isPU3;
        unsigned i0max_new=0, i1max_new=0;
        unsigned j0max_new=0, j1max_new=0;
        for(unsigned i4=i3; i4<nPart; ++i4){
            T partial4 = partial3 * Es[i4];

            uvec dRlist = {DR3, dRs[i0][i4], 
                                dRs[i1][i4], 
                                dRs[i2][i4], 
                                dRs[i3][i4]};
            auto maxel = max_element(dRlist.begin(), dRlist.end());
            unsigned DR4 = *maxel;
            switch(std::distance(dRlist.begin(), maxel)){
                case 0:
                    i0max_new = i0max;
                    i1max_new = i1max;
                    break;
                case 1:
                    i0max_new = 0;
                    i1max_new = 4;
                    break;
                case 2:
                    i0max_new = 1;
                    i1max_new = 4;
                    break;
                case 3:
                    i0max_new = 2;
                    i1max_new = 4;
                    break;
                case 4:
                    i0max_new = 3;
                    i1max_new = 4;
                    break;
            };

            T symfac;
            if (i0==i1){
                if(i1==i2){
                    if(i2==i3){
                        symfac = (i3==i4) ? 1 : 5; //(5) vs (4, 1)
                    } else {
                        symfac = (i3==i4) ? 10 : 20; //(3, 2) vs (3, 1, 1)
                    }
                } else {
                    if(i2==i3){
                        symfac = (i3==i4) ? 10 : 30; //(2, 3) vs (2, 2, 1)
                    } else {
                        symfac = (i3==i4) ? 30 : 60; //(2, 1, 2) vs (2, 1, 1, 1)
                    }
                }
            } else {
                if(i1==i2){
                    if(i2==i3){
                        symfac = (i3==i4) ? 5 : 20; //(1, 4) vs (1, 3, 1)
                    } else{
                        symfac = (i3==i4) ? 30 : 60; //(1, 2, 2) vs (1, 2, 1, 1)
                    }
                } else {
                    if(i2==i3){
                        symfac = (i3==i4) ? 20 : 60; //(1, 1, 3) vs (1, 1, 2, 1)
                    } else {
                        symfac = (i3==i4) ? 60 : 120; //(1, 1, 1, 2) vs (1, 1, 1, 1, 1)
                    }
                }
            }
            weight5 = symfac * partial4;

            if constexpr(nontransfer){
                //accumulate
                ans.wts5[DR4] += weight5;
                if constexpr(doPU){
                    isPU4 = isPU3 || PU->at(i4);
                    if(isPU4){
                        ans.wts5_PU[DR4] += weight5;
                    }
                }
            }

            if constexpr(doTransfer){
                const uvec& adj4 = tin->adj.at(i4);
                if(adj4.empty() || partialtrans3 == 0){
                    if constexpr (maxOrder >=6 && nontransfer){
                        do6<T, doPU, doTransfer, maxOrder, nontransfer>(
                            dRs, Es, nPart, ans,
                            i0, i1, i2, i3, i4, partial4, 
                            DR4, i0max_new, i1max_new,
                            isPU4,
                            0, 0, 0, 0, 0, 0,
                            0, 0, 0,
                            PU, tin
                        );
                    }
                } else {
                    unsigned j4 = adj4[0];
                    T partialtrans4 = partialtrans3 * tin->ptrans[i4][j4];

                    uvec dRlist_reco = {DR3_Reco, tin->dRs[j0][j4], 
                                                  tin->dRs[j1][j4], 
                                                  tin->dRs[j2][j4], 
                                                  tin->dRs[j3][j4]};
                    auto maxel_reco = max_element(dRlist_reco.begin(), dRlist_reco.end());
                    unsigned DR4_Reco = *maxel_reco;
                    switch(std::distance(dRlist_reco.begin(), maxel_reco)){
                        case 0:
                            j0max_new = j0max;
                            j1max_new = j1max;
                            break;
                        case 1:
                            j0max_new = 0;
                            j1max_new = 4;
                            break;
                        case 2:
                            j0max_new = 1;
                            j1max_new = 4;
                            break;
                        case 3:
                            j0max_new = 2;
                            j1max_new = 4;
                            break;
                        case 4:
                            j0max_new = 3;
                            j1max_new = 4;
                            break;
                    };

                    ans.transfer5[DR4][DR4_Reco] += partialtrans4 * weight5;

                    if constexpr (maxOrder >=6){
                        do6<T, doPU, doTransfer, maxOrder, nontransfer>(
                            dRs, Es, nPart, ans,
                            i0, i1, i2, i3, i4, partial4, 
                            DR4, i0max_new, i1max_new,
                            isPU4,
                            j0, j1, j2, j3, j4, partialtrans4, 
                            DR4_Reco, j0max_new, j1max_new,
                            PU, tin
                        );
                    }

                    for(unsigned j=1; j<adj4.size(); ++j){
                        j4 = adj4[j];
                        partialtrans4 = partialtrans3 * tin->ptrans[i4][j4];

                        dRlist_reco[1] = tin->dRs[j0][j4];
                        dRlist_reco[2] = tin->dRs[j1][j4];
                        dRlist_reco[3] = tin->dRs[j2][j4];
                        dRlist_reco[4] = tin->dRs[j3][j4];
                        maxel_reco = max_element(dRlist_reco.begin(), dRlist_reco.end());
                        unsigned DR4_Reco = *maxel_reco;
                        switch(std::distance(dRlist_reco.begin(), maxel_reco)){
                            case 0:
                                j0max_new = j0max;
                                j1max_new = j1max;
                                break;
                            case 1:
                                j0max_new = 0;
                                j1max_new = 4;
                                break;
                            case 2:
                                j0max_new = 1;
                                j1max_new = 4;
                                break;
                            case 3:
                                j0max_new = 2;
                                j1max_new = 4;
                                break;
                            case 4:
                                j0max_new = 3;
                                j1max_new = 4;
                                break;
                        };

                        ans.transfer5[DR4][DR4_Reco] += partialtrans4 * weight5;

                        if constexpr (maxOrder >= 6){
                            do6<T, doPU, doTransfer, maxOrder, false>(
                                dRs, Es, nPart, ans,
                                i0, i1, i2, i3, i4, partial4, 
                                DR4, i0max_new, i1max_new,
                                isPU4,
                                j0, j1, j2, j3, j4, partialtrans4, 
                                DR4_Reco, j0max_new, j1max_new,
                                PU, tin
                            );
                        }
                    }
                }
            } else {
                if constexpr(maxOrder >= 6){
                    do6<T, doPU, doTransfer, maxOrder, true>(
                        dRs, Es, nPart, ans,
                        i0, i1, i2, i3, i4, partial4, 
                        DR4, i0max_new, i1max_new,
                        isPU4,
                        0, 0, 0, 0, 0, 0, 
                        0, 0, 0,
                        PU, tin
                    );
                }
            }
        }
        if(isPU4){
            isPU4 = true;
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
