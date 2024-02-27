#ifndef EECS_FAST3_H
#define EECS_FAST3_H

#include "SRothman/EECs/src/fast4.h"

namespace fastEEC{
    template <typename T, bool doPU, bool doTransfer, unsigned maxOrder, bool nontransfer>
    void do3(const umat& dRs,
             const vector<T>& Es,
             const unsigned nPart,
             result<T>& ans,

             const unsigned i0,
             const unsigned i1,
             const T partial1,
             const unsigned DR1,
             bool isPU1,

             const unsigned j0,
             const unsigned j1,
             const T partialtrans1,
             const unsigned DR1_Reco,

             const vector<bool>* const PU = nullptr,
             const transferInputs<T>* const tin = nullptr) {

        T weight3;
        bool isPU2=isPU1;
        unsigned i0max=0, i1max=0;
        unsigned j0max=0, j1max=0;
        for(unsigned i2=i1; i2<nPart; ++i2){
            T partial2 = partial1 * Es[i2];

            uvec dRlist = {DR1, dRs[i0][i2], dRs[i1][i2]};
            auto maxel = max_element(dRlist.begin(), dRlist.end());
            unsigned DR2 = *maxel;
            switch(std::distance(dRlist.begin(), maxel)){
                case 0:
                    i0max = 0;
                    i1max = 1;
                    break;
                case 1:
                    i0max = 0;
                    i1max = 2;
                    break;
                case 2:
                    i0max = 1;
                    i1max = 2;
                    break;
            };

            T symfac;
            if(i0==i1){
                symfac = (i1==i2) ? 1 : 3;
            } else{
                symfac = (i1==i2) ? 3 : 6;
            }
            weight3 = symfac * partial2;

            if constexpr(nontransfer){
                //accumulate
                ans.wts3[DR2] += weight3;
                if constexpr(doPU){
                    isPU2 = isPU1 || PU->at(i2);
                    if(isPU2){
                        ans.wts3_PU[DR2] += weight3;
                    }
                }
            }

            if constexpr(doTransfer){
                const uvec& adj2 = tin->adj.at(i2);
                if(adj2.empty() || partialtrans1==0){
                    if constexpr (maxOrder >=4 && nontransfer){
                        do4<T, doPU, doTransfer, maxOrder, nontransfer>(
                            dRs, Es, nPart, ans,
                            i0, i1, i2, partial2, 
                            DR2, i0max, i1max,
                            isPU2,
                            0, 0, 0, 0, 
                            0, 0, 0,
                            PU, tin
                        );
                    }
                } else {
                    unsigned j2 = adj2[0];
                    T partialtrans2 = partialtrans1 * tin->ptrans[i2][j2];

                    uvec dRlist_reco = {DR1_Reco, tin->dRs[j0][j2], tin->dRs[j1][j2]};
                    auto maxel_reco = max_element(dRlist_reco.begin(), dRlist_reco.end());
                    unsigned DR2_Reco = *maxel_reco;
                    switch(std::distance(dRlist_reco.begin(), maxel_reco)){
                        case 0:
                            j0max = 0;
                            j1max = 1;
                            break;
                        case 1:
                            j0max = 0;
                            j1max = 2;
                            break;
                        case 2:
                            j0max = 1;
                            j1max = 2;
                            break;
                    };

                    ans.transfer3[DR2][DR2_Reco] += partialtrans2 * weight3;


                    if constexpr (maxOrder >=4){
                        do4<T, doPU, doTransfer, maxOrder, nontransfer>(
                            dRs, Es, nPart, ans,
                            i0, i1, i2, partial2, 
                            DR2, i0max, i1max, 
                            isPU2,
                            j0, j1, j2, partialtrans2, 
                            DR2_Reco, j0max, j1max,
                            PU, tin
                        );
                    }

                    for(unsigned j=1; j<adj2.size(); ++j){
                        j2 = adj2[j];
                        partialtrans2 = partialtrans1 * tin->ptrans[i2][j2];

                        dRlist_reco[1] = tin->dRs[j0][j2];
                        dRlist_reco[2] = tin->dRs[j1][j2];
                        maxel_reco = max_element(dRlist_reco.begin(), dRlist_reco.end());
                        DR2_Reco = *maxel_reco;
                        switch(std::distance(dRlist_reco.begin(), maxel_reco)){
                            case 0:
                                j0max = 0;
                                j1max = 1;
                                break;
                            case 1:
                                j0max = 0;
                                j1max = 2;
                                break;
                            case 2:
                                j0max = 1;
                                j1max = 2;
                                break;
                        };

                        ans.transfer3[DR2][DR2_Reco] += partialtrans2 * weight3;

                        if constexpr (maxOrder >= 4){
                            do4<T, doPU, doTransfer, maxOrder, false>(
                                dRs, Es, nPart, ans,
                                i0, i1, i2, partial2, 
                                DR2, i0max, i1max, 
                                isPU2,
                                j0, j1, j2, partialtrans2, 
                                DR2_Reco, j0max, j1max,
                                PU, tin
                            );
                        }
                    }
                }
            } else {
                if constexpr(maxOrder >= 4){
                    do4<T, doPU, doTransfer, maxOrder, true>(
                        dRs, Es, nPart, ans,
                        i0, i1, i2, partial2, 
                        DR2, i0max, i1max,
                        isPU2,
                        0, 0, 0, 0, 
                        0, 0, 0,
                        PU, tin
                    );
                }
            }
        }
        if(isPU2){
            isPU2 = true;
        } else if(i0max){
            i0max = 0;
        } else if(i1max){
            i1max = 0;
        } else if(j0max){
            j0max = 0;
        } else if(j1max){
            j1max = 0;
        }
    }
};

#endif
