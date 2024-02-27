#ifndef EECS_FAST6_H
#define EECS_FAST6_H

namespace fastEEC{
    template <typename T, bool doPU, bool doTransfer, unsigned maxOrder, bool nontransfer>
    void do6(const umat& dRs,
             const vector<T>& Es,
             const unsigned nPart,
             result<T>& ans,

             const unsigned i0,
             const unsigned i1,
             const unsigned i2,
             const unsigned i3,
             const unsigned i4,
             const T partial4,
             const unsigned DR4,
             const unsigned i0max,
             const unsigned i1max,
             bool isPU4,

             const unsigned j0,
             const unsigned j1,
             const unsigned j2,
             const unsigned j3,
             const unsigned j4,
             const T partialtrans4,
             const unsigned DR4_Reco,
             const unsigned j0max,
             const unsigned j1max,

             const vector<bool>* const PU = nullptr,
             const transferInputs<T>* const tin = nullptr) {

        T weight6;
        bool isPU5=isPU4;
        unsigned i0max_new=0, i1max_new=0;
        unsigned j0max_new=0, j1max_new=0;
        for(unsigned i5=i4; i5<nPart; ++i5){
            T partial5 = partial4 * Es[i5];

            uvec dRlist = {DR4, dRs[i0][i5],
                                dRs[i1][i5],
                                dRs[i2][i5],
                                dRs[i3][i5],
                                dRs[i4][i5]};
            auto maxel = max_element(dRlist.begin(), dRlist.end());
            unsigned DR5 = *maxel;
            switch(std::distance(dRlist.begin(), maxel)){
                case 0:
                    i0max_new = i0max;
                    i1max_new = i1max;
                    break;
                case 1:
                    i0max_new = 0;
                    i1max_new = 5;
                    break;
                case 2:
                    i0max_new = 1;
                    i1max_new = 5;
                    break;
                case 3:
                    i0max_new = 2;
                    i1max_new = 5;
                    break;
                case 4:
                    i0max_new = 3;
                    i1max_new = 5;
                    break;
                case 5:
                    i0max_new = 4;
                    i1max_new = 5;
                    break;
            }

            T symfac;
            if (i0==i1){
                if(i1==i2){
                    if(i2==i3){
                        if(i3==i4){
                            symfac = (i4==i5) ? 1 : 6; //(6) vs (5, 1)
                        } else {
                            symfac = (i4==i5) ? 15 : 30; //(4, 2) vs (4, 1, 1)
                        }
                    } else {
                        if(i3==i4){
                            symfac = (i4==i5) ? 20 : 60; //(3, 3) vs (3, 2, 1)
                        } else {
                            symfac = (i4==i5) ? 60 : 120; //(3, 1, 2) vs (3, 1, 1, 1)
                        }
                    }
                } else {
                    if(i2==i3){
                        if(i3==i4){
                            symfac = (i4==i5) ? 15 : 60; //(2, 4) vs (2, 3, 1)
                        } else {
                            symfac = (i4==i5) ? 90 : 180; //(2, 2, 2) vs (2, 2, 1, 1)
                        }
                    } else {
                        if(i3==i4){
                            symfac = (i4==i5) ? 60 : 180; //(2, 1, 3) vs (2, 1, 2, 1)
                        } else {
                            symfac = (i4==i5) ? 180 : 360; //(2, 1, 1, 2) vs (2, 1, 1, 1, 1)
                        }
                    }
                }
            } else {
                if(i1==i2){
                    if(i2==i3){
                        if(i3==i4){
                            symfac = (i4==i5) ? 6 : 30; //(1, 5) vs (1, 4, 1)
                        } else {
                            symfac = (i4==i5) ? 60 : 120; //(1, 3, 2) vs (1, 3, 1, 1)
                        }
                    } else {
                        if(i3==i4){
                            symfac = (i4==i5) ? 60 : 180; //(1, 2, 3) vs (1, 2, 2, 1)
                        } else {
                            symfac = (i4==i5) ? 180 : 360; //(1, 2, 1, 2) vs (1, 2, 1, 1, 1)
                        }
                    }
                } else {
                    if(i2==i3){
                        if(i3==i4){
                            symfac = (i4==i5) ? 30 : 120; //(1, 1, 4) vs (1, 1, 3, 1)
                        } else {
                            symfac = (i4==i5) ? 180 : 360; //(1, 1, 2, 2) vs (1, 1, 2, 1, 1)
                        }
                    } else {
                        if(i3==i4){
                            symfac = (i4==i5) ? 120 : 360; //(1, 1, 1, 3) vs (1, 1, 1, 2, 1)
                        } else {
                            symfac = (i4==i5) ? 360 : 720; //(1, 1, 1, 1, 2) vs (1, 1, 1, 1, 1, 1)
                        }
                    }
                }
            }


            weight6 = symfac * partial5;

            if constexpr(nontransfer){
                //accumulate
                ans.wts6[DR5] += weight6;
                if constexpr(doPU){
                    isPU5 = isPU4 || PU->at(i5);
                    if(isPU5){
                        ans.wts6_PU[DR5] += weight6;
                    }
                }
            }

            if constexpr(doTransfer){
                const uvec& adj5 = tin->adj.at(i5);
                if (adj5.size()==0 || partialtrans4==0){
                    continue;
                } else {
                    unsigned j5 = adj5[0];
                    T partialtrans5 = partialtrans4 * tin->ptrans[i5][j5];

                    uvec dRlist_reco = {DR4_Reco, tin->dRs[j0][j5],
                                                  tin->dRs[j1][j5],
                                                  tin->dRs[j2][j5],
                                                  tin->dRs[j3][j5],
                                                  tin->dRs[j4][j5]};
                    auto maxel_reco = max_element(dRlist_reco.begin(), dRlist_reco.end());
                    unsigned DR5_Reco = *maxel_reco;
                    switch(std::distance(dRlist_reco.begin(), maxel_reco)){
                        case 0:
                            j0max_new = j0max;
                            j1max_new = j1max;
                            break;
                        case 1:
                            j0max_new = 0;
                            j1max_new = 5;
                            break;
                        case 2:
                            j0max_new = 1;
                            j1max_new = 5;
                            break;
                        case 3:
                            j0max_new = 2;
                            j1max_new = 5;
                            break;
                        case 4:
                            j0max_new = 3;
                            j1max_new = 5;
                            break;
                        case 5:
                            j0max_new = 4;
                            j1max_new = 5;
                            break;
                    };

                    ans.transfer6[DR5][DR5_Reco] += partialtrans5 * weight6;

                    for(unsigned j=1; j<adj5.size(); ++j){
                        j5 = adj5[j];
                        partialtrans5 = partialtrans4 * tin->ptrans[i5][j5];

                        dRlist_reco[1] = tin->dRs[j0][j5];
                        dRlist_reco[2] = tin->dRs[j1][j5];
                        dRlist_reco[3] = tin->dRs[j2][j5];
                        dRlist_reco[4] = tin->dRs[j3][j5];
                        dRlist_reco[5] = tin->dRs[j4][j5];
                        maxel_reco = max_element(dRlist_reco.begin(), dRlist_reco.end());
                        DR5_Reco = *maxel_reco;
                        switch(std::distance(dRlist_reco.begin(), maxel_reco)){
                            case 0:
                                j0max_new = j0max;
                                j1max_new = j1max;
                                break;
                            case 1:
                                j0max_new = 0;
                                j1max_new = 5;
                                break;
                            case 2:
                                j0max_new = 1;
                                j1max_new = 5;
                                break;
                            case 3:
                                j0max_new = 2;
                                j1max_new = 5;
                                break;
                            case 4:
                                j0max_new = 3;
                                j1max_new = 5;
                                break;
                            case 5:
                                j0max_new = 4;
                                j1max_new = 5;
                                break;
                        };
                        ans.transfer6[DR5][DR5_Reco] += partialtrans5 * weight6;
                    }
                }
            }
        }
        if(isPU5){
            isPU5 = true;
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
