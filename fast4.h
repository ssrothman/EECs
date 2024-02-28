#ifndef EECS_FAST4_H
#define EECS_FAST4_H

#include "fast5.h"
#include "resolved4.h"

namespace fastEEC{
    template <typename T, bool doPU, bool doTransfer, unsigned maxOrder, bool nontransfer>
    void do4(const umat& dRs,
             const vector<T>& Es,
             const unsigned nPart,
             const resolvedInputs<T>& rin,
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
            unsigned qi0=0, qi1=0, qi2=0, qi3=0;
            switch(std::distance(dRlist.begin(), maxel)){
                case 0:
                    i0max_new = i0max;
                    i1max_new = i1max;
                    switch(i0max){
                        case 0:
                            switch(i1max){
                                case 1:
                                    qi0 = i0;
                                    qi1 = i1;
                                    qi2 = i2;
                                    qi3 = i3;
                                    break;
                                case 2:
                                    qi0 = i0;
                                    qi1 = i2;
                                    qi2 = i1;
                                    qi3 = i3;
                                    break;
                            };
                            break;
                        case 1:
                            qi0 = i1;
                            qi1 = i2;
                            qi2 = i0;
                            qi3 = i3;
                            break;
                    };
                    break;
                case 1:
                    i0max_new = 0;
                    i1max_new = 3;
                    qi0 = i0;
                    qi1 = i3;
                    qi2 = i1;
                    qi3 = i2;
                    break;
                case 2:
                    i0max_new = 1;
                    i1max_new = 3;
                    qi0 = i1;
                    qi1 = i3;
                    qi2 = i0;
                    qi3 = i2;
                    break;
                case 3:
                    i0max_new = 2;
                    i1max_new = 3;
                    qi0 = i2;
                    qi1 = i3;
                    qi2 = i0;
                    qi3 = i1;
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

            /*
            unsigned fixedshape_idx;
            fixedshape4(qi0, qi1, qi2, qi3, rin, fixedshape_idx);
            */

            if constexpr(nontransfer){
                //accumulate
                ans.wts4[DR3] += weight4;
                //ans.resolved4_fixed[fixedshape_idx][DR3] += weight4;
                if constexpr(doPU){
                    isPU3 = isPU2 || PU->at(i3);
                    if(isPU3){
                        ans.wts4_PU[DR3] += weight4;
                        //ans.resolved4_fixed_PU[fixedshape_idx][DR3] += weight4;
                    }
                }
            }

            unsigned shape_idx, RL_idx, r_idx, ct_idx;
            resolved4(qi0, qi1, qi2, qi3, rin, shape_idx, RL_idx, r_idx, ct_idx); 
            if constexpr(nontransfer){
                ans.resolved4_shapes[shape_idx][RL_idx][r_idx][ct_idx] += weight4;
                if constexpr(doPU){
                    if(isPU3){
                        ans.resolved4_shapes_PU[shape_idx][RL_idx][r_idx][ct_idx] += weight4;
                    }
                }
            }

            if constexpr(doTransfer){
                const uvec& adj3 = tin->adj.at(i3);
                if(adj3.empty() || partialtrans2==0){
                    if constexpr (maxOrder >=5 && nontransfer){
                        do5<T, doPU, doTransfer, maxOrder, nontransfer>(
                            dRs, Es, nPart, rin, ans,
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
                    unsigned qj0=0, qj1=0, qj2=0, qj3=0;
                    switch(std::distance(dRlist_reco.begin(), maxel_reco)){
                        case 0:
                            j0max_new = j0max;
                            j1max_new = j1max;
                            switch(j0max){
                                case 0:
                                    switch(j1max){
                                        case 1:
                                            qj0 = j0;
                                            qj1 = j1;
                                            qj2 = j2;
                                            qj3 = j3;
                                            break;
                                        case 2:
                                            qj0 = j0;
                                            qj1 = j2;
                                            qj2 = j1;
                                            qj3 = j3;
                                            break;
                                    };
                                    break;
                                case 1:
                                    qj0 = j1;
                                    qj1 = j2;
                                    qj2 = j0;
                                    qj3 = j3;
                                    break;
                            };
                            break;
                        case 1:
                            j0max_new = 0;
                            j1max_new = 3;
                            qj0 = j0;
                            qj1 = j3;
                            qj2 = j1;
                            qj3 = j2;
                            break;
                        case 2:
                            j0max_new = 1;
                            j1max_new = 3;
                            qj0 = j1;
                            qj1 = j3;
                            qj2 = j0;
                            qj3 = j2;
                            break;
                        case 3:
                            j0max_new = 2;
                            j1max_new = 3;
                            qj0 = j2;
                            qj1 = j3;
                            qj2 = j0;
                            qj3 = j1;
                            break;
                    };
                    ans.transfer4[DR3][DR3_Reco] += partialtrans3 * weight4;

                    /*
                    unsigned shapeidx_fixed_reco;
                    fixedshape4(qj0, qj1, qj2, qj3, rin, shapeidx_fixed_reco);
                    ans.transfer_res4_fixed[fixedshape_idx][DR3][shapeidx_fixed_reco][DR3_Reco] += partialtrans3 * weight4;
                    */

                    unsigned shape_idx_reco, RL_idx_reco, r_idx_reco, ct_idx_reco;
                    resolved4(qj0, qj1, qj2, qj3, tin->rin, shape_idx_reco, RL_idx_reco, r_idx_reco, ct_idx_reco);
                    ans.transfer_res4_shapes[shape_idx][RL_idx][r_idx][ct_idx][shape_idx_reco][RL_idx_reco][r_idx_reco][ct_idx_reco] += partialtrans3 * weight4;

                    if constexpr (maxOrder >=5){
                        do5<T, doPU, doTransfer, maxOrder, nontransfer>(
                            dRs, Es, nPart, rin, ans,
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
                                switch(j0max){
                                    case 0:
                                        switch(j1max){
                                            case 1:
                                                qj0 = j0;
                                                qj1 = j1;
                                                qj2 = j2;
                                                qj3 = j3;
                                                break;
                                            case 2:
                                                qj0 = j0;
                                                qj1 = j2;
                                                qj2 = j1;
                                                qj3 = j3;
                                                break;
                                        };
                                        break;
                                    case 1:
                                        qj0 = j1;
                                        qj1 = j2;
                                        qj2 = j0;
                                        qj3 = j3;
                                        break;
                                };
                                break;
                            case 1:
                                j0max_new = 0;
                                j1max_new = 3;
                                qj0 = j0;
                                qj1 = j3;
                                qj2 = j1;
                                qj3 = j2;
                                break;
                            case 2:
                                j0max_new = 1;
                                j1max_new = 3;
                                qj0 = j1;
                                qj1 = j3;
                                qj2 = j0;
                                qj3 = j2;
                                break;
                            case 3:
                                j0max_new = 2;
                                j1max_new = 3;
                                qj0 = j2;
                                qj1 = j3;
                                qj2 = j0;
                                qj3 = j1;
                                break;
                        };

                        ans.transfer4[DR3][DR3_Reco] += partialtrans3 * weight4;

                        /*
                        fixedshape4(qj0, qj1, qj2, qj3, rin, shapeidx_fixed_reco);
                        ans.transfer_res4_fixed[fixedshape_idx][DR3][shapeidx_fixed_reco][DR3_Reco] += partialtrans3 * weight4;
                        */

                        resolved4(qj0, qj1, qj2, qj3, tin->rin, shape_idx_reco, RL_idx_reco, r_idx_reco, ct_idx_reco);
                        ans.transfer_res4_shapes[shape_idx][RL_idx][r_idx][ct_idx][shape_idx_reco][RL_idx_reco][r_idx_reco][ct_idx_reco] += partialtrans3 * weight4;

                        if constexpr (maxOrder >= 5){
                            do5<T, doPU, doTransfer, maxOrder, false>(
                                dRs, Es, nPart, rin, ans,
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
                        dRs, Es, nPart, rin, ans,
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
