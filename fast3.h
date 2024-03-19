#ifndef EECS_FAST3_H
#define EECS_FAST3_H

#include "SRothman/EECs/src/fast4.h"

namespace fastEEC{
    template <typename T, bool doPU, bool doTransfer, unsigned maxOrder, bool nontransfer, bool doRes3, bool doRes4, bool doRes4Fixed>
    void do3(const umat& dRs,
             const vector<T>& Es,
             const unsigned nPart,
             const resolvedInputs<T>& rin,
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
        bool isPU2 __attribute__((unused));
        unsigned i0max __attribute__((unused));
        unsigned i1max __attribute__((unused));
        unsigned j0max __attribute__((unused));
        unsigned j1max __attribute__((unused));
        for(unsigned i2=i1; i2<nPart; ++i2){
            T partial2 = partial1 * Es[i2];

            uvec dRlist = {DR1, dRs[i0][i2], 
                                dRs[i1][i2]};
            auto maxel = max_element(dRlist.begin(), dRlist.end());
            unsigned DR2 = *maxel;
            unsigned qi0 __attribute__((unused));
            unsigned qi1 __attribute__((unused));
            unsigned qi2 __attribute__((unused)); //used for resolved
            switch(std::distance(dRlist.begin(), maxel)){
                case 0:
                    i0max = 0;
                    i1max = 1;
                    if constexpr(doRes3){
                        qi0 = i0;
                        qi1 = i1;
                        qi2 = i2;
                    }
                    break;
                case 1:
                    i0max = 0;
                    i1max = 2;
                    if constexpr (doRes3){
                        qi0 = i0;
                        qi1 = i2;
                        qi2 = i1;
                    }
                    break;
                case 2:
                    i0max = 1;
                    i1max = 2;
                    if constexpr (doRes3){
                        qi0 = i1;
                        qi1 = i2;
                        qi2 = i0;
                    }
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
                (*ans.wts3)[DR2] += weight3;
                if constexpr(doPU){
                    isPU2 = isPU1 || PU->at(i2);
                    if(isPU2){
                        (*ans.wts3_PU)[DR2] += weight3;
                    }
                }
            }

            unsigned RL_idx __attribute__((unused));
            unsigned xi_idx __attribute__((unused));
            unsigned phi_idx __attribute__((unused));
            if constexpr (doRes3){
                T RL = rin.floatDRs[qi0][qi1];
                T RM = rin.floatDRs[qi0][qi2];
                T RS = rin.floatDRs[qi1][qi2];
                if (RS > RM){
                    std::swap(RS, RM);
                }
                
                T xi;
                if (RM==0){
                    xi = 0;
                } else {
                    xi = RS/RM;
                }
                T phi;
                if (RS==0){
                    phi = 0;
                } else {
                    phi = std::abs(std::asin(std::sqrt(1 - square((RL-RM)/RS))));
                }

                RL_idx = static_cast<unsigned>(rin.coarseRL->index(RL) + 1);
                xi_idx = static_cast<unsigned>(rin.xi->index(xi) + 1);
                phi_idx = static_cast<unsigned>(rin.phi->index(phi) + 1);
                if constexpr(nontransfer){
                    //accumulate
                    (*ans.resolved3)[RL_idx][xi_idx][phi_idx] += weight3;
                    if constexpr(doPU){
                        if(isPU2){
                            (*ans.resolved3_PU)[RL_idx][xi_idx][phi_idx] += weight3;
                        }
                    }
                }
            }

            if constexpr(doTransfer){
                const uvec& adj2 = tin->adj.at(i2);
                if(adj2.empty() || partialtrans1==0){
                    if constexpr (maxOrder >=4 && nontransfer){
                        do4<T, doPU, doTransfer, maxOrder, nontransfer, doRes4, doRes4Fixed>(
                            dRs, Es, nPart, rin, ans,
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
                    unsigned qj0 __attribute__((unused));
                    unsigned qj1 __attribute__((unused));
                    unsigned qj2 __attribute__((unused)); //used for resolved
                    switch(std::distance(dRlist_reco.begin(), maxel_reco)){
                        case 0:
                            j0max = 0;
                            j1max = 1;
                            if constexpr (doRes3){
                                qj0 = j0;
                                qj1 = j1;
                                qj2 = j2;
                            }
                            break;
                        case 1:
                            j0max = 0;
                            j1max = 2;
                            if constexpr (doRes3){
                                qj0 = j0;
                                qj1 = j2;
                                qj2 = j1;
                            }
                            break;
                        case 2:
                            j0max = 1;
                            j1max = 2;
                            if constexpr (doRes3){
                                qj0 = j1;
                                qj1 = j2;
                                qj2 = j0;
                            }
                            break;
                    };

                    (*ans.transfer3)[DR2][DR2_Reco] += partialtrans2 * weight3;

                    if constexpr (doRes3){
                        T RL_reco = tin->rin.floatDRs[qj0][qj1]; 
                        T RM_reco = tin->rin.floatDRs[qj0][qj2];
                        T RS_reco = tin->rin.floatDRs[qj1][qj2];
                        if (RS_reco > RM_reco){
                            std::swap(RS_reco, RM_reco);
                        }

                        T xi_reco;
                        if (RM_reco==0){
                            xi_reco = 0;
                        } else {
                            xi_reco = RS_reco/RM_reco;
                        }
                        T phi_reco;
                        if (RS_reco==0){
                            phi_reco = 0;
                        } else {
                            phi_reco = std::abs(std::asin(std::sqrt(1 - square((RL_reco-RM_reco)/RS_reco))));
                        }
                        unsigned RL_reco_idx = static_cast<unsigned>(rin.coarseRL->index(RL_reco) + 1);
                        unsigned xi_reco_idx = static_cast<unsigned>(rin.xi->index(xi_reco) + 1);
                        unsigned phi_reco_idx = static_cast<unsigned>(rin.phi->index(phi_reco) + 1);

                        (*ans.transfer_res3)[RL_idx][xi_idx][phi_idx][RL_reco_idx][xi_reco_idx][phi_reco_idx] 
                            += partialtrans2 * weight3;
                    }

                    if constexpr (maxOrder >=4){
                        do4<T, doPU, doTransfer, maxOrder, nontransfer, doRes4, doRes4Fixed>(
                            dRs, Es, nPart, rin, ans,
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
                                if constexpr (doRes3){
                                    qj0 = j0;
                                    qj1 = j1;
                                    qj2 = j2;
                                }
                                break;
                            case 1:
                                j0max = 0;
                                j1max = 2;
                                if constexpr (doRes3){
                                    qj0 = j0;
                                    qj1 = j2;
                                    qj2 = j1;
                                }
                                break;
                            case 2:
                                j0max = 1;
                                j1max = 2;
                                if constexpr (doRes3){
                                    qj0 = j1;
                                    qj1 = j2;
                                    qj2 = j0;
                                }
                                break;
                        };

                        (*ans.transfer3)[DR2][DR2_Reco] += partialtrans2 * weight3;

                        if constexpr (doRes3){
                            T RL_reco = tin->rin.floatDRs[qj0][qj1];
                            T RM_reco = tin->rin.floatDRs[qj0][qj2];
                            T RS_reco = tin->rin.floatDRs[qj1][qj2];
                            if (RS_reco > RM_reco){
                                std::swap(RS_reco, RM_reco);
                            }

                            T xi_reco;
                            if (RM_reco==0){
                                xi_reco = 0;
                            } else {
                                xi_reco = RS_reco/RM_reco;
                            }

                            T phi_reco;
                            if (RS_reco==0){
                                phi_reco = 0;
                            } else {
                                phi_reco = std::abs(std::asin(std::sqrt(1 - square((RL_reco-RM_reco)/RS_reco))));
                            }

                            unsigned RL_reco_idx = static_cast<unsigned>(rin.coarseRL->index(RL_reco) + 1);
                            unsigned xi_reco_idx = static_cast<unsigned>(rin.xi->index(xi_reco) + 1);
                            unsigned phi_reco_idx = static_cast<unsigned>(rin.phi->index(phi_reco) + 1);

                            (*ans.transfer_res3)[RL_idx][xi_idx][phi_idx][RL_reco_idx][xi_reco_idx][phi_reco_idx] 
                                += partialtrans2 * weight3;
                        }

                        if constexpr (maxOrder >= 4){
                            do4<T, doPU, doTransfer, maxOrder, false, doRes4, doRes4Fixed>(
                                dRs, Es, nPart, rin, ans,
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
                    do4<T, doPU, doTransfer, maxOrder, true, doRes4, doRes4Fixed>(
                        dRs, Es, nPart, rin, ans,
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
    }
};

#endif
