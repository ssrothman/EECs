#ifndef EECS_FAST_H

#include <boost/histogram.hpp>
#include "boost/multi_array.hpp"

#include <cassert>

#include "SRothman/SimonTools/src/jets.h"

#include "fastStructs.h"
#include "fastTransfer6_1.h"
#include "fastTransfer6_2.h"
#include "fastTransfer6_3.h"
#include "fastTransfer6_4.h"
#include "fastTransfer6_5.h"
#include "fastTransfer6_6.h"

#include "faststart.h"

namespace fastEEC{
    double getNormFact(const jet& J, const normType& nt){
        switch (nt){
            case RAWPT:
                return J.rawpt;
            case SUMPT:
                return J.sumpt;
            case CORRPT:
                return J.pt;
            default:
                throw std::invalid_argument("Invalid normType");
        }
    }

    template <typename T>
    void getEs(vector<T>& ans, const jet& J, const normType nt){
        ans.clear();
        ans.reserve(J.nPart);

        double normFact = getNormFact(J, nt);
        for(unsigned i=0; i<J.nPart; ++i){
            ans.push_back(J.particles.at(i).pt/normFact);
        }
    }

    template <typename T>
    void getEtasPhis(vector<T>& etas, vector<T>& phis, const jet& J){
        etas.clear();
        phis.clear();

        etas.reserve(J.nPart);
        phis.reserve(J.nPart);

        for(unsigned i=0; i<J.nPart; ++i){
            etas.push_back(J.particles.at(i).eta);
            phis.push_back(J.particles.at(i).phi);
        }
    }

    void getDRs(umat& ans, const jet& J, const axisptr& ax){
        ans.resize(extents[J.nPart][J.nPart]);
        for (unsigned i0=0; i0<J.nPart; ++i0){
            for (unsigned i1=0; i1<J.nPart; ++i1){
                double deltaR = dR(J.particles[i0], J.particles[i1]);
                unsigned idx = static_cast<unsigned>(ax->index(deltaR) + 1);
                ans[i0][i1] = idx;
            }
        }
    }

    template <typename T>
    void getFloatDRs(multi_array<T, 2>& ans, const jet& J){
        ans.resize(extents[J.nPart][J.nPart]);
        for (unsigned i0=0; i0<J.nPart; ++i0){
            for (unsigned i1=0; i1<J.nPart; ++i1){
                T deltaR = dR(J.particles[i0], J.particles[i1]);
                ans[i0][i1] = deltaR;
            }
        }
    }

    template <typename T>
    void getPtrans(multi_array<T, 2>& ans, const arma::mat& ptrans){//NB we transpose for faster iteration
        ans.resize(extents[ptrans.n_cols][ptrans.n_rows]);
        for(unsigned i=0; i<ptrans.n_rows; ++i){
            for(unsigned j=0; j<ptrans.n_cols; ++j){
                ans[j][i] = ptrans(i,j);
            }
        }
    }

    template <typename T, bool doPU>
    void fastSecondOrder(const umat& dRs, 
                         const vector<T>& Es,
                         const unsigned NDR,
                         result<T>& ans,
                         const vector<bool>* const PU = nullptr){
        ans.wts2.resize(NDR);
        fill(ans.wts2.begin(), ans.wts2.end(), 0);
        if constexpr (doPU){
            ans.wts2_PU.resize(NDR);
            fill(ans.wts2_PU.begin(), ans.wts2_PU.end(), 0);
        }

        const unsigned nPart = Es.size();

        T weight = 1;
        bool isPU=false;

        for(unsigned i0=0; i0<nPart; ++i0){
            T partial0 = Es[i0];
            unsigned DR0 = 0;

            //partition (2)
            weight = square(partial0); 
            ans.wts2[DR0] += weight;

            if constexpr (doPU){
                isPU = PU->at(i0);
                if(isPU){
                    ans.wts2_PU[DR0] += weight;
                }
            }

            for(unsigned i1=i0+1; i1<nPart; ++i1){
                //partition (1,1)
                weight = partial0 * Es[i1] * 2; 
                unsigned DR1 = dRs[i0][i1];
                ans.wts2[DR1] += weight;

                if constexpr (doPU){
                    isPU = isPU || PU->at(i1);
                    if(isPU){
                        ans.wts2_PU[DR1] += weight;
                    }
                }
            }
        }
    }

    template <typename T, bool doPU>
    void fastThirdOrder(const umat& dRs,
                        const vector<T>& Es,
                        const unsigned NDR,
                        result<T>& ans,
                         const vector<bool>* const PU = nullptr){
        ans.wts2.resize(NDR);
        ans.wts3.resize(NDR);
        fill(ans.wts2.begin(), ans.wts2.end(), 0);
        fill(ans.wts3.begin(), ans.wts3.end(), 0);

        if constexpr (doPU){
            ans.wts2_PU.resize(NDR);
            ans.wts3_PU.resize(NDR);
            fill(ans.wts2_PU.begin(), ans.wts2_PU.end(), 0);
            fill(ans.wts3_PU.begin(), ans.wts3_PU.end(), 0);
        }

        const unsigned nPart = Es.size();

        T weight2, weight3;
        bool isPU=false;
        for(unsigned i0=0; i0<nPart; ++i0){
            T E0 = Es[i0];
            T partial0 = E0;

            unsigned DR0 = 0;

            //partition (2)
            weight2 = square(partial0);
            ans.wts2[DR0] += weight2;
            //partition (3)
            weight3 = weight2 * partial0;
            ans.wts3[DR0] += weight3;

            if constexpr (doPU){
                isPU = PU->at(i0);
                if(isPU){
                    ans.wts2_PU[DR0] += weight2;
                    ans.wts3_PU[DR0] += weight3;
                }
            }

            for(unsigned i1=i0+1; i1<nPart; ++i1){
                T E1 = Es[i1];
                T partial1 = partial0 * E1;

                unsigned DR1 = dRs[i0][i1];

                //partition (1,1)
                weight2 = partial1 * 2;
                ans.wts2[DR1] += weight2;
                //partition (2,1)
                weight3 = 3 * partial1 * (E0 + E1);//case i1==i2
                ans.wts3[DR1] += weight3;

                if constexpr (doPU){
                    isPU = isPU || PU->at(i1);
                    if(isPU){
                        ans.wts2_PU[DR1] += weight2;
                        ans.wts3_PU[DR1] += weight3;
                    }
                }

                for(unsigned i2=i1+1; i2<nPart; ++i2){
                    //partition (1,1,1)
                    weight3 = partial1 * Es[i2] * 6;
                    unsigned DR2 = max({DR1, dRs[i0][i2], 
                                             dRs[i1][i2]});
                    ans.wts3[DR2] += weight3;

                    if constexpr (doPU){
                        isPU = isPU || PU->at(i2);
                        if(isPU){
                            ans.wts3_PU[DR2] += weight3;
                        }
                    }
                }
            }
        }
    }

    template <typename T, bool doPU>
    void fastFourthOrder(const umat& dRs,
                         const vector<T>& Es,
                         const unsigned NDR,
                         result<T>& ans,
                         const vector<bool>* const PU = nullptr){
        ans.wts2.resize(NDR);
        ans.wts3.resize(NDR);
        ans.wts4.resize(NDR);
        fill(ans.wts2.begin(), ans.wts2.end(), 0);
        fill(ans.wts3.begin(), ans.wts3.end(), 0);
        fill(ans.wts4.begin(), ans.wts4.end(), 0);

        if constexpr (doPU){
            ans.wts2_PU.resize(NDR);
            ans.wts3_PU.resize(NDR);
            ans.wts4_PU.resize(NDR);
            fill(ans.wts2_PU.begin(), ans.wts2_PU.end(), 0);
            fill(ans.wts3_PU.begin(), ans.wts3_PU.end(), 0);
            fill(ans.wts4_PU.begin(), ans.wts4_PU.end(), 0);
        }

        const unsigned nPart = Es.size();

        T weight2, weight3, weight4;
        bool isPU=false;

        for(unsigned i0=0; i0<nPart; ++i0){
            T E0 = Es[i0];
            T sqE0 = square(E0);
            T partial0 = E0;

            unsigned DR0 = 0;

            //partition (2)
            weight2 = square(partial0);
            ans.wts2[DR0] += weight2;
            //partition (3)
            weight3 = weight2 * partial0;
            ans.wts3[DR0] += weight3;
            //partition (4)
            weight4 = weight3 * partial0;
            ans.wts4[DR0] += weight4;

            if constexpr (doPU){
                isPU = PU->at(i0);
                if(isPU){
                    ans.wts2_PU[DR0] += weight2;
                    ans.wts3_PU[DR0] += weight3;
                    ans.wts4_PU[DR0] += weight4;
                }
            }

            for(unsigned i1=i0+1; i1<nPart; ++i1){
                T E1 = Es[i1];
                T partial1 = partial0 * E1;
                
                unsigned DR1 = dRs[i0][i1];

                //partition (1,1)
                weight2 = 2 * partial1;
                ans.wts2[DR1] += weight2;
                //partition (2,1)
                weight3 = 3 * partial1 * (E0 + E1);
                ans.wts3[DR1] += weight3;
                //partition (2,2)
                weight4 = 6 * square(partial1);
                //partition (3,1)
                weight4 += 4 * partial1 * (sqE0 + square(E1));
                ans.wts4[DR1] += weight4;

                if constexpr (doPU){
                    isPU = isPU || PU->at(i1);
                    if(isPU){
                        ans.wts2_PU[DR1] += weight2;
                        ans.wts3_PU[DR1] += weight3;
                        ans.wts4_PU[DR1] += weight4;
                    }
                }

                for(unsigned i2=i1+1; i2<nPart; ++i2){
                    T E2 = Es[i2];
                    T partial2 = partial1 * E2;

                    unsigned DR2 = max({DR1, dRs[i0][i2], 
                                             dRs[i1][i2]});

                    //partition (1,1,1)
                    weight3 = partial2 * 6;
                    ans.wts3[DR2] += weight3;
                    //partition (2,1,1)
                    weight4 = 12 * partial2 * (E0 + E1 + E2);
                    ans.wts4[DR2] += weight4;

                    if constexpr (doPU){
                        isPU = isPU || PU->at(i2);
                        if(isPU){
                            ans.wts3_PU[DR2] += weight3;
                            ans.wts4_PU[DR2] += weight4;
                        }
                    }

                    for(unsigned i3=i2+1; i3<nPart; ++i3){
                        unsigned DR3 = max({DR2, dRs[i0][i3],
                                                 dRs[i1][i3], 
                                                 dRs[i2][i3]});
                        //partition (1,1,1,1)
                        weight4 = partial2 * Es[i3] * 24;
                        ans.wts4[DR3] += weight4;

                        if constexpr (doPU){
                            isPU = isPU || PU->at(i3);
                            if(isPU){
                                ans.wts4_PU[DR3] += weight4;
                            }
                        }
                    }
                }
            }
        }
    }

    template <typename T, bool doPU>
    void fastFifthOrder(const umat& dRs,
                        const vector<T>& Es,
                        const unsigned NDR,
                        result<T>& ans,
                         const vector<bool>* const PU = nullptr){
        ans.wts2.resize(NDR);
        ans.wts3.resize(NDR);
        ans.wts4.resize(NDR);
        ans.wts5.resize(NDR);

        fill(ans.wts2.begin(), ans.wts2.end(), 0);
        fill(ans.wts3.begin(), ans.wts3.end(), 0);
        fill(ans.wts4.begin(), ans.wts4.end(), 0);
        fill(ans.wts5.begin(), ans.wts5.end(), 0);

        if constexpr (doPU){
            ans.wts2_PU.resize(NDR);
            ans.wts3_PU.resize(NDR);
            ans.wts4_PU.resize(NDR);
            ans.wts5_PU.resize(NDR);
            fill(ans.wts2_PU.begin(), ans.wts2_PU.end(), 0);
            fill(ans.wts3_PU.begin(), ans.wts3_PU.end(), 0);
            fill(ans.wts4_PU.begin(), ans.wts4_PU.end(), 0);
            fill(ans.wts5_PU.begin(), ans.wts5_PU.end(), 0);
        }

        const unsigned nPart = Es.size();

        T weight2, weight3, weight4, weight5;
        bool isPU=false;

        for(unsigned i0=0; i0<nPart; ++i0){
            T E0 = Es[i0];
            T sqE0 = square(E0);

            T partial0 = E0;

            unsigned DR0 = 0;
            //partition (2)
            weight2 = square(partial0);
            ans.wts2[DR0] += weight2;
            //partition (3)
            weight3 = weight2 * partial0;
            ans.wts3[DR0] += weight3;
            //partition (4)
            weight4 = weight3 * partial0;
            ans.wts4[DR0] += weight4;
            //partition (5)
            weight5 = weight4 * partial0;
            ans.wts5[DR0] += weight5;

            if constexpr (doPU){
                isPU = isPU || PU->at(i0);
                if(isPU){
                    ans.wts2_PU[DR0] += weight2;
                    ans.wts3_PU[DR0] += weight3;
                    ans.wts4_PU[DR0] += weight4;
                    ans.wts5_PU[DR0] += weight5;
                }
            }

            for(unsigned i1=i0+1; i1<nPart; ++i1){
                T E1 = Es[i1];
                T sqE1 = square(E1);

                T partial1 = partial0 * E1;

                //case two distinct values
                unsigned DR1 = dRs[i0][i1];
                
                //partition (1, 1)
                weight2 = 2 * partial1;
                ans.wts2[DR1] += weight2;

                //partition (2, 1)
                weight3 = 3 * partial1 * (E0 + E1);
                ans.wts3[DR1] += weight3;

                //partition (2, 2)
                weight4 = 6 * square(partial1);
                //partition (3, 1)
                weight4 += 4 * partial1 * (sqE0 + sqE1);
                ans.wts4[DR1] += weight4;

                //partition (4, 1)
                weight5 = 5 * partial1 * (sqE1 * E1 +
                                          sqE0 * E0);
                //partition (3, 2)
                weight5 += 10 * square(partial1) * (E0 + E1);
                ans.wts5[DR1] += weight5;

                if constexpr (doPU){
                    isPU = isPU || PU->at(i1);
                    if(isPU){
                        ans.wts2_PU[DR1] += weight2;
                        ans.wts3_PU[DR1] += weight3;
                        ans.wts4_PU[DR1] += weight4;
                        ans.wts5_PU[DR1] += weight5;
                    }
                }

                for(unsigned i2=i1+1; i2<nPart; ++i2){
                    T E2 = Es[i2];
                    T sqE2 = square(E2);

                    T partial2 = partial1 * E2;

                    unsigned DR2 = max({DR1, dRs[i0][i2], 
                                             dRs[i1][i2]});

                    //partition (1, 1, 1)
                    weight3 = 6 * partial2;
                    ans.wts3[DR2] += weight3;
                    //partition (2, 1, 1)
                    weight4 = 12 * partial2 * (E0 + E1 + E2);
                    ans.wts4[DR2] += weight4;
                    //partition (3, 1, 1)
                    weight5 = 20 * partial2 * (sqE0 + sqE1 + sqE2);
                    //partition (2, 2, 1)
                    weight5 += 30 * partial2 * (E0 * E1 +
                                                E0 * E2 +
                                                E1 * E2);
                    ans.wts5[DR2] += weight5;

                    if constexpr (doPU){
                        isPU = isPU || PU->at(i2);
                        if(isPU){
                            ans.wts3_PU[DR2] += weight3;
                            ans.wts4_PU[DR2] += weight4;
                            ans.wts5_PU[DR2] += weight5;
                        }
                    }

                    for(unsigned i3=i2+1; i3<nPart; ++i3){
                        T E3 = Es[i3];

                        T partial3 = partial2 * E3;

                        //case four distinct values
                        unsigned DR3 = max({DR2, dRs[i0][i3], 
                                                 dRs[i1][i3], 
                                                 dRs[i2][i3]});
                        //partition (1, 1, 1, 1)
                        weight4 = 24 * partial3;
                        ans.wts4[DR3] += weight4;
                        //partition (2, 1, 1, 1)
                        weight5 = 60 * partial3 * (E0 + E1 + E2 + E3);
                        ans.wts5[DR3] += weight5;

                        if constexpr (doPU){
                            isPU = isPU || PU->at(i3);
                            if(isPU){
                                ans.wts4_PU[DR3] += weight4;
                                ans.wts5_PU[DR3] += weight5;
                            }
                        }

                        for(unsigned i4=i3+1; i4<nPart; ++i4){
                            //case five distinct values
                            weight5 = partial3 * Es[i4] * 120;
                            unsigned DR4 = max({DR3, dRs[i0][i4], 
                                                     dRs[i1][i4], 
                                                     dRs[i2][i4],
                                                     dRs[i3][i4]});
                            //partition (1, 1, 1, 1, 1)
                            ans.wts5[DR4] += weight5;

                            if constexpr (doPU){
                                isPU = isPU || PU->at(i4);
                                if(isPU){
                                    ans.wts5_PU[DR4] += weight5;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    template <typename T, bool doPU, bool doTransfer>
    void fastSixthOrder2(const umat& dRs,
                         const vector<T>& Es,
                         const unsigned NDR,
                         result<T>& ans,
                         const vector<bool>* const PU = nullptr,
                         const transferInputs<T>* const tin = nullptr){
        //clear<T, doPU, doTransfer, 6>(NDR, ans);

        unsigned nPart = Es.size();

        T weight2, weight3, weight4, weight5, weight6;
        bool isPU;

        for(unsigned i0=0; i0<nPart; ++i0){
            T E0 = Es[i0];
            if constexpr (doPU){
                isPU = PU->at(i0);
            }

            for(unsigned i1=i0; i1<nPart; ++i1){
                T E1 = Es[i1];

                T partial2 = E0 * E1;
                unsigned DR2 = dRs[i0][i1];

                T symfac = (i0==i1) ? 1 : 2;
                weight2 = symfac * partial2;
                ans.wts2[DR2] += weight2;

                if constexpr(doPU){
                    isPU = isPU || PU->at(i1);
                    if(isPU){
                        ans.wts2_PU[DR2] += weight2;
                    }
                }

                for(unsigned i2=i1; i2<nPart; ++i2){
                    T E2 = Es[i2];

                    T partial3 = partial2 * E2;
                    unsigned DR3 = max({DR2, dRs[i0][i2],
                                             dRs[i1][i2]});

                    T symfac;
                    if(i0==i1){
                        symfac = (i1==i2) ? 1 : 3;
                    } else{
                        symfac = (i1==i2) ? 3 : 6;
                    }

                    weight3 = symfac * partial3;
                    ans.wts3[DR3] += weight3;

                    if constexpr (doPU){
                        isPU = isPU || PU->at(i2);
                        if(isPU){
                            ans.wts3_PU[DR3] += weight3;
                        }
                    }

                    for(unsigned i3=i2; i3<nPart; ++i3){
                        T E3 = Es[i3];

                        T partial4 = partial3 * E3;
                        unsigned DR4 = max({DR3, dRs[i0][i3],
                                                 dRs[i1][i3],
                                                 dRs[i2][i3]});

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
                        weight4 = symfac * partial4;
                        ans.wts4[DR4] += weight4;

                        if constexpr (doPU){
                            isPU = isPU || PU->at(i3);
                            if(isPU){
                                ans.wts4_PU[DR4] += weight4;
                            }
                        }

                        for(unsigned i4=i3; i4<nPart; ++i4){
                            T E4 = Es[i4];

                            T partial5 = partial4 * E4;
                            unsigned DR5 = max({DR4, dRs[i0][i4],
                                                     dRs[i1][i4],
                                                     dRs[i2][i4],
                                                     dRs[i3][i4]});

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
                            weight5 = symfac * partial5;
                            ans.wts5[DR5] += weight5;

                            if constexpr (doPU){
                                isPU = isPU || PU->at(i4);
                                if(isPU){
                                    ans.wts5_PU[DR5] += weight5;
                                }
                            }

                            for(unsigned i5=i4; i5<nPart; ++i5){
                                T E5 = Es[i5];

                                T partial6 = partial5 * E5;
                                unsigned DR6 = max({DR5, dRs[i0][i5],
                                                         dRs[i1][i5],
                                                         dRs[i2][i5],
                                                         dRs[i3][i5],
                                                         dRs[i4][i5]});

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

                                weight6 = symfac * partial6;
                                ans.wts6[DR6] += weight6;

                                if constexpr (doPU){
                                    isPU = isPU || PU->at(i5);
                                    if(isPU){
                                        ans.wts6_PU[DR6] += weight6;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

    }


    template <typename T, bool doPU, bool doTransfer>
    void fastSixthOrder(const umat& dRs,
                        const vector<T>& Es,
                        const unsigned NDR,
                        result<T>& ans,
                        const vector<bool>* const PU = nullptr,
                        const transferInputs<T>* const tin = nullptr){
        ans.wts2.resize(NDR);
        ans.wts3.resize(NDR);
        ans.wts4.resize(NDR);
        ans.wts5.resize(NDR);
        ans.wts6.resize(NDR);

        fill(ans.wts2.begin(), ans.wts2.end(), 0);
        fill(ans.wts3.begin(), ans.wts3.end(), 0);
        fill(ans.wts4.begin(), ans.wts4.end(), 0);
        fill(ans.wts5.begin(), ans.wts5.end(), 0);
        fill(ans.wts6.begin(), ans.wts6.end(), 0);

        if constexpr (doPU){
            ans.wts2_PU.resize(NDR);
            ans.wts3_PU.resize(NDR);
            ans.wts4_PU.resize(NDR);
            ans.wts5_PU.resize(NDR);
            ans.wts6_PU.resize(NDR);

            fill(ans.wts2_PU.begin(), ans.wts2_PU.end(), 0);
            fill(ans.wts3_PU.begin(), ans.wts3_PU.end(), 0);
            fill(ans.wts4_PU.begin(), ans.wts4_PU.end(), 0);
            fill(ans.wts5_PU.begin(), ans.wts5_PU.end(), 0);
            fill(ans.wts6_PU.begin(), ans.wts6_PU.end(), 0);
        }

        if constexpr (doTransfer){
            ans.transfer2.resize(extents[NDR][NDR]);
            ans.transfer3.resize(extents[NDR][NDR]);
            ans.transfer4.resize(extents[NDR][NDR]);
            ans.transfer5.resize(extents[NDR][NDR]);
            ans.transfer6.resize(extents[NDR][NDR]);

            for(unsigned i=0; i<NDR; ++i){
                for(unsigned j=0; j<NDR; ++j){
                    ans.transfer2[i][j] = 0;
                    ans.transfer3[i][j] = 0;
                    ans.transfer4[i][j] = 0;
                    ans.transfer5[i][j] = 0;
                    ans.transfer6[i][j] = 0;
                }
            }
        }

        const unsigned nPart = Es.size();

        T weight2, weight3, weight4, weight5, weight6;
        bool isPU=false;

        for(unsigned i0=0; i0<nPart; ++i0){
            T E0 = Es[i0]; 
            T sqE0 = square(E0);

            T partial0 = E0;

            unsigned DR0 = 0;

            //partition (2)
            weight2 = square(partial0);
            ans.wts2[DR0] += weight2;
            //partition (3)
            weight3 = weight2 * partial0;
            ans.wts3[DR0] += weight3;
            //partition (4)
            weight4 = weight3 * partial0;
            ans.wts4[DR0] += weight4;
            //partition (5)
            weight5 = weight4 * partial0;
            ans.wts5[DR0] += weight5;
            //partition (6)
            weight6 = weight5 * partial0;
            ans.wts6[DR0] += weight6;

            if constexpr (doPU){
                isPU = PU->at(i0);
                if(isPU){
                    ans.wts2_PU[DR0] += weight2;
                    ans.wts3_PU[DR0] += weight3;
                    ans.wts4_PU[DR0] += weight4;
                    ans.wts5_PU[DR0] += weight5;
                    ans.wts6_PU[DR0] += weight6;
                }
            }//end if doPU
            
            if constexpr (doTransfer){
                transfer6_1(i0, DR0, 
                            weight2, 
                            weight3, 
                            weight4,
                            weight5,
                            weight6,
                            tin, ans);
            }

            for(unsigned i1=i0+1; i1<nPart; ++i1){
                T E1 = Es[i1];
                T sqE1 = square(E1);

                T partial1 = partial0 * E1;
                T sqpartial1 = square(partial1);

                unsigned DR1 = dRs[i0][i1];

                //(1, 1)
                T weight2_11 = 2 * partial1;
                weight2 = weight2_11;
                ans.wts2[DR1] += weight2;

                T weight3_prefac = 3*partial1;
                //(2, 1)
                T weight3_21 = weight3_prefac * E0;
                //(1, 2)
                T weight3_12 = weight3_prefac * E1;
                weight3 = weight3_12 + weight3_21;
                ans.wts3[DR1] += weight3;

                T weight4_prefac = 4*partial1;
                //(3, 1)
                T weight4_31 = weight4_prefac * sqE0;
                //(1, 3)
                T weight4_13 = weight4_prefac * sqE1;
                //(2, 2)
                T weight4_22 = 6*sqpartial1;
                weight4 = weight4_13 + weight4_31 + weight4_22;
                ans.wts4[DR1] += weight4;

                T weight5_prefac_1 = 5 * partial1;
                //(4, 1)
                T weight5_41 = weight5_prefac_1 * (sqE0 * E0);
                //(1, 4)
                T weight5_14 = weight5_prefac_1 * (sqE1 * E1);
                T weight5_prefac_2 = 10 * sqpartial1;
                //(3, 2)
                T weight5_32 = weight5_prefac_2 * E0;
                //(2, 3)
                T weight5_23 = weight5_prefac_2 * E1;
                weight5 = weight5_14 + weight5_41 + 
                          weight5_23 + weight5_32;
                ans.wts5[DR1] += weight5;

                T weight6_prefac_1 = 6 * partial1;
                //(5, 1)
                T weight6_51 = weight6_prefac_1 * (sqE0 * sqE0);
                //(1, 5)
                T weight6_15 = weight6_prefac_1 * (sqE1 * sqE1);
                T weight6_prefac_2 = 15 * sqpartial1;
                //(4, 2)
                T weight6_42 = weight6_prefac_2 * sqE0;
                //(2, 4)
                T weight6_24 = weight6_prefac_2 * sqE1;
                //(3, 3)
                T weight6_33 = 20 * sqpartial1 * partial1;
                weight6 = weight6_15 + weight6_51 + 
                          weight6_24 + weight6_42 + 
                          weight6_33;
                ans.wts6[DR1] += weight6;

                if constexpr (doPU){
                    isPU = isPU || PU->at(i1);
                    if(isPU){
                        ans.wts2_PU[DR1] += weight2;
                        ans.wts3_PU[DR1] += weight3;
                        ans.wts4_PU[DR1] += weight4;
                        ans.wts5_PU[DR1] += weight5;
                        ans.wts6_PU[DR1] += weight6;
                    }
                }

                if constexpr (doTransfer){
                    transfer6_2(i0, i1, DR1,
                                weight2_11,

                                weight3_21, 
                                weight3_12,

                                weight4_31,
                                weight4_22,
                                weight4_13,

                                weight5_41,
                                weight5_32,
                                weight5_23,
                                weight5_14,

                                weight6_51,
                                weight6_42,
                                weight6_33,
                                weight6_24,
                                weight6_15,

                                tin, ans);
                }

                for(unsigned i2=i1+1; i2<nPart; ++i2){
                    T E2 = Es[i2];
                    T sqE2 = square(E2);

                    T partial2 = partial1 * E2;

                    unsigned DR2 = max({DR1, dRs[i0][i2], 
                                             dRs[i1][i2]});

                    //(1, 1, 1)
                    weight3 = 6 * partial2;
                    ans.wts3[DR2] += weight3;

                    T weight4_prefac = 12 * partial2;
                    //(2, 1, 1)
                    T weight4_211 = weight4_prefac * E0;
                    //(1, 2, 1)
                    T weight4_121 = weight4_prefac * E1;
                    //(1, 1, 2)
                    T weight4_112 = weight4_prefac * E2;
                    weight4 = weight4_211 + 
                              weight4_121 +
                              weight4_112;
                    ans.wts4[DR2] += weight4;

                    T weight5_prefac_1 = 20 * partial2;
                    //(3, 1, 1)
                    T weight5_311 = weight5_prefac_1 * sqE0;
                    //(1, 3, 1)
                    T weight5_131 = weight5_prefac_1 * sqE1;
                    //(1, 1, 3)
                    T weight5_113 = weight5_prefac_1 * sqE2;
                    T weight5_prefac_2 = 30 * partial2;
                    //(2, 2, 1)
                    T weight5_221 = weight5_prefac_2 * (E0 * E1);
                    //(2, 1, 2)
                    T weight5_212 = weight5_prefac_2 * (E0 * E2);
                    //(1, 2, 2)
                    T weight5_122 = weight5_prefac_2 * (E1 * E2);
                    weight5 = weight5_311 + weight5_131 + 
                              weight5_113 + weight5_221 + 
                              weight5_212 + weight5_122;
                    ans.wts5[DR2] += weight5;

                    T weight6_prefac_1 = 30 * partial2;
                    //(4, 1, 1)
                    T weight6_411 = weight6_prefac_1 * (sqE0 * E0);
                    //(1, 4, 1)
                    T weight6_141 = weight6_prefac_1 * (sqE1 * E1);
                    //(1, 1, 4)
                    T weight6_114 = weight6_prefac_1 * (sqE2 * E2);
                    T weight6_prefac_2 = 60 * partial2;
                    //(3, 2, 1)
                    T weight6_321 = weight6_prefac_2 * (sqE0 * E1);
                    //(3, 1, 2)
                    T weight6_312 = weight6_prefac_2 * (sqE0 * E2);
                    //(1, 3, 2)
                    T weight6_132 = weight6_prefac_2 * (sqE1 * E2);
                    //(2, 3, 1)
                    T weight6_231 = weight6_prefac_2 * (E0 * sqE1);
                    //(2, 1, 3)
                    T weight6_213 = weight6_prefac_2 * (E0 * sqE2);
                    //(1, 2, 3)
                    T weight6_123 = weight6_prefac_2 * (E1 * sqE2);
                    //(2, 2, 2)
                    T weight6_222 = 90 * square(partial2);
                    weight6 = weight6_411 + weight6_141 + 
                              weight6_114 + weight6_321 + 
                              weight6_312 + weight6_132 + 
                              weight6_231 + weight6_213 + 
                              weight6_123 + weight6_222;
                    ans.wts6[DR2] += weight6;

                    if constexpr (doPU){
                        isPU = isPU || PU->at(i2);
                        if(isPU){
                            ans.wts3_PU[DR2] += weight3;
                            ans.wts4_PU[DR2] += weight4;
                            ans.wts5_PU[DR2] += weight5;
                            ans.wts6_PU[DR2] += weight6;
                        }
                    }

                    if constexpr (doTransfer){
                        transfer6_3(i0, i1, i2, DR2,
                                     weight3,

                                     weight4_211,
                                     weight4_121,
                                     weight4_112,

                                     weight5_311,
                                     weight5_131,
                                     weight5_113,
                                     weight5_221,
                                     weight5_212,
                                     weight5_122,

                                     weight6_411,
                                     weight6_141,
                                     weight6_114,
                                     weight6_321,
                                     weight6_312,
                                     weight6_231,
                                     weight6_132,
                                     weight6_213,
                                     weight6_123,
                                     weight6_222,

                                     tin, ans);
                    }

                    for(unsigned i3=i2+1; i3<nPart; ++i3){
                        T E3 = Es[i3];
                        T sqE3 = square(E3);

                        T partial3 = partial2 * E3;

                        unsigned DR3 = max({DR2, dRs[i0][i3], 
                                                 dRs[i1][i3], 
                                                 dRs[i2][i3]});

                        //(1, 1, 1, 1)
                        weight4 = 24 * partial3;
                        ans.wts4[DR3] += weight4;

                        T weight5_prefac = 60 * partial3;
                        //(2, 1, 1, 1)
                        T weight5_2111 = weight5_prefac * E0;
                        //(1, 2, 1, 1)
                        T weight5_1211 = weight5_prefac * E1;
                        //(1, 1, 2, 1)
                        T weight5_1121 = weight5_prefac * E2;
                        //(1, 1, 1, 2)
                        T weight5_1112 = weight5_prefac * E3;
                        weight5 = weight5_2111 + weight5_1211 + 
                                  weight5_1121 + weight5_1112;
                        ans.wts5[DR3] += weight5;

                        T weight6_prefac_1 = 120 * partial3;
                        //(3, 1, 1, 1)
                        T weight6_3111 = weight6_prefac_1 * sqE0;
                        //(1, 3, 1, 1)
                        T weight6_1311 = weight6_prefac_1 * sqE1;
                        //(1, 1, 3, 1)
                        T weight6_1131 = weight6_prefac_1 * sqE2;
                        //(1, 1, 1, 3)
                        T weight6_1113 = weight6_prefac_1 * sqE3;
                        T weight6_prefac_2 = 180 * partial3;
                        //(2, 2, 1, 1)
                        T weight6_2211 = weight6_prefac_2 * E0 * E1;
                        //(2, 1, 2, 1)
                        T weight6_2121 = weight6_prefac_2 * E0 * E2;
                        //(2, 1, 1, 2)
                        T weight6_2112 = weight6_prefac_2 * E0 * E3;
                        //(1, 2, 2, 1)
                        T weight6_1221 = weight6_prefac_2 * E1 * E2;
                        //(1, 2, 1, 2)
                        T weight6_1212 = weight6_prefac_2 * E1 * E3;
                        //(1, 1, 2, 2)
                        T weight6_1122 = weight6_prefac_2 * E2 * E3;
                        weight6 = weight6_3111 + weight6_1311 + 
                                  weight6_1131 + weight6_1113 + 
                                  weight6_2211 + weight6_2121 + 
                                  weight6_2112 + weight6_1221 + 
                                  weight6_1212 + weight6_1122;
                        ans.wts6[DR3] += weight6;

                        if constexpr (doPU){
                            isPU = isPU || PU->at(i3);
                            if(isPU){
                                ans.wts4_PU[DR3] += weight4;
                                ans.wts5_PU[DR3] += weight5;
                                ans.wts6_PU[DR3] += weight6;
                            }
                        }

                        if constexpr (doTransfer){
                            transfer6_4(i0, i1, i2, i3, DR3,

                                        weight4,

                                        weight5_2111,
                                        weight5_1211,
                                        weight5_1121,
                                        weight5_1112,

                                        weight6_3111,
                                        weight6_1311,
                                        weight6_1131,
                                        weight6_1113,
                                        weight6_2211,
                                        weight6_2121,
                                        weight6_2112,
                                        weight6_1221,
                                        weight6_1212,
                                        weight6_1122,

                                        tin, ans);
                        }

                        for(unsigned i4=i3+1; i4<nPart; ++i4){
                            T E4 = Es[i4];

                            T partial4 = partial3 * E4;

                            unsigned DR4 = max({DR3, dRs[i0][i4], 
                                                     dRs[i1][i4], 
                                                     dRs[i2][i4],
                                                     dRs[i3][i4]});

                            //(1, 1, 1, 1, 1)
                            weight5 = 120 * partial4;
                            ans.wts5[DR4] += weight5;

                            T weight6_prefac = 360 * partial4;
                            //(2, 1, 1, 1, 1)
                            T weight6_21111 = weight6_prefac * E0;
                            //(1, 2, 1, 1, 1)
                            T weight6_12111 = weight6_prefac * E1;
                            //(1, 1, 2, 1, 1)
                            T weight6_11211 = weight6_prefac * E2;
                            //(1, 1, 1, 2, 1)
                            T weight6_11121 = weight6_prefac * E3;
                            //(1, 1, 1, 1, 2)
                            T weight6_11112= weight6_prefac * E4;
                            weight6 = weight6_21111 + weight6_12111 + 
                                      weight6_11211 + weight6_11121 + 
                                      weight6_11112;
                            ans.wts6[DR4] += weight6;

                            if constexpr (doPU){
                                isPU = isPU || PU->at(i4);
                                if(isPU){
                                    ans.wts5_PU[DR4] += weight5;
                                    ans.wts6_PU[DR4] += weight6;
                                }
                            }

                            if constexpr (doTransfer){
                                transfer6_5(i0, i1, i2, i3, i4, DR4,
                                            
                                            weight5,

                                            weight6_21111,
                                            weight6_12111,
                                            weight6_11211,
                                            weight6_11121,
                                            weight6_11112,

                                            tin, ans);
                            }

                            for(unsigned i5=i4+1; i5<nPart; ++i5){
                                unsigned DR5 = max({DR4, dRs[i0][i5], 
                                                         dRs[i1][i5], 
                                                         dRs[i2][i5],
                                                         dRs[i3][i5],
                                                         dRs[i4][i5]});
                                //partition (1, 1, 1, 1, 1, 1)
                                weight6 = partial4 * Es[i5] * 720;
                                ans.wts6[DR5] += weight6;

                                if constexpr (doPU){
                                    isPU = isPU || PU->at(i5);
                                    if(isPU){
                                        ans.wts6_PU[DR5] += weight6;
                                    }
                                }

                                if constexpr (doTransfer){
                                    transfer6_6(i0, i1, i2, i3, i4, i5, DR5,
                                                
                                                weight6,

                                                tin, ans);
                                }
                            }
                        }
                    }
                }
            }
        }
    }   

    template <typename T, bool doPU, bool doTransfer, bool doRes3, bool doRes4, bool doRes4Fixed>
    result<T> fastEEC(result<T>& ans,

                      const jet& J, const axisptr& ax, 
                      const int order, const normType nt,

                      axisptr& coarseRLax,
                      axisptr& xiax,
                      axisptr& phiax,

                      axisptr& rax_dipole,
                      axisptr& ctax_dipole,

                      axisptr& rax_tee,
                      axisptr& ctax_tee,

                      axisptr& rax_triangle,
                      axisptr& ctax_triangle,
                      T shapetol,

                      const std::vector<bool>* const PU = nullptr,
                      const jet * const J_Reco = nullptr,
                      const arma::mat* ptrans = nullptr){
        assert(order >= 2 && order <= 6);

        umat dRs; 
        std::vector<T> Es;

        getDRs(dRs, J, ax);
        getEs<T>(Es, J, nt);

        unsigned NDR = histogram::axis::traits::extent(*ax);

        struct transferInputs<T> tin;
        if constexpr (doTransfer){
            getDRs(tin.dRs, *J_Reco, ax);
            tin.adj = adjacency(*ptrans);
            getPtrans<T>(tin.ptrans, *ptrans);

            getFloatDRs(tin.rin.floatDRs, *J_Reco);
            getEtasPhis(tin.rin.etas, tin.rin.phis, *J_Reco);

            tin.rin.coarseRL = coarseRLax;

            tin.rin.xi = xiax;
            tin.rin.phi = phiax;

            tin.rin.r_dipole = rax_dipole;
            tin.rin.ct_dipole = ctax_dipole;

            tin.rin.r_tee = rax_tee;
            tin.rin.ct_tee = ctax_tee;

            tin.rin.r_triangle = rax_triangle;
            tin.rin.ct_triangle = ctax_triangle;

            tin.rin.shapetol = shapetol;
        }

        struct resolvedInputs<T> rin;
        getFloatDRs(rin.floatDRs, J);
        getEtasPhis(rin.etas, rin.phis, J);
        rin.coarseRL = coarseRLax;

        rin.xi = xiax;
        rin.phi = phiax;

        rin.r_dipole = rax_dipole;
        rin.ct_dipole = ctax_dipole;

        rin.r_tee = rax_tee;
        rin.ct_tee = ctax_tee;

        rin.r_triangle = rax_triangle;
        rin.ct_triangle = ctax_triangle;

        rin.shapetol = shapetol;

        if (order == 2){
            start<T, doPU, doTransfer, 2, doRes3, doRes4, doRes4Fixed>(
                    dRs, Es, NDR, rin, ans, PU, &tin
            );
        } else if(order == 3){
            start<T, doPU, doTransfer, 3, doRes3, doRes4, doRes4Fixed>(
                    dRs, Es, NDR, rin, ans, PU, &tin
            );
        } else if(order == 4){
            start<T, doPU, doTransfer, 4, doRes3, doRes4, doRes4Fixed>(
                    dRs, Es, NDR, rin, ans, PU, &tin
            );
        } else if(order == 5){
            start<T, doPU, doTransfer, 5, doRes3, doRes4, doRes4Fixed>(
                    dRs, Es, NDR, rin, ans, PU, &tin
            );
        } else if(order == 6){
            start<T, doPU, doTransfer, 6, doRes3, doRes4, doRes4Fixed>(
                    dRs, Es, NDR, rin, ans, PU, &tin
            );
        }

        return ans;
    }
};

#endif
