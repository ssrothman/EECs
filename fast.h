#include <boost/histogram.hpp>
#include "boost/multi_array.hpp"

#include <cassert>

#include "SRothman/SimonTools/src/jets.h"

namespace fastEEC{
    using namespace boost;
    using namespace std;

    using axis_t = histogram::axis::variable<double>;
    using axisptr = std::shared_ptr<axis_t>;

    using umat = multi_array<unsigned, 2>;
    using umatptr = std::shared_ptr<umat>;
    
    template <typename T>
    struct result{
        vector<T> wts2;
        vector<T> wts3;
        vector<T> wts4;
        vector<T> wts5;
        vector<T> wts6;

        vector<T> wts2_PU;
        vector<T> wts3_PU;
        vector<T> wts4_PU;
        vector<T> wts5_PU;
        vector<T> wts6_PU;

        multi_array<T, 2> transfer2;
        multi_array<T, 2> transfer3;
        multi_array<T, 2> transfer4;
        multi_array<T, 2> transfer5;
        multi_array<T, 2> transfer6;
    };

    enum normType {
        RAWPT, 
        CORRPT,
        SUMPT, 
    };

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
    vector<T> getEs(const jet& J, const normType nt){
        vector<T> ans;
        ans.reserve(J.nPart);

        double normFact = getNormFact(J, nt);
        for(unsigned i=0; i<J.nPart; ++i){
            ans.push_back(J.particles.at(i).pt/normFact);
        }
        return ans;
    }

    umat getDRs(const jet& J, const axisptr& ax){
        umat ans(extents[J.nPart][J.nPart]);

        for (unsigned i0=0; i0<J.nPart; ++i0){
            for (unsigned i1=0; i1<J.nPart; ++i1){
                double deltaR = dR(J.particles[i0], J.particles[i1]);
                unsigned idx = static_cast<unsigned>(ax->index(deltaR) + 1);
                ans[i0][i1] = idx;
            }
        }
        return ans;
    }

    template <typename T>
    void fastTransfer2(const T wt,
                       const unsigned i0,
                       const unsigned i1,
                       const unsigned DR_Gen,
                       const unsigned maxOrder,
                       const umat* dRs_Reco, 
                       const adjacency* adj,
                       const multi_array<T, 2>* const ptrans,
                       result<T>& ans){
        const auto& adj0 = adj->at(i0);
        const auto& adj1 = adj->at(i1);

        T weight0, weight1;
        for (const auto& j0 : adj0){
            weight0 = wt * (*ptrans)[i0][j0];
            for(const auto& j1 : adj1){
                weight1 = weight0 * (*ptrans)[i1][j1];
                unsigned DR_Reco = (*dRs_Reco)[j0][j1];
                ans.transfer2[DR_Gen][DR_Reco] += weight1;
            }
        }
    }

    template <typename T>
    void fastTransfer3(const T wt,
                       const unsigned i0,
                       const unsigned i1,
                       const unsigned i2,
                       const unsigned DR_Gen,
                       const unsigned maxOrder,
                       const umat* dRs_Reco, 
                       const adjacency* adj,
                       const multi_array<T, 2>* const ptrans,
                       result<T>& ans){
        const auto& adj0 = adj->at(i0);
        const auto& adj1 = adj->at(i1);
        const auto& adj2 = adj->at(i2);

        T weight0, weight1, weight2;
        for (const auto& j0 : adj0){
            weight0 = wt * (*ptrans)[i0][j0];
            for(const auto& j1 : adj1){
                weight1 = weight0 * (*ptrans)[i1][j1];
                unsigned DR1 = (*dRs_Reco)[j0][j1];
                for(const auto& j2 : adj2){
                    weight2 = weight1 * (*ptrans)[i2][j2];
                    unsigned DR_Reco = max(DR1, (*dRs_Reco)[j0][j2],
                                                (*dRs_Reco)[j1][j2]);
                    ans.transfer3[DR_Gen][DR_Reco] += weight2;
                }
            }
        }
    }

    template <typename T>
    void fastTransfer4(const T wt,
                       const unsigned i0,
                       const unsigned i1,
                       const unsigned i2,
                       const unsigned i3,
                       const unsigned DR_Gen,
                       const unsigned maxOrder,
                       const umat* dRs_Reco, 
                       const adjacency* adj,
                       const multi_array<T, 2>* const ptrans,
                       result<T>& ans){
        const auto& adj0 = adj->at(i0);
        const auto& adj1 = adj->at(i1);
        const auto& adj2 = adj->at(i2);
        const auto& adj3 = adj->at(i3);

        T weight0, weight1, weight2, weight3;
        for (const auto& j0 : adj0){
            weight0 = wt * (*ptrans)[i0][j0];
            for(const auto& j1 : adj1){
                weight1 = weight0 * (*ptrans)[i1][j1];
                unsigned DR1 = (*dRs_Reco)[j0][j1];
                for(const auto& j2 : adj2){
                    weight2 = weight1 * (*ptrans)[i2][j2];
                    unsigned DR2 = max(DR1, (*dRs_Reco)[j0][j2],
                                            (*dRs_Reco)[j1][j2]);
                    for(const auto& j3 : adj3){
                        weight3 = weight2 * (*ptrans)[i3][j3];
                        unsigned DR_Reco = max(DR2, (*dRs_Reco)[j0][j3],
                                                    (*dRs_Reco)[j1][j3],
                                                    (*dRs_Reco)[j2][j3]);
                        ans.transfer4[DR_Gen][DR_Reco] += weight3;
                    }
                }
            }
        }
    }

    template <typename T>
    void fastTransfer5(const T wt,
                       const unsigned i0,
                       const unsigned i1,
                       const unsigned i2,
                       const unsigned i3,
                       const unsigned i4,
                       const unsigned DR_Gen,
                       const unsigned maxOrder,
                       const umat* dRs_Reco, 
                       const adjacency* adj,
                       const multi_array<T, 2>* const ptrans,
                       result<T>& ans){
        const auto& adj0 = adj->at(i0);
        const auto& adj1 = adj->at(i1);
        const auto& adj2 = adj->at(i2);
        const auto& adj3 = adj->at(i3);
        const auto& adj4 = adj->at(i4);

        T weight0, weight1, weight2, weight3, weight4;
        for (const auto& j0 : adj0){
            weight0 = wt * (*ptrans)[i0][j0];
            for(const auto& j1 : adj1){
                weight1 = weight0 * (*ptrans)[i1][j1];
                unsigned DR1 = (*dRs_Reco)[j0][j1];
                for(const auto& j2 : adj2){
                    weight2 = weight1 * (*ptrans)[i2][j2];
                    unsigned DR2 = max(DR1, (*dRs_Reco)[j0][j2],
                                            (*dRs_Reco)[j1][j2]);
                    for(const auto& j3 : adj3){
                        weight3 = weight2 * (*ptrans)[i3][j3];
                        unsigned DR3 = max(DR2, (*dRs_Reco)[j0][j3],
                                                (*dRs_Reco)[j1][j3],
                                                (*dRs_Reco)[j2][j3]);
                        for(const auto& j4 : adj4){
                            weight4 = weight3 * (*ptrans)[i4][j4];
                            unsigned DR_Reco = max(DR3, (*dRs_Reco)[j0][j4],
                                                        (*dRs_Reco)[j1][j4],
                                                        (*dRs_Reco)[j2][j4],
                                                        (*dRs_Reco)[j3][j4]);
                            ans.transfer5[DR_Gen][DR_Reco] += weight4;
                        }
                    }
                }
            }
        }
    }

    template <typename T>
    void fastTransfer6(const T wt,
                       const unsigned i0,
                       const unsigned i1,
                       const unsigned i2,
                       const unsigned i3,
                       const unsigned i4,
                       const unsigned i5,
                       const unsigned DR_Gen,
                       const unsigned maxOrder,
                       const umat* dRs_Reco, 
                       const adjacency* adj,
                       const multi_array<T, 2>* const ptrans,
                       result<T>& ans){
        const auto& adj0 = adj->at(i0);
        const auto& adj1 = adj->at(i1);
        const auto& adj2 = adj->at(i2);
        const auto& adj3 = adj->at(i3);
        const auto& adj4 = adj->at(i4);
        const auto& adj5 = adj->at(i5);

        T weight0, weight1, weight2, weight3, weight4, weight5;
        for (const auto& j0 : adj0){
            weight0 = wt * (*ptrans)[i0][j0];
            for(const auto& j1 : adj1){
                weight1 = weight0 * (*ptrans)[i1][j1];
                unsigned DR1 = (*dRs_Reco)[j0][j1];
                for(const auto& j2 : adj2){
                    weight2 = weight1 * (*ptrans)[i2][j2];
                    unsigned DR2 = max(DR1, (*dRs_Reco)[j0][j2],
                                            (*dRs_Reco)[j1][j2]);
                    for(const auto& j3 : adj3){
                        weight3 = weight2 * (*ptrans)[i3][j3];
                        unsigned DR3 = max(DR2, (*dRs_Reco)[j0][j3],
                                                (*dRs_Reco)[j1][j3],
                                                (*dRs_Reco)[j2][j3]);
                        for(const auto& j4 : adj4){
                            weight4 = weight3 * (*ptrans)[i4][j4];
                            unsigned DR4 = max(DR3, (*dRs_Reco)[j0][j4],
                                                    (*dRs_Reco)[j1][j4],
                                                    (*dRs_Reco)[j2][j4],
                                                    (*dRs_Reco)[j3][j4]);
                            for(const auto& j5 : adj5){
                                weight5 = weight4 * (*ptrans)[i5][j5];
                                unsigned DR_Reco = max(DR4, (*dRs_Reco)[j0][j5],
                                                            (*dRs_Reco)[j1][j5],
                                                            (*dRs_Reco)[j2][j5],
                                                            (*dRs_Reco)[j3][j5],
                                                            (*dRs_Reco)[j4][j5]);
                                ans.transfer6[DR_Gen][DR_Reco] += weight5;
                            }
                        }
                    }
                }
            }
        }
    }

    /*
     * Idea for tomorrow
     *
     *  1. implement helper functions for each partition
     *      eg parittion (4,2) for that particular case of 6th order 
     *  2. do transfer inside the helper functions,
     *      a) compute prefactor 
     *      b) handle each permutation separately 
     *      c) test the overhead of this when not doing transfer
     *  3. need to think a bit more about the resolved case. 
     *      Should have the same amount of symmetry as transfer calculations
     *      So maybe we can also do an if constexpr in the helpers? 
     *      Need to make sure everything has a nice uniform interface
     */

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

        double weight = 1;
        bool isPU=false;

        for(unsigned i0=0; i0<nPart; ++i0){
            double partial0 = Es[i0];
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

        double weight2, weight3;
        bool isPU=false;
        for(unsigned i0=0; i0<nPart; ++i0){
            double E0 = Es[i0];
            double partial0 = E0;

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
                double E1 = Es[i1];
                double partial1 = partial0 * E1;

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

        double weight2, weight3, weight4;
        bool isPU=false;

        for(unsigned i0=0; i0<nPart; ++i0){
            double E0 = Es[i0];
            double sqE0 = square(E0);
            double partial0 = E0;

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
                double E1 = Es[i1];
                double partial1 = partial0 * E1;
                
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
                    double E2 = Es[i2];
                    double partial2 = partial1 * E2;

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

        double weight2, weight3, weight4, weight5;
        bool isPU=false;

        for(unsigned i0=0; i0<nPart; ++i0){
            double E0 = Es[i0];
            double sqE0 = square(E0);

            double partial0 = E0;

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
                double E1 = Es[i1];
                double sqE1 = square(E1);

                double partial1 = partial0 * E1;

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
                    double E2 = Es[i2];
                    double sqE2 = square(E2);

                    double partial2 = partial1 * E2;

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
                        double E3 = Es[i3];

                        double partial3 = partial2 * E3;

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
 
    template <typename T, bool doPU>
    void fastSixthOrder(const umat& dRs,
                        const vector<T>& Es,
                        const unsigned NDR,
                        result<T>& ans,
                         const vector<bool>* const PU = nullptr){
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

        const unsigned nPart = Es.size();

        double weight2, weight3, weight4, weight5, weight6;
        bool isPU=false;

        for(unsigned i0=0; i0<nPart; ++i0){
            double E0 = Es[i0]; 
            double sqE0 = square(E0);

            double partial0 = E0;

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
                if(PU->at(i0)){
                    ans.wts2_PU[DR0] += weight2;
                    ans.wts3_PU[DR0] += weight3;
                    ans.wts4_PU[DR0] += weight4;
                    ans.wts5_PU[DR0] += weight5;
                    ans.wts6_PU[DR0] += weight6;
                }
            }

            for(unsigned i1=i0+1; i1<nPart; ++i1){
                double E1 = Es[i1];
                double sqE1 = square(E1);

                double partial1 = partial0 * E1;
                double sqpartial1 = square(partial1);

                unsigned DR1 = dRs[i0][i1];

                //partition (1,1)
                weight2 = 2 * partial1;
                ans.wts2[DR1] += weight2;
                //partition (2, 1)
                weight3 = 3 * partial1 * (E0 + E1);
                ans.wts3[DR1] += weight3;
                //partition (2, 2)
                weight4 = 6 * sqpartial1;
                //partition (3, 1)
                weight4 += 4 * partial1 * (sqE0 + sqE1);
                ans.wts4[DR1] += weight4;
                //partition (4, 1)
                weight5 = 5 * partial1 * (sqE1 * E1 +
                                          sqE0 * E0);
                //partition (3, 2)
                weight5 += 10 * sqpartial1 * (E0 + E1);
                ans.wts5[DR1] += weight5;

                //partition (5, 1)
                weight6 = 6 * partial1 * (sqE1 * sqE1 + 
                                          sqE0 * sqE0);
                //partition (4, 2)
                weight6 += 15 * sqpartial1 * (sqE0 + sqE1);
                //partition (3,3)
                weight6 += 20 * sqpartial1 * partial1;
                ans.wts6[DR1] += weight6;

                if constexpr (doPU){
                    if(PU->at(i1)){
                        ans.wts2_PU[DR1] += weight2;
                        ans.wts3_PU[DR1] += weight3;
                        ans.wts4_PU[DR1] += weight4;
                        ans.wts5_PU[DR1] += weight5;
                        ans.wts6_PU[DR1] += weight6;
                    }
                }

                for(unsigned i2=i1+1; i2<nPart; ++i2){
                    double E2 = Es[i2];
                    double sqE2 = square(E2);

                    double partial2 = partial1 * E2;

                    unsigned DR2 = max({DR1, dRs[i0][i2], 
                                             dRs[i1][i2]});

                    //partition (1,1,1)
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
                    //partition (4, 1, 1)
                    weight6 = 30 * partial2 * (sqE0 * E0 + 
                                              sqE1 * E1 + 
                                              sqE2 * E2);
                    //partition (3, 2, 1)
                    weight6 += 60 * partial2 * (sqE0 * E1 +
                                               sqE0 * E2 + 
                                               sqE1 * E2 +
                                               sqE1 * E0 +
                                               sqE2 * E0 +
                                               sqE2 * E1);
                    //partition (2, 2, 2)
                    weight6 += 90 * square(partial2);
                    ans.wts6[DR2] += weight6;

                    if constexpr (doPU){
                        if(PU->at(i2)){
                            ans.wts3_PU[DR2] += weight3;
                            ans.wts4_PU[DR2] += weight4;
                            ans.wts5_PU[DR2] += weight5;
                            ans.wts6_PU[DR2] += weight6;
                        }
                    }

                    for(unsigned i3=i2+1; i3<nPart; ++i3){
                        double E3 = Es[i3];
                        double sqE3 = square(E3);

                        double partial3 = partial2 * E3;

                        unsigned DR3 = max({DR2, dRs[i0][i3], 
                                                 dRs[i1][i3], 
                                                 dRs[i2][i3]});

                        //partition (1, 1, 1, 1)
                        weight4 = 24 * partial3;
                        ans.wts4[DR3] += weight4;
                        //partition (2, 1, 1, 1)
                        weight5 = 60 * partial3 * (E0 + E1 + E2 + E3);
                        ans.wts5[DR3] += weight5;
                        //partition (3, 1, 1, 1)
                        weight6 = 120 * partial3 * (sqE0 +
                                                    sqE1 +
                                                    sqE2 +
                                                    sqE3);
                        //partition (2, 2, 1, 1)
                        weight6 += 180 * partial3 * (E0 * E1 +
                                                     E0 * E2 +
                                                     E0 * E3 +
                                                     E1 * E2 +
                                                     E1 * E3 +
                                                     E2 * E3);
                        ans.wts6[DR3] += weight6;

                        if constexpr (doPU){
                            if(PU->at(i3)){
                                ans.wts4_PU[DR3] += weight4;
                                ans.wts5_PU[DR3] += weight5;
                                ans.wts6_PU[DR3] += weight6;
                            }
                        }

                        for(unsigned i4=i3+1; i4<nPart; ++i4){
                            double E4 = Es[i4];

                            double partial4 = partial3 * E4;

                            unsigned DR4 = max({DR3, dRs[i0][i4], 
                                                     dRs[i1][i4], 
                                                     dRs[i2][i4],
                                                     dRs[i3][i4]});

                            //partition (1, 1, 1, 1, 1)
                            weight5 = 120 * partial4;
                            ans.wts5[DR4] += weight5;

                            //partition (2, 1, 1, 1, 1)
                            weight6 = 360 * partial4 * (E0 + 
                                                        E1 + 
                                                        E2 + 
                                                        E3 + 
                                                        E4);
                            ans.wts6[DR4] += weight6;

                            if constexpr (doPU){
                                if(PU->at(i4)){
                                    ans.wts5_PU[DR4] += weight5;
                                    ans.wts6_PU[DR4] += weight6;
                                }
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
                                    if(PU->at(i5)){
                                        ans.wts6_PU[DR5] += weight6;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }   

    template <typename T, bool doPU>
    result<T> fastEEC(const jet& J, const axisptr& ax, 
                      const int order, const normType nt,
                      const std::vector<bool>* const PU = nullptr){

        assert(order >= 2 && order <= 6);

        result<T> ans;

        umat dRs = getDRs(J, ax);
        std::vector<T> Es = getEs<T>(J, nt);
        unsigned NDR = histogram::axis::traits::extent(*ax);

        if(order == 2){
            fastSecondOrder<T, doPU>(dRs, Es, NDR, ans, PU); 
        } else if(order == 3){
            fastThirdOrder<T, doPU>(dRs, Es, NDR, ans, PU);
        } else if(order == 4){
            fastFourthOrder<T, doPU>(dRs, Es, NDR, ans, PU);
        } else if(order == 5){
            fastFifthOrder<T, doPU>(dRs, Es, NDR, ans, PU);
        } else if(order == 6){
            fastSixthOrder<T, doPU>(dRs, Es, NDR, ans, PU);
        }
        return ans;
    }
};
