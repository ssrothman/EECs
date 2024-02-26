#ifndef EECS_FASTTRANSFER6_3_H
#define EECS_FASTTRANSFER6_3_H

namespace fastEEC{
    template <typename T>
    void transfer6_3(const unsigned i0,
                     const unsigned i1,
                     const unsigned i2,

                     const unsigned DR_Gen,

                     const T weight3_111,

                     const T weight4_211,
                     const T weight4_121,
                     const T weight4_112,

                     const T weight5_311,
                     const T weight5_131,
                     const T weight5_113,
                     const T weight5_221,
                     const T weight5_212,
                     const T weight5_122,

                     const T weight6_411,
                     const T weight6_141,
                     const T weight6_114,
                     const T weight6_321,
                     const T weight6_312,
                     const T weight6_231,
                     const T weight6_132,
                     const T weight6_213,
                     const T weight6_123,
                     const T weight6_222,

                     const transferInputs<T>* tin,
                     result<T>& ans){

        const uvec& adj0 = tin->adj.at(i0);
        const uvec& adj1 = tin->adj.at(i1);
        const uvec& adj2 = tin->adj.at(i2);

        for(const unsigned j0 : adj0){
            for(const unsigned j1 : adj1){
                unsigned DR1_Reco = tin->dRs[j0][j1];
                for(const unsigned j2 : adj2){
                    T partial2 = tin->ptrans[i0][j0] * tin->ptrans[i1][j1] * tin->ptrans[i2][j2]; 

                    unsigned DR2_Reco = max({DR1_Reco, tin->dRs[j0][j2],
                                                       tin->dRs[j1][j2]});

                    //transfer (1,1,1)
                    ans.transfer3[DR_Gen][DR2_Reco] += partial2 * weight3_111;

                    for(const unsigned j3 : adj2){
                        T partial3 = partial2 * tin->ptrans[i2][j3];

                        unsigned DR3_Reco = max({DR2_Reco, tin->dRs[j0][j3],
                                                           tin->dRs[j1][j3],
                                                           tin->dRs[j2][j3]});
                        //transfer (1,1,2)
                        ans.transfer4[DR_Gen][DR3_Reco] += partial3 * weight4_112;

                        for(const unsigned j4 : adj2){
                            T partial4 = partial3 * tin->ptrans[i2][j4];

                            unsigned DR4_Reco = max({DR3_Reco, tin->dRs[j0][j4],
                                                               tin->dRs[j1][j4],
                                                               tin->dRs[j2][j4],
                                                               tin->dRs[j3][j4]});
                            //transfer (1,1,3)
                            ans.transfer5[DR_Gen][DR4_Reco] += partial4 * weight5_113;

                            for(const unsigned j5 : adj2){
                                T partial5 = partial4 * tin->ptrans[i2][j5];

                                unsigned DR5_Reco = max({DR4_Reco, tin->dRs[j0][j5],
                                                                   tin->dRs[j1][j5],
                                                                   tin->dRs[j2][j5],
                                                                   tin->dRs[j3][j5],
                                                                   tin->dRs[j4][j5]});
                                //transfer (1,1,4)
                                ans.transfer6[DR_Gen][DR5_Reco] += partial5 * weight6_114;
                            }

                            for(const unsigned j5 : adj0){
                                T partial5 = partial4 * tin->ptrans[i0][j5];

                                unsigned DR5_Reco = max({DR4_Reco, tin->dRs[j0][j5],
                                                                   tin->dRs[j1][j5],
                                                                   tin->dRs[j2][j5],
                                                                   tin->dRs[j3][j5],
                                                                   tin->dRs[j4][j5]});
                                //transfer (2,1,3)
                                ans.transfer6[DR_Gen][DR5_Reco] += partial5 * weight6_213;
                            }
                        }

                        for(const unsigned j4 : adj1){
                            T partial4 = partial3 * tin->ptrans[i1][j4];

                            unsigned DR4_Reco = max({DR3_Reco, tin->dRs[j0][j4],
                                                               tin->dRs[j1][j4],
                                                               tin->dRs[j2][j4],
                                                               tin->dRs[j3][j4]});
                            //transfer (1,2,2)
                            ans.transfer5[DR_Gen][DR4_Reco] += partial4 * weight5_122;

                            for(const unsigned j5 : adj2){
                                T partial5 = partial4 * tin->ptrans[i2][j5];

                                unsigned DR5_Reco = max({DR4_Reco, tin->dRs[j0][j5],
                                                                   tin->dRs[j1][j5],
                                                                   tin->dRs[j2][j5],
                                                                   tin->dRs[j3][j5],
                                                                   tin->dRs[j4][j5]});
                                //transfer (1,2,3)
                                ans.transfer6[DR_Gen][DR5_Reco] += partial5 * weight6_123;
                            }

                            for(const unsigned j5 : adj1){
                                T partial5 = partial4 * tin->ptrans[i1][j5];

                                unsigned DR5_Reco = max({DR4_Reco, tin->dRs[j0][j5],
                                                                   tin->dRs[j1][j5],
                                                                   tin->dRs[j2][j5],
                                                                   tin->dRs[j3][j5],
                                                                   tin->dRs[j4][j5]});
                                //transfer (1,3,2)
                                ans.transfer6[DR_Gen][DR5_Reco] += partial5 * weight6_132;
                            }

                            for(const unsigned j5 : adj0){
                                T partial5 = partial4 * tin->ptrans[i0][j5];

                                unsigned DR5_Reco = max({DR4_Reco, tin->dRs[j0][j5],
                                                                   tin->dRs[j1][j5],
                                                                   tin->dRs[j2][j5],
                                                                   tin->dRs[j3][j5],
                                                                   tin->dRs[j4][j5]});

                                //transfer (2,2,2)
                                ans.transfer6[DR_Gen][DR5_Reco] += partial5 * weight6_222;
                            }
                        }

                        for(const unsigned j4 : adj0){
                            T partial4 = partial3 * tin->ptrans[i0][j4];

                            unsigned DR4_Reco = max({DR3_Reco, tin->dRs[j0][j4],
                                                               tin->dRs[j1][j4],
                                                               tin->dRs[j2][j4],
                                                               tin->dRs[j3][j4]});
                            //transfer (2,1,2)
                            ans.transfer5[DR_Gen][DR4_Reco] += partial4 * weight5_212;

                            for(const unsigned j5 : adj0){
                                T partial5 = partial4 * tin->ptrans[i0][j5];

                                unsigned DR5_Reco = max({DR4_Reco, tin->dRs[j0][j5],
                                                                   tin->dRs[j1][j5],
                                                                   tin->dRs[j2][j5],
                                                                   tin->dRs[j3][j5],
                                                                   tin->dRs[j4][j5]});
                                //transfer (3,1,2)
                                ans.transfer6[DR_Gen][DR5_Reco] += partial5 * weight6_312;
                            }
                        }
                    }

                    for(const unsigned j3 : adj1){
                        T partial3 = partial2 * tin->ptrans[i1][j3];

                        unsigned DR3_Reco = max({DR2_Reco, tin->dRs[j0][j3],
                                                           tin->dRs[j1][j3],
                                                           tin->dRs[j2][j3]});
                        //transfer (1,2,1)
                        ans.transfer4[DR_Gen][DR3_Reco] += partial3 * weight4_121;

                        for(const unsigned j4 : adj1){
                            T partial4 = partial3 * tin->ptrans[i1][j4];

                            unsigned DR4_Reco = max({DR3_Reco, tin->dRs[j0][j4],
                                                               tin->dRs[j1][j4],
                                                               tin->dRs[j2][j4],
                                                               tin->dRs[j3][j4]});
                            //transfer (1,3,1)
                            ans.transfer5[DR_Gen][DR4_Reco] += partial4 * weight5_131;

                            for(const unsigned j5 : adj1){
                                T partial5 = partial4 * tin->ptrans[i1][j5];

                                unsigned DR5_Reco = max({DR4_Reco, tin->dRs[j0][j5],
                                                                   tin->dRs[j1][j5],
                                                                   tin->dRs[j2][j5],
                                                                   tin->dRs[j3][j5],
                                                                   tin->dRs[j4][j5]});
                                //transfer (1,4,1)
                                ans.transfer6[DR_Gen][DR5_Reco] += partial5 * weight6_141;
                            }

                            for(const unsigned j5 : adj0){
                                T partial5 = partial4 * tin->ptrans[i0][j5];

                                unsigned DR5_Reco = max({DR4_Reco, tin->dRs[j0][j5],
                                                                   tin->dRs[j1][j5],
                                                                   tin->dRs[j2][j5],
                                                                   tin->dRs[j3][j5],
                                                                   tin->dRs[j4][j5]});
                                //transfer (2,3,1)
                                ans.transfer6[DR_Gen][DR5_Reco] += partial5 * weight6_231;
                            }
                        }

                        for(const unsigned j4 : adj0){
                            T partial4 = partial3 * tin->ptrans[i0][j4];

                            unsigned DR4_Reco = max({DR3_Reco, tin->dRs[j0][j4],
                                                               tin->dRs[j1][j4],
                                                               tin->dRs[j2][j4],
                                                               tin->dRs[j3][j4]});
                            //transfer (2,2,1)
                            ans.transfer5[DR_Gen][DR4_Reco] += partial4 * weight5_221;

                            for(const unsigned j5 : adj0){
                                T partial5 = partial4 * tin->ptrans[i0][j5];

                                unsigned DR5_Reco = max({DR4_Reco, tin->dRs[j0][j5],
                                                                   tin->dRs[j1][j5],
                                                                   tin->dRs[j2][j5],
                                                                   tin->dRs[j3][j5],
                                                                   tin->dRs[j4][j5]});
                                //transfer (3,2,1)
                                ans.transfer6[DR_Gen][DR5_Reco] += partial5 * weight6_321;
                            }
                        }
                    }

                    for(const unsigned j3 : adj0){
                        T partial3 = partial2 * tin->ptrans[i0][j3];

                        unsigned DR3_Reco = max({DR2_Reco, tin->dRs[j0][j3],
                                                           tin->dRs[j1][j3],
                                                           tin->dRs[j2][j3]});

                        //transfer (2,1,1)
                        ans.transfer4[DR_Gen][DR3_Reco] += partial3 * weight4_211;

                        for(const unsigned j4 : adj0){
                            T partial4 = partial3 * tin->ptrans[i0][j4];

                            unsigned DR4_Reco = max({DR3_Reco, tin->dRs[j0][j4],
                                                               tin->dRs[j1][j4],
                                                               tin->dRs[j2][j4],
                                                               tin->dRs[j3][j4]});
                            //transfer (3,1,1)
                            ans.transfer5[DR_Gen][DR4_Reco] += partial4 * weight5_311;

                            for(const unsigned j5 : adj0){
                                T partial5 = partial4 * tin->ptrans[i0][j5];

                                unsigned DR5_Reco = max({DR4_Reco, tin->dRs[j0][j5],
                                                                   tin->dRs[j1][j5],
                                                                   tin->dRs[j2][j5],
                                                                   tin->dRs[j3][j5],
                                                                   tin->dRs[j4][j5]});
                                //transfer (4,1,1)
                                ans.transfer6[DR_Gen][DR5_Reco] += partial5 * weight6_411;
                            }
                        }
                    }
                }
            }
        }
    }
};

#endif
