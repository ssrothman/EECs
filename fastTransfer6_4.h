#ifndef EECS_FASTTRANSFER6_4_H
#define EECS_FASTTRANSFER6_4_H

namespace fastEEC{
    template <typename T>
    void transfer6_4(const unsigned i0,
                     const unsigned i1,
                     const unsigned i2,
                     const unsigned i3,

                     const unsigned DR_Gen,

                     const T weight4_1111,

                     const T weight5_2111,
                     const T weight5_1211,
                     const T weight5_1121,
                     const T weight5_1112,

                     const T weight6_3111,
                     const T weight6_1311,
                     const T weight6_1131,
                     const T weight6_1113,
                     const T weight6_2211,
                     const T weight6_2121,
                     const T weight6_2112,
                     const T weight6_1221,
                     const T weight6_1212,
                     const T weight6_1122,

                     const transferInputs<T>* tin,
                     result<T>& ans){

        uvec adj0 = tin->adj.at(i0);
        uvec adj1 = tin->adj.at(i1);
        uvec adj2 = tin->adj.at(i2);
        uvec adj3 = tin->adj.at(i3);

        for(const unsigned j0 : adj0){
            T partial0 = tin->ptrans[i0][j0];
            for(const unsigned j1 : adj1){
                T partial1 = partial0 * tin->ptrans[i1][j1];

                unsigned DR1_Reco = tin->dRs[j0][j1];

                for(const unsigned j2 : adj2){
                    T partial2 = partial1 * tin->ptrans[i2][j2];

                    unsigned DR2_Reco = max({DR1_Reco, tin->dRs[j0][j2],
                                                       tin->dRs[j1][j2]});

                    for(const unsigned j3 : adj3){
                        T partial3 = partial2 * tin->ptrans[i3][j3];
                        
                        unsigned DR3_Reco = max({DR2_Reco, tin->dRs[j0][j3],
                                                           tin->dRs[j1][j3],
                                                           tin->dRs[j2][j3]});

                        //transfer (1,1,1,1)
                        ans.transfer4[DR_Gen][DR3_Reco] += partial3 * weight4_1111;

                        for(const unsigned j4 : adj3){
                            T partial4 = partial3 * tin->ptrans[i3][j4];

                            unsigned DR4_Reco = max({DR3_Reco, tin->dRs[j0][j4],
                                                               tin->dRs[j1][j4],
                                                               tin->dRs[j2][j4],
                                                               tin->dRs[j3][j4]});
                            //transfer (1,1,1,2)
                            ans.transfer5[DR_Gen][DR4_Reco] += partial4 * weight5_1112;

                            for(const unsigned j5 : adj3){
                                T partial5 = partial4 * tin->ptrans[i3][j5];

                                unsigned DR5_Reco = max({DR4_Reco, tin->dRs[j0][j5],
                                                                   tin->dRs[j1][j5],
                                                                   tin->dRs[j2][j5],
                                                                   tin->dRs[j3][j5],
                                                                   tin->dRs[j4][j5]});
                                //transfer (1,1,1,3)
                                ans.transfer6[DR_Gen][DR5_Reco] += partial5 * weight6_1113;
                            }
                        }

                        for(const unsigned j4 : adj2){
                            T partial4 = partial3 * tin->ptrans[i2][j4];

                            unsigned DR4_Reco = max({DR3_Reco, tin->dRs[j0][j4],
                                                               tin->dRs[j1][j4],
                                                               tin->dRs[j2][j4],
                                                               tin->dRs[j3][j4]});
                            //transfer (1,1,2,1)
                            ans.transfer5[DR_Gen][DR4_Reco] += partial4 * weight5_1121;
                            
                            for(const unsigned j5 : adj3){
                                T partial5 = partial4 * tin->ptrans[i3][j5];

                                unsigned DR5_Reco = max({DR4_Reco, tin->dRs[j0][j5],
                                                                   tin->dRs[j1][j5],
                                                                   tin->dRs[j2][j5],
                                                                   tin->dRs[j3][j5],
                                                                   tin->dRs[j4][j5]});
                                //transfer (1,1,2,2)
                                ans.transfer6[DR_Gen][DR5_Reco] += partial5 * weight6_1122;
                            }

                            for(const unsigned j5 : adj2){
                                T partial5 = partial4 * tin->ptrans[i2][j5];

                                unsigned DR5_Reco = max({DR4_Reco, tin->dRs[j0][j5],
                                                                   tin->dRs[j1][j5],
                                                                   tin->dRs[j2][j5],
                                                                   tin->dRs[j3][j5],
                                                                   tin->dRs[j4][j5]});
                                //transfer (1,1,3,1)
                                ans.transfer6[DR_Gen][DR5_Reco] += partial5 * weight6_1131;
                            }
                        }

                        for(const unsigned j4 : adj1){
                            T partial4 = partial3 * tin->ptrans[i1][j4];

                            unsigned DR4_Reco = max({DR3_Reco, tin->dRs[j0][j4],
                                                               tin->dRs[j2][j4],
                                                               tin->dRs[j2][j4],
                                                               tin->dRs[j3][j4]});
                            //transfer (1,2,1,1)
                            ans.transfer5[DR_Gen][DR4_Reco] += partial4 * weight5_1211;

                            for(const unsigned j5 : adj3){
                                T partial5 = partial4 * tin->ptrans[i3][j5];

                                unsigned DR5_Reco = max({DR4_Reco, tin->dRs[j0][j5],
                                                                   tin->dRs[j1][j5],
                                                                   tin->dRs[j2][j5],
                                                                   tin->dRs[j3][j5],
                                                                   tin->dRs[j4][j5]});
                                //transfer (1,2,1,2)
                                ans.transfer6[DR_Gen][DR5_Reco] += partial5 * weight6_1212;
                            }

                            for(const unsigned j5 : adj2){
                                T partial5 = partial4 * tin->ptrans[i1][j5];

                                unsigned DR5_Reco = max({DR4_Reco, tin->dRs[j0][j5],
                                                                   tin->dRs[j1][j5],
                                                                   tin->dRs[j2][j5],
                                                                   tin->dRs[j3][j5],
                                                                   tin->dRs[j4][j5]});
                                //transfer (1,2,2,1)
                                ans.transfer6[DR_Gen][DR5_Reco] += partial5 * weight6_1221;
                            }

                            for(const unsigned j5 : adj0){
                                T partial5 = partial4 * tin->ptrans[i0][j5];

                                unsigned DR5_Reco = max({DR4_Reco, tin->dRs[j0][j5],
                                                                   tin->dRs[j1][j5],
                                                                   tin->dRs[j2][j5],
                                                                   tin->dRs[j3][j5],
                                                                   tin->dRs[j4][j5]});
                                //transfer (2,2,1,1)
                                ans.transfer6[DR_Gen][DR5_Reco] += partial5 * weight6_2211;
                            }

                            for(const unsigned j5 : adj1){
                                T partial5 = partial4 * tin->ptrans[i1][j5];

                                unsigned DR5_Reco = max({DR4_Reco, tin->dRs[j0][j5],
                                                                   tin->dRs[j1][j5],
                                                                   tin->dRs[j2][j5],
                                                                   tin->dRs[j3][j5],
                                                                   tin->dRs[j4][j5]});
                                //transfer (1,3,1,1)
                                ans.transfer6[DR_Gen][DR5_Reco] += partial5 * weight6_1311;
                            }
                        }

                        for(const unsigned j4 : adj0){
                            T partial4 = partial3 * tin->ptrans[i0][j4];

                            unsigned DR4_Reco = max({DR3_Reco, tin->dRs[j0][j4],
                                                               tin->dRs[j1][j4],
                                                               tin->dRs[j2][j4],
                                                               tin->dRs[j3][j4]});
                            //transfer (2,1,1,1)
                            ans.transfer5[DR_Gen][DR4_Reco] += partial4 * weight5_2111;

                            for(const unsigned j5 : adj2){
                                T partial5 = partial4 * tin->ptrans[i2][j5];

                                unsigned DR5_Reco = max({DR4_Reco, tin->dRs[j0][j5],
                                                                   tin->dRs[j1][j5],
                                                                   tin->dRs[j2][j5],
                                                                   tin->dRs[j3][j5],
                                                                   tin->dRs[j4][j5]});
                                //transfer (2,1,2,1)
                                ans.transfer6[DR_Gen][DR5_Reco] += partial5 * weight6_2121;
                            }

                            for(const unsigned j5 : adj3){
                                T partial5 = partial4 * tin->ptrans[i3][j5];

                                unsigned DR5_Reco = max({DR4_Reco, tin->dRs[j0][j5],
                                                                   tin->dRs[j1][j5],
                                                                   tin->dRs[j2][j5],
                                                                   tin->dRs[j3][j5],
                                                                   tin->dRs[j4][j5]});
                                //transfer (2,1,1,2)
                                ans.transfer6[DR_Gen][DR5_Reco] += partial5 * weight6_2112;
                            }

                            for(const unsigned j5 : adj0){
                                T partial5 = partial4 * tin->ptrans[i0][j5];

                                unsigned DR5_Reco = max({DR4_Reco, tin->dRs[j0][j5],
                                                                   tin->dRs[j1][j5],
                                                                   tin->dRs[j2][j5],
                                                                   tin->dRs[j3][j5],
                                                                   tin->dRs[j4][j5]});
                                //transfer (3,1,1,1)
                                ans.transfer6[DR_Gen][DR5_Reco] += partial5 * weight6_3111;
                            }
                        }
                    }
                }
            }
        }
    }
};

#endif
