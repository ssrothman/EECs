#ifndef EECS_FASTTRANSFER6_2_H
#define EECS_FASTTRANSFER6_2_H

namespace fastEEC{
    template <typename T>
    void transfer6_2(const unsigned i0,
                     const unsigned i1,

                     const unsigned DR_Gen,

                     const T weight2_11, 
                     
                     const T weight3_21, 
                     const T weight3_12,

                     const T weight4_31, 
                     const T weight4_22,
                     const T weight4_13,

                     const T weight5_41, 
                     const T weight5_32, 
                     const T weight5_23, 
                     const T weight5_14,

                     const T weight6_51, 
                     const T weight6_42,
                     const T weight6_33,
                     const T weight6_24,
                     const T weight6_15,

                     const transferInputs<T>* tin,
                     result<T>& ans){

        const uvec& adj0 = tin->adj.at(i0);
        const uvec& adj1 = tin->adj.at(i1);

        for(const unsigned j0 : adj0){
            T partial0 = tin->ptrans[i0][j0];
            for (const unsigned j1 : adj1){
                T partial1 = partial0 * tin->ptrans[i1][j1];

                unsigned DR1_Reco = tin->dRs[j0][j1];

                //transfer (1,1)
                ans.transfer2[DR_Gen][DR1_Reco] += partial1 * weight2_11;

                for (const unsigned j2 : adj1){
                    T partial2 = partial1 * tin->ptrans[i1][j2];

                    unsigned DR2_Reco = max({DR1_Reco, tin->dRs[j0][j2],
                                                       tin->dRs[j1][j2]});
                    //transfer (1,2)
                    ans.transfer3[DR_Gen][DR2_Reco] += partial2 * weight3_12;

                    for(const unsigned j3 : adj1){
                        T partial3 = partial2 * tin->ptrans[i1][j3];

                        unsigned DR3_Reco = max({DR2_Reco, tin->dRs[j0][j3],
                                                           tin->dRs[j1][j3],
                                                           tin->dRs[j2][j3]});
                        //transfer (1,3)
                        ans.transfer4[DR_Gen][DR3_Reco] += partial3 * weight4_13;

                        for(const unsigned j4 : adj1){
                            T partial4 = partial3 * tin->ptrans[i1][j4];

                            unsigned DR4_Reco = max({DR3_Reco, tin->dRs[j0][j4],
                                                               tin->dRs[j1][j4],
                                                               tin->dRs[j2][j4],
                                                               tin->dRs[j3][j4]});

                            //transfer (1,4)
                            ans.transfer5[DR_Gen][DR4_Reco] += partial4 * weight5_14;

                            for(const unsigned j5 : adj1){
                                T partial5 = partial4 * tin->ptrans[i1][j5];

                                unsigned DR5_Reco = max({DR4_Reco, tin->dRs[j0][j5],
                                                                   tin->dRs[j1][j5],
                                                                   tin->dRs[j2][j5],
                                                                   tin->dRs[j3][j5],
                                                                   tin->dRs[j4][j5]});
                                //transfer (1,5)
                                ans.transfer6[DR_Gen][DR5_Reco] += partial5 * weight6_15;
                            }

                            for(const unsigned j5 : adj0){
                                T partial5 = partial4 * tin->ptrans[i0][j5];

                                unsigned DR5_Reco = max({DR4_Reco, tin->dRs[j0][j5],
                                                                   tin->dRs[j1][j5],
                                                                   tin->dRs[j2][j5],
                                                                   tin->dRs[j3][j5],
                                                                   tin->dRs[j4][j5]});
                                //transfer (2,4)
                                ans.transfer6[DR_Gen][DR5_Reco] += partial5 * weight6_24;
                            }
                        }

                        for(const unsigned j4 : adj0){
                            T partial4 = partial3 * tin->ptrans[i0][j4];

                            unsigned DR4_Reco = max({DR3_Reco, tin->dRs[j0][j4],
                                                               tin->dRs[j1][j4],
                                                               tin->dRs[j2][j4],
                                                               tin->dRs[j3][j4]});
                            //transfer (2,3)
                            ans.transfer5[DR_Gen][DR4_Reco] += partial4 * weight5_23;

                            for(const unsigned j5 : adj0){
                                T partial5 = partial4 * tin->ptrans[i0][j5];

                                unsigned DR5_Reco = max({DR4_Reco, tin->dRs[j0][j5],
                                                                   tin->dRs[j1][j5],
                                                                   tin->dRs[j2][j5],
                                                                   tin->dRs[j3][j5],
                                                                   tin->dRs[j4][j5]});
                                //transfer (3,3)
                                ans.transfer6[DR_Gen][DR5_Reco] += partial5 * weight6_33;
                            }
                        }
                    }

                    for(const unsigned j3:adj0){
                        T partial3 = partial2 * tin->ptrans[i0][j3];

                        unsigned DR3_Reco = max({DR2_Reco, tin->dRs[j0][j3],
                                                           tin->dRs[j1][j3],
                                                           tin->dRs[j2][j3]});
                        //transfer (2,2)
                        ans.transfer4[DR_Gen][DR3_Reco] += partial3 * weight4_22;

                        for(const unsigned j4:adj0){
                            T partial4 = partial2 * tin->ptrans[i0][j4];

                            unsigned DR4_Reco = max({DR3_Reco, tin->dRs[j0][j4],
                                                               tin->dRs[j1][j4],
                                                               tin->dRs[j2][j4],
                                                               tin->dRs[j3][j4]});
                            //transfer (3,2)
                            ans.transfer5[DR_Gen][DR4_Reco] += partial4 * weight5_32;

                            for(const unsigned j5:adj0){
                                T partial5 = partial2 * tin->ptrans[i0][j5];

                                unsigned DR5_Reco = max({DR4_Reco, tin->dRs[j0][j5],
                                                                   tin->dRs[j1][j5],
                                                                   tin->dRs[j2][j5],
                                                                   tin->dRs[j3][j5],
                                                                   tin->dRs[j4][j5]});
                                //transfer (4,2)
                                ans.transfer6[DR_Gen][DR5_Reco] += partial5 * weight6_42;
                            }
                        }
                    }
                } 

                for(const unsigned j2 : adj0){
                    T partial2 = partial1 * tin->ptrans[i0][j2];

                    unsigned DR2_Reco = max({DR1_Reco, tin->dRs[j0][j2],
                                                       tin->dRs[j1][j2]});

                    //transfer (2,1)
                    ans.transfer3[DR_Gen][DR2_Reco] += partial2 * weight3_21;

                    for(const unsigned j3 : adj0){
                        T partial3 = partial1 * tin->ptrans[i0][j3];

                        unsigned DR3_Reco = max({DR2_Reco, tin->dRs[j0][j3],
                                                           tin->dRs[j1][j3],
                                                           tin->dRs[j2][j3]});
                        //transfer (3,1)
                        ans.transfer3[DR_Gen][DR3_Reco] += partial3 * weight4_31;

                        for(const unsigned j4 : adj0){
                            T partial4 = partial1 * tin->ptrans[i0][j4];

                            unsigned DR4_Reco = max({DR3_Reco, tin->dRs[j0][j4],
                                                               tin->dRs[j1][j4],
                                                               tin->dRs[j2][j4],
                                                               tin->dRs[j3][j4]});
                            //transfer (4,1)
                            ans.transfer4[DR_Gen][DR4_Reco] += partial4 * weight5_41;

                            for(const unsigned j5 : adj0){
                                T partial5 = partial1 * tin->ptrans[i0][j5];

                                unsigned DR5_Reco = max({DR4_Reco, tin->dRs[j0][j5],
                                                                   tin->dRs[j1][j5],
                                                                   tin->dRs[j2][j5],
                                                                   tin->dRs[j3][j5],
                                                                   tin->dRs[j4][j5]});
                                //transfer (5,1)
                                ans.transfer5[DR_Gen][DR5_Reco] += partial5 * weight6_51;
                            }//end for j5
                        }//end for j4
                    }//end for j3
                }//end for j2
            }//end for j1
        }//end for j0
    }
};

#endif
