#ifndef EECS_FASTTRANSFER6_1_H
#define EECS_FASTTRANSFER6_1_H

#include "fastStructs.h"

namespace fastEEC{
    template <typename T>
    void transfer6_1(const unsigned i0,
                     const unsigned DR_Gen,
                     const T weight2, 
                     const T weight3, 
                     const T weight4, 
                     const T weight5, 
                     const T weight6,
                     const transferInputs<T>* tin,
                     result<T>& ans){
        uvec adj0 = tin->adj.at(i0);
        for (const auto& j0 : adj0){
            T partialtrans0 = tin->ptrans[i0][j0];
            for (const auto& j1 : adj0){
                T partialtrans1 = tin->ptrans[i0][j1] * partialtrans0;
                unsigned DR1_Reco = tin->dRs[j0][j1];
                ans.transfer2[DR_Gen][DR1_Reco] += partialtrans1 * weight2;
                for(const auto& j2 : adj0){
                    T partialtrans2 = tin->ptrans[i0][j2] * partialtrans1;
                    unsigned DR2_Reco = max({DR1_Reco, tin->dRs[j0][j2],
                                                       tin->dRs[j1][j2]});
                    ans.transfer3[DR_Gen][DR2_Reco] += partialtrans2 * weight3;
                    for(const auto& j3 : adj0){
                        T partialtrans3 = tin->ptrans[i0][j3] * partialtrans2;
                        unsigned DR3_Reco = max({DR2_Reco, tin->dRs[j0][j3],
                                                           tin->dRs[j1][j3],
                                                           tin->dRs[j2][j3]});
                        ans.transfer4[DR_Gen][DR3_Reco] += partialtrans3 * weight4;
                        for(const auto& j4 : adj0){
                            T partialtrans4 = tin->ptrans[i0][j4] * partialtrans3;
                            unsigned DR4_Reco = max({DR3_Reco, tin->dRs[j0][j4],
                                                               tin->dRs[j1][j4],
                                                               tin->dRs[j2][j4],
                                                               tin->dRs[j3][j4]});
                            ans.transfer5[DR_Gen][DR4_Reco] += partialtrans4 * weight5;
                            for(const auto& j5 : adj0){
                                T partialtrans5 = tin->ptrans[i0][j5] * partialtrans4;
                                unsigned DR5_Reco = max({DR4_Reco, tin->dRs[j0][j5],
                                                                   tin->dRs[j1][j5],
                                                                   tin->dRs[j2][j5],
                                                                   tin->dRs[j3][j5],
                                                                   tin->dRs[j4][j5]});
                                ans.transfer6[DR_Gen][DR5_Reco] += partialtrans5 * weight6;
                            }//end for j5
                        }//end for j4
                    }//end for j3
                }//end for j2
            }//end for j1
        }//end for j0
    }
};

#endif
