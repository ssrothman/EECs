#ifndef EECS_FASTTRANFER6_6_H
#define EECS_FASTTRANSFER6_6_H

namespace fastEEC{
    template <typename T>
    void transfer6_6(const unsigned i0,
                     const unsigned i1,
                     const unsigned i2,
                     const unsigned i3,
                     const unsigned i4,
                     const unsigned i5,

                     const unsigned DR_Gen,

                     const T weight6_111111,

                     const transferInputs<T>* tin,
                     result<T>& ans){
                     
        uvec adj0 = tin->adj.at(i0);
        uvec adj1 = tin->adj.at(i1);
        uvec adj2 = tin->adj.at(i2);
        uvec adj3 = tin->adj.at(i3);
        uvec adj4 = tin->adj.at(i4);
        uvec adj5 = tin->adj.at(i5);

        for(const unsigned j0 : adj0){
            for(const unsigned j1 : adj1){
                for(const unsigned j2 : adj2){
                    for(const unsigned j3 : adj3){
                        for(const unsigned j4 : adj4){
                            for(const unsigned j5 : adj5){
                                T partial = tin->ptrans[i0][j0] * 
                                            tin->ptrans[i1][j1] * 
                                            tin->ptrans[i2][j2] * 
                                            tin->ptrans[i3][j3] * 
                                            tin->ptrans[i4][j4] *
                                            tin->ptrans[i5][j5];
                                unsigned DR_Reco = max({tin->dRs[j0][j1], 
                                                        tin->dRs[j0][j2], 
                                                        tin->dRs[j0][j3], 
                                                        tin->dRs[j0][j4], 
                                                        tin->dRs[j0][j5],
                                                        tin->dRs[j1][j2], 
                                                        tin->dRs[j1][j3], 
                                                        tin->dRs[j1][j4], 
                                                        tin->dRs[j1][j5],
                                                        tin->dRs[j2][j3], 
                                                        tin->dRs[j2][j4], 
                                                        tin->dRs[j2][j5],
                                                        tin->dRs[j3][j4], 
                                                        tin->dRs[j3][j5],
                                                        tin->dRs[j4][j5]});
                                ans.transfer6[DR_Gen][DR_Reco] += partial * weight6_111111;
                            }
                        }
                    }
                }
            }
        }
    }
};

#endif
