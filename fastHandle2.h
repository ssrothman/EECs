#ifndef EECS_FASTHANDLE2_H
#define EECS_FASTHANDLE2_H

namespace fastEEC{
    template <unsigned maxOrder, typename T, bool doPU, bool doTransfer>
    inline void handle2(const unsigned i0, 
                 const unsigned i1,
                 const T E0,
                 const T E1,
                 const T sqE0,
                 const T sqE1,
                 const T partial0,
                 const unsigned DR0,
                 const umat& dRs, 

                 T& partial1,
                 result<T>& ans,
                 T& partial1, bool& isPU, unsigned& DR1,
                 const vector<bool> *const PU = nullptr,
                 const transferInputs<T>* const tin = nullptr){
        static_assert(maxOrder >=2 && maxOrder <=6);

        DR1 = max({DR0, dRs(i0,i1)});
        partial1 = partial0 * E1;
        T sqpartial1 = square(partial1);
        if constexpr(doPU){
            isPU = isPU || PU->at(i1);
        }

        T weight2, weight3, weight4, weight5, weight6;

        T weight3_21, weight3_12;
        T weight4_31, weight4_22, weight4_13;
        T weight5_41, weight5_32, weight5_23, weight5_14;
        T weight6_51, weight6_42, weight6_33, weight6_24, weight6_15;

        if constexpr (maxOrder >= 2){
            //(1, 1)
            weight2 = 2 * partial1;
            ans.wts2[DR1] += weight2;
            if constexpr(doPU){
                if (isPU){
                    ans.wts2_PU[DR1] += weight2;
                }
            }
        }

        if constexpr (maxOrder >= 3){
            T weight3_prefac = 3*partial1;
            //(2, 1)
            T weight3_21 = weight3_prefac * E0;
            //(1, 2)
            T weight3_12 = weight3_prefac * E1;
            weight3 = weight3_12 + weight3_21;
            ans.wts3[DR1] += weight3;
            if constexpr (doPU){
                if (isPU){
                    ans.wts3_PU[DR1] += weight3;
                }
            }
        }

        if constexpr (maxOrder >=4 ){
            T weight4_prefac = 4*partial1;
            //(3, 1)
            T weight4_31 = weight4_prefac * sqE0;
            //(1, 3)
            T weight4_13 = weight4_prefac * sqE1;
            //(2, 2)
            T weight4_22 = 6*sqpartial1;
            weight4 = weight4_13 + weight4_31 + weight4_22;
            ans.wts4[DR1] += weight4;
            if constexpr (doPU){
                if (isPU){
                    ans.wts4_PU[DR1] += weight4;
                }
            }
        }

        if constexpr (maxOrder >= 5){
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
            if constexpr (doPU){
                if (isPU){
                    ans.wts5_PU[DR1] += weight5;
                }
            }
        }

        if constexpr (maxOrder >= 6){
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
    }
}

#endif
