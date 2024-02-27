#ifndef EECS_FASTHANDLE1_H
#define EECS_FASTHANDLE1_H

namespace fastEEC{
    template <unsigned maxOrder, typename T, bool doPU, bool doTransfer>
    inline void handle1(const unsigned i0,
                 const T E0,
                 result<T>& ans,
                 T& partial0, bool& isPU, unsigned& DR0,
                 const vector<bool> *const PU = nullptr,
                 const transferInputs<T>* const tin = nullptr){
        static_assert(maxOrder >=2 && maxOrder <=6);

        DR0 = 0;
        partial0 = E0;
        if constexpr(doPU){
            isPU = PU->at(i0);
        }

        T weight2, weight3, weight4, weight5, weight6;
        
        if constexpr (maxOrder >= 2){
            //partition (2)
            weight2 = square(partial0);
            ans.wts2[DR0] += weight2;
        }

        if constexpr (maxOrder >= 3){
            //partition (3)
            weight3 = partial0 * weight2;
            ans.wts3[DR0] += weight3;
            if constexpr(doPU){
                if (isPU){
                    ans.wts3_PU[DR0] += weight3;
                }
            }
        }

        if constexpr (maxOrder >= 4){
            //partition (4)
            weight4 = square(weight2);
            ans.wts4[DR0] += weight4;
            if constexpr(doPU){
                if (isPU){
                    ans.wts4_PU[DR0] += weight4;
                }
            }
        }

        if constexpr (maxOrder >= 5){
            //partition (5)
            weight5 = partial0 * weight4;
            ans.wts5[DR0] += weight5;
            if constexpr(doPU){
                if (isPU){
                    ans.wts5_PU[DR0] += weight5;
                }
            }
        }

        if constexpr (maxOrder >= 6){
            //partition (6)
            weight6 = square(weight3);
            ans.wts6[DR0] += weight6;
            if constexpr(doPU){
                if (isPU){
                    ans.wts6_PU[DR0] += weight6;
                }
            }
        }
        
        if constexpr (doTransfer){
            transfer6_1(i0, DR0,
                        weight2, weight3, weight4, weight5, weight6,
                        tin, ans);
        }
    }
};

#endif
