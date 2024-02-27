#ifndef EECS_CLEAR_H
#define EECS_CLEAR_H

namespace fastEEC{
    template <typename T, bool doPU, bool doTransfer, unsigned maxOrder>
    void clear(result<T>& ans, const unsigned NDR){
        ans.wts2.resize(NDR);
        fill(ans.wts2.begin(), ans.wts2.end(), 0);

        if constexpr(maxOrder >= 3){
            ans.wts3.resize(NDR);
            fill(ans.wts3.begin(), ans.wts3.end(), 0);
        }
        if constexpr(maxOrder >= 4){
            ans.wts4.resize(NDR);
            fill(ans.wts4.begin(), ans.wts4.end(), 0);
        }
        if constexpr(maxOrder >= 5){
            ans.wts5.resize(NDR);
            fill(ans.wts5.begin(), ans.wts5.end(), 0);
        }
        if constexpr(maxOrder >= 6){
            ans.wts6.resize(NDR);
            fill(ans.wts6.begin(), ans.wts6.end(), 0);
        }

        if constexpr (doPU){
            ans.wts2_PU.resize(NDR);
            fill(ans.wts2_PU.begin(), ans.wts2_PU.end(), 0);
            if constexpr(maxOrder >= 3){
                ans.wts3_PU.resize(NDR);
                fill(ans.wts3_PU.begin(), ans.wts3_PU.end(), 0);
            }
            if constexpr(maxOrder >= 4){
                ans.wts4_PU.resize(NDR);
                fill(ans.wts4_PU.begin(), ans.wts4_PU.end(), 0);
            }
            if constexpr(maxOrder >= 5){
                ans.wts5_PU.resize(NDR);
                fill(ans.wts5_PU.begin(), ans.wts5_PU.end(), 0);
            }
            if constexpr(maxOrder >= 6){
                ans.wts6_PU.resize(NDR);
                fill(ans.wts6_PU.begin(), ans.wts6_PU.end(), 0);
            }
        }

        if constexpr (doTransfer){
            ans.transfer2.resize(extents[NDR][NDR]);
            if constexpr(maxOrder >= 3){
                ans.transfer3.resize(extents[NDR][NDR]);
            }
            if constexpr(maxOrder >= 4){
                ans.transfer4.resize(extents[NDR][NDR]);
            }
            if constexpr(maxOrder >= 5){
                ans.transfer5.resize(extents[NDR][NDR]);
            }
            if constexpr(maxOrder >= 6){
                ans.transfer6.resize(extents[NDR][NDR]);
            }

            for(unsigned i=0; i<NDR; ++i){
                for(unsigned j=0; j<NDR; ++j){
                    ans.transfer2[i][j] = 0;
                    if constexpr(maxOrder >= 3){
                        ans.transfer3[i][j] = 0;
                    }
                    if constexpr(maxOrder >= 4){
                        ans.transfer4[i][j] = 0;
                    }
                    if constexpr(maxOrder >= 5){
                        ans.transfer5[i][j] = 0;
                    }
                    if constexpr(maxOrder >= 6){
                        ans.transfer6[i][j] = 0;
                    }
                }
            }
        }
    }
};

#endif
