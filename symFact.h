#ifndef EECS_FAST_SYMFAC_H
#define EECS_FAST_SYMFAC_H

#include "usings.h"
#include "fastStructs.h"
#include "prev.h"

namespace fastEEC{
    template <typename>
    constexpr auto always_false = false;

    //TODO: simplify getSymFac with something constexpr??
    template <int order>
    constexpr unsigned symfacTableSize(){
        static_assert(always_false<order>, "symfacLookup not implemented for this order");
        return 0;
    }

    template <typename T, int order>
    constexpr auto symfacLookup(){
        static_assert(always_false<T>, "symfacLookup not implemented for this order");
        return std::array<T, 0>();
    }

    template <>
    constexpr auto symfacLookup<double, 1>(){
        return std::array<double, 1>({{1}});
    }

    template <>
    constexpr auto symfacLookup<double, 2>(){
        return std::array<double, 2>({{1, 2}});
    }

    template <>
    constexpr auto symfacLookup<double, 3>(){
        return std::array<double, 4>({{1, 3, 3, 6}});
    }

    template <>
    constexpr auto symfacLookup<double, 4>(){
        return std::array<double, 8>({{1, 4, 6, 12, 4, 12, 12, 24}});
    }

    template <>
    constexpr auto symfacLookup<double, 5>(){
        return std::array<double, 16>({{1, 5, 10, 20,
                                       10, 30, 30, 60,
                                       5, 20, 30, 60,
                                       20, 60, 60, 120}});
    }

    template <>
    constexpr auto symfacLookup<double, 6>(){
        return std::array<double, 32>({{1, 6, 15, 30,
                                       20, 60, 60, 120,
                                       15, 60, 90, 180,
                                       60, 180, 180, 360,
                                       6, 30, 60, 120,
                                       60, 180, 180, 360,
                                       30, 120, 180, 360,
                                       120, 360, 360, 720}});
    }

    template <typename T, int order>
    T getSymfac(const prev_t<T, order>& prev){
        static_assert(always_false<T>, "getSymfac not implemented for this order");
    }

    template <typename T>
    T getSymfac([[maybe_unused]] const prev_t<T, 2>& prev){
        return 1;
    }

    template <typename T>
    T getSymfac(const prev_t<T, 3>& prev){
        return (prev.is[0] == prev.is[1]) ? 1 : 2;;
    }

    template <typename T>
    T getSymfac(const prev_t<T, 4>& prev){
        if(prev.is[0] == prev.is[1]){
            return (prev.is[1] == prev.is[2]) ? 1 : 3;
        } else {
            return (prev.is[1] == prev.is[2]) ? 3 : 6;
        }
    }

    template <typename T>
    T getSymfac(const prev_t<T, 5>& prev){
        if(prev.is[0] == prev.is[1]){
            if(prev.is[1] == prev.is[2]){
                return (prev.is[2] == prev.is[3]) ? 1 : 4;
            } else {
                return (prev.is[2] == prev.is[3]) ? 6 : 12;
            }
        } else {
            if(prev.is[1] == prev.is[2]){
                return (prev.is[2] == prev.is[3]) ? 4 : 12;
            } else {
                return (prev.is[2] == prev.is[3]) ? 12 : 24;
            }
        }
    }

    template <typename T>
    T getSymfac(const prev_t<T, 6>& prev){
        if(prev.is[0] == prev.is[1]){
            if(prev.is[1] == prev.is[2]){
                if(prev.is[2] == prev.is[3]){
                    return (prev.is[3] == prev.is[4]) ? 1 : 5;
    } else {
                    return (prev.is[3] == prev.is[4]) ? 10 : 20;
                }
            } else {
                if(prev.is[2] == prev.is[3]){
                    return (prev.is[3] == prev.is[4]) ? 10 : 30;
                } else {
                    return (prev.is[3] == prev.is[4]) ? 30 : 60;
                }
            }
        } else {
            if(prev.is[1] == prev.is[2]){
                if(prev.is[2] == prev.is[3]){
                    return (prev.is[3] == prev.is[4]) ? 5 : 20;
                } else {
                    return (prev.is[3] == prev.is[4]) ? 30 : 60;
                }
            } else {
                if(prev.is[2] == prev.is[3]){
                    return (prev.is[3] == prev.is[4]) ? 20 : 60;
                } else {
                    return (prev.is[3] == prev.is[4]) ? 60 : 120;
                }
            }
        }
    }

    template <typename T>
    T getSymfac(const prev_t<T, 7>& prev){
        unsigned i0 = prev.is[0];
        unsigned i1 = prev.is[1];
        unsigned i2 = prev.is[2];
        unsigned i3 = prev.is[3];
        unsigned i4 = prev.is[4];
        unsigned i5 = prev.is[5];

    if (i0==i1){
            if(i1==i2){
                if(i2==i3){
                    if(i3==i4){
                        return (i4==i5) ? 1 : 6; //(6) vs (5, 1)
                    } else {
                        return (i4==i5) ? 15 : 30; //(4, 2) vs (4, 1, 1)
                    }
                } else {
                    if(i3==i4){
                        return (i4==i5) ? 20 : 60; //(3, 3) vs (3, 2, 1)
                    } else {
    return (i4==i5) ? 60 : 120; //(3, 1, 2) vs (3, 1, 1, 1)
                    }
                }
            } else {
                if(i2==i3){
                    if(i3==i4){
    return (i4==i5) ? 15 : 60; //(2, 4) vs (2, 3, 1)
                    } else {
                        return (i4==i5) ? 90 : 180; //(2, 2, 2) vs (2, 2, 1, 1)
                    }
                } else {
                    if(i3==i4){
                        return (i4==i5) ? 60 : 180; //(2, 1, 3) vs (2, 1, 2, 1)
                    } else {
                        return (i4==i5) ? 180 : 360; //(2, 1, 1, 2) vs (2, 1, 1, 1, 1)
                    }
                }
            }
        } else {
            if(i1==i2){
                if(i2==i3){
                    if(i3==i4){
                        return (i4==i5) ? 6 : 30; //(1, 5) vs (1, 4, 1)
                    } else {
    return (i4==i5) ? 60 : 120; //(1, 3, 2) vs (1, 3, 1, 1)
                    }
                } else {
                    if(i3==i4){
                        return (i4==i5) ? 60 : 180; //(1, 2, 3) vs (1, 2, 2, 1)
                    } else {
                        return (i4==i5) ? 180 : 360; //(1, 2, 1, 2) vs (1, 2, 1, 1, 1)
                    }
                }
            } else {
                if(i2==i3){
                    if(i3==i4){
                        return (i4==i5) ? 30 : 120; //(1, 1, 4) vs (1, 1, 3, 1)
                    } else {
                        return (i4==i5) ? 180 : 360; //(1, 1, 2, 2) vs (1, 1, 2, 1, 1)
                    }
                } else {
                    if(i3==i4){
                        return (i4==i5) ? 120 : 360; //(1, 1, 1, 3) vs (1, 1, 1, 2, 1)
                    } else {
                        return (i4==i5) ? 360 : 720; //(1, 1, 1, 1, 2) vs (1, 1, 1, 1, 1, 1)
                    }
                }
            }
        }
    }
}

#endif
