#ifndef EECS_FAST_PREV_H
#define EECS_FAST_PREV_H

#include "fastStructs.h"
#include "usings.h"

namespace fastEEC{
template <typename T, int order>
    struct prev_t{
        T partial;
        T maxDR;
        unsigned maxDRbin;
        bool isPU = false;

        std::array<unsigned, order-1> is;

        /*for res3*/
        unsigned RL_res3_idx;
        unsigned xi_res3_idx;
        unsigned phi_res3_idx;

        /*for res4*/
        unsigned shape_res4_idx;
        unsigned RL_res4_idx;
        unsigned r_res4_idx;
        unsigned ct_res4_idx;
    };
}

#endif
