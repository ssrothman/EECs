#ifndef EECS_FAST_PREV_H
#define EECS_FAST_PREV_H

#include "fastStructs.h"
#include "usings.h"

namespace fastEEC{
template <typename T, int order>
    struct prev_t{
        T partial; //product of all the particle energies so far
        T maxDR;   //floating point max DR so far
        unsigned maxDRbin; //bin of the max DR so far
        bool isPU = false; //are any of the particles so far PU?

        std::array<unsigned, order-1> is; //indices of the particles so far

        std::array<T, order-2> wts; //weights at each order so far
        std::array<unsigned, order-2> dRbins; //DR bins at each order so far

        /*for res3*/
        unsigned RL_res3_idx;
        unsigned xi_res3_idx;
        unsigned phi_res3_idx;

        /*for res4*/
        std::array<unsigned, 3> shape_res4_idx = {{0, 0, 0}};
        std::array<unsigned, 3> RL_res4_idx;
        std::array<unsigned, 3> r_res4_idx;
        std::array<unsigned, 3> ct_res4_idx;
        unsigned ires4;

        std::array<bool, 4> istri_res4 = {{false, false, false, false}};
        std::array<unsigned, 4> RL_res4tri_idx;
        std::array<unsigned, 4> r_res4tri_idx;
        std::array<unsigned, 4> ct_res4tri_idx;
    };

    template <typename T>
    struct prev_t<T, 1>{
        T partial;
        T maxDR;
        unsigned maxDRbin;
        bool isPU = false;

        std::array<unsigned, 0> is;
        std::array<T, 0> wts;

        unsigned RL_res3_idx;
        unsigned xi_res3_idx;
        unsigned phi_res3_idx;

        unsigned shape_res4_idx;
        unsigned RL_res4_idx;
        unsigned r_res4_idx;
    };
}

#endif
