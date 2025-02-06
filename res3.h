#ifndef EECS_FAST_RES3_H
#define EECS_FAST_RES3_H

#include "util.h"
#include "fastStructs.h"
#include "SRothman/SimonTools/src/util.h"

namespace fastEEC{

    template <typename T>
    void runRes3(const jetDetails_t<T>& jetDetails,
              const res3axes_t& res3ax,

              prev_t<T, 4>& next){
        T RL = jetDetails.floatDRs[next.is[0]][next.is[1]];
        T RM = jetDetails.floatDRs[next.is[0]][next.is[2]];
        T RS = jetDetails.floatDRs[next.is[1]][next.is[2]];

        if (RS > RM){
            std::swap(RS, RM);
        }
        if (RM > RL){
            std::swap(RM, RL);
        }
        if (RS > RM){
            std::swap(RS, RM);
        }
        //they're now sorted

        T xi, phi;

        if (RM == 0){
            xi = 0;
            phi = 0;
        } else {
            xi = RS/RM;
            if (RS == 0){
                phi = 0;
            } else {
                phi = std::abs(std::asin(std::sqrt(1 - square((RL-RM)/RS))));
            }
        }
#endif

        next.RL_res3_idx = getIndex(RL, res3ax.RL);
        next.xi_res3_idx = getIndex(xi, res3ax.xi);
        next.phi_res3_idx = getIndex(phi, res3ax.phi);
    }
}

#endif

