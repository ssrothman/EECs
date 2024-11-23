#ifndef SROTHMAN_EECS_RES4MINR_H
#define SROTHMAN_EECS_RES4MINR_H

#include "util.h"
#include "fastStructs.h"
#include "SRothman/SimonTools/src/deltaR.h"

namespace fastEEC{
    template <typename T>
    bool runRes4_minR(const jetDetails_t<T>& jetDetails,
                      [[maybe_unused]] const res4shapesAxes_t& res4ax,

                      prev_t<T, 5>& next) noexcept {
        std::array<std::pair<char, T>, 6> dRs;
        dRs[0].first = 0;
        dRs[0].second = jetDetails.floatDRs[next.is[0]][next.is[1]];
        dRs[1].first = 1;
        dRs[1].second = jetDetails.floatDRs[next.is[0]][next.is[2]];
        dRs[2].first = 2;
        dRs[2].second = jetDetails.floatDRs[next.is[0]][next.is[3]];
        dRs[3].first = 3;
        dRs[3].second = jetDetails.floatDRs[next.is[1]][next.is[2]];
        dRs[4].first = 4;
        dRs[4].second = jetDetails.floatDRs[next.is[1]][next.is[3]];
        dRs[5].first = 5;
        dRs[5].second = jetDetails.floatDRs[next.is[2]][next.is[3]];

        std::sort(dRs.begin(), dRs.end(), [](const auto& a, const auto& b) {
            return a.second < b.second;
        });

        if (dRs[0].second == 0){
            return false;
        }


        /*
         * These are set such that 
         *  d(A,B) <= d(C,D) <= all other distances 
         *  and A != B != C != D
         *  if the two smallest distances share a point, we break and return
         */
        unsigned A=0, B=0, C=0, D=0;
        switch (dRs[0].first){
            case 0:
                if(dRs[1].first == 5){
                    A = next.is[0];
                    B = next.is[1];
                    C = next.is[2];
                    D = next.is[3];
                } else {
                    return false;
                }
                break;
            case 1:
                if (dRs[1].first == 4){
                    A = next.is[0];
                    B = next.is[2];
                    C = next.is[1];
                    D = next.is[3];
                } else {
                    return false;
                }
                break;
            case 2:
                if (dRs[1].first == 3){
                    A = next.is[0];
                    B = next.is[3];
                    C = next.is[1];
                    D = next.is[2];
                } else {
                    return false;
                }
                break;
            case 3:
                if (dRs[1].first == 2){
                    A = next.is[1];
                    B = next.is[2];
                    C = next.is[0];
                    D = next.is[3];
                } else {
                    return false;
                }
                break;
            case 4:
                if (dRs[1].first == 1){
                    A = next.is[1];
                    B = next.is[3];
                    C = next.is[0];
                    D = next.is[2];
                } else {
                    return false;
                }
                break;
            case 5:
                if (dRs[1].first == 0){
                    A = next.is[2];
                    B = next.is[3];
                    C = next.is[0];
                    D = next.is[1];
                } else {
                    return false;
                }
                break;
            default:
                assert(false);
                return false;
        } //end switch(dRs[0].first)


        T etaA = jetDetails.etas[A];
        T etaB = jetDetails.etas[B];
        T etaC = jetDetails.etas[C];
        T etaD = jetDetails.etas[D];

        T phiA = jetDetails.phis[A];
        T phiB = jetDetails.phis[B];
        T phiC = jetDetails.phis[C];
        T phiD = jetDetails.phis[D];

        T midAB_eta = (etaA + etaB) / 2;
        T midAB_phi = (phiA + phiB) / 2;

        T midCD_eta = (etaC + etaD) / 2;
        T midCD_phi = (phiC + phiD) / 2;

        T dR_mids = dR(midAB_eta, midAB_phi, midCD_eta, midCD_phi);
        T r1 = dRs[1].second;
        T r2 = dRs[0].second;

        T AB_eta = etaA - etaB;
        T AB_phi = deltaphi(phiA, phiB);

        T CD_eta = etaC - etaD;
        T CD_phi = deltaphi(phiC, phiD);

        T dot = AB_eta * CD_eta + AB_phi * CD_phi;
        T cosPhi = dot / (r1 * r2);
        T phi = acos(cosPhi);
        if (phi > M_PI / 2){
            phi = M_PI - phi;
        }

        //printf("%u, %u, %u, %u\n", A, B, C, D);
        //printf("dRs:\n");
        //for (const auto& dR : dRs) {
        //    printf("%d: %f\n", dR.first, dR.second);
        //}
        //printf("\n");

        //printf("dR_mids: %f\n", dR_mids);
        //printf("r1: %f\n", r1);
        //printf("r2: %f\n", r2);
        //printf("phi: %f\n", phi);
        //printf("\n");

        next.isRes4_minR = true;
        next.minR_R_idx = getIndex(dR_mids, res4ax.RL);
        next.minR_r1_idx = getIndex(r1/dR_mids, res4ax.r_minR);
        next.minR_r2_idx = getIndex(r2/dR_mids, res4ax.r_minR);
        next.minR_phi_idx = getIndex(phi, res4ax.phi_minR);

        return true;
    }//end runRes4_minR
}//end namespace fastEEC

#endif
