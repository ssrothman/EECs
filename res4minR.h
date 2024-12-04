#ifndef SROTHMAN_EECS_RES4MINR_H
#define SROTHMAN_EECS_RES4MINR_H

#include "util.h"
#include "fastStructs.h"
#include "SRothman/SimonTools/src/deltaR.h"

namespace fastEEC{
    /*
     * Run the experimntal/new "minR" configuration for the res4
     *
     * This looks for configurations where:
     *      The two smallest distnces are between AB and CD
     *      AB and CD are disjoint
     *
     * The interesting quantities are
     *      R = distance between dR(midpoint AB, midpoint CD)
     *      r1 = max{dR(A,B)/R, dR(C,D)/R}
     *      r2 = min{dR(A,B)/R, dR(C,D)/R}
     *      theta = angle betwen AB and CD
     *      phi1 = angle between r1 and R
     *      phi2 = angle between r2 and R
     *
     * We can't possibly bin in all of these variables at once
     * I need to think a bit about what to do here...
     *
     * It's also slightly overconstrained, because phi1, phi2 are related
     *
     * I think it depends on what we are looking for.
     *
     * To look for entanglement:
     *      I think the most interesting quantity here is either theta, 
     *      or phi1-phi2
     * To look for spin interference:
     *      I think the most interesting quantity here is phi1 and phi2
     *
     * In principle the effect size will be a function of r1/R, r2/R
     * It's not completely obvious to me the size of that effect
     */
    template <typename T>
    bool runRes4_minR(const jetDetails_t<T>& jetDetails,
                      const res4shapesAxes_t& res4ax,

                      prev_t<T, 5>& next) noexcept {
        std::array<std::pair<char, T>, 6> dRs;
        dRs[0].second = jetDetails.floatDRs[next.is[0]][next.is[1]];
        dRs[1].second = jetDetails.floatDRs[next.is[0]][next.is[2]];
        dRs[2].second = jetDetails.floatDRs[next.is[0]][next.is[3]];
        dRs[3].second = jetDetails.floatDRs[next.is[1]][next.is[2]];
        dRs[4].second = jetDetails.floatDRs[next.is[1]][next.is[3]];
        dRs[5].second = jetDetails.floatDRs[next.is[2]][next.is[3]];

        dRs[0].first = 0;
        dRs[1].first = 1;
        dRs[2].first = 2;
        dRs[3].first = 3;
        dRs[4].first = 4;
        dRs[5].first = 5;

        //sort in increasing order
        std::sort(dRs.begin(), dRs.end(), [](const auto& a, const auto& b) {
            return a.second < b.second;
        });

        //if the smallest distance is 0, we can't do anything
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

        //distances
        T dR_mids = dR(midAB_eta, midAB_phi, midCD_eta, midCD_phi);
        T r1 = dRs[1].second;
        T r2 = dRs[0].second;

        //angle between AB and CD
        T AB_eta = etaB - etaA;
        T AB_phi = deltaPhi(phiA, phiB);

        T CD_eta = etaD - etaC;
        T CD_phi = deltaPhi(phiC, phiD);

        T dot = AB_eta * CD_eta + AB_phi * CD_phi;
        T cosTheta = dot / (r1 * r2);
        T theta = acos(cosTheta);
        if (theta > M_PI / 2){
            theta = M_PI - theta;
        }

        //angle between r1 and R
        T mid_eta = midCD_eta - midAB_eta;
        T mid_phi = deltaPhi(midAB_phi, midCD_phi);

        dot = CD_eta * mid_eta + CD_phi * mid_phi;
        T cosPhi1 = dot / (r1 * dR_mids);
        T phi1 = acos(cosPhi1);
        if (phi1 > M_PI / 2){
            phi1 = M_PI - phi1;
        }

        //angle between r2 and R
        dot = AB_eta * mid_eta + AB_phi * mid_phi;
        T cosPhi2 = dot / (r2 * dR_mids);
        T phi2 = acos(cosPhi2);
        if (phi2 > M_PI / 2){
            phi2 = M_PI - phi2;
        }

        double phidiff = std::abs(phi1 - phi2);

        next.isRes4_minR = true;
        next.minR_R_idx = getIndex(dR_mids, res4ax.RL);
        next.minR_rmax_idx = getIndex(r1/dR_mids, res4ax.r_minR);
        next.minR_phi1_idx = getIndex(phi1, res4ax.phi_minR);
        next.minR_phi2_idx = getIndex(phi2, res4ax.phi_minR);
        next.minR_phidiff_idx = getIndex(phidiff, res4ax.phi_minR);
        next.minR_theta_idx = getIndex(theta, res4ax.theta_minR);

        return true;
    }//end runRes4_minR
}//end namespace fastEEC

#endif
