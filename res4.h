#ifndef EECS_FAST_RES4_H
#define EECS_FAST_RES4_H

#include "util.h"
#include "fastStructs.h"

namespace fastEEC{
    template <typename T>
    unsigned linesCross(const float shapetol,
                    const T etaA,
                    const T etaB,
                    const T etaC,
                    const T etaD,
                    const T phiA,
                    const T phiB,
                    const T phiC,
                    const T phiD){
        /* 
         * Assumes that lines are already sorted
         * phiA < phiB
         * phiC < phiD
         */

        //printf("Line 1: (%g, %g) -> (%g, %g)\n", etaA, phiA, etaB, phiB);
        //printf("Line 2: (%g, %g) -> (%g, %g)\n", etaC, phiC, etaD, phiD);

        //centroid of CD
        T etaCD = 0.5 * (etaC + etaD);
        T phiCD = 0.5 * (phiC + phiD);

        //centroid of AB
        T etaAB = 0.5 * (etaA + etaB);
        T phiAB = 0.5 * (phiA + phiB);

        /*
         * Check for dipole configuration
         */
        T detaCentroid = etaCD - etaAB;
        T dphiCentroid = phiCD - phiAB;
        if(dphiCentroid > M_PI){
            dphiCentroid -= 2*M_PI;
        }
        if(dphiCentroid < -M_PI){
            dphiCentroid += 2*M_PI;
        }
        const T dR2Centroid = detaCentroid*detaCentroid + dphiCentroid*dphiCentroid;

        if (dR2Centroid < shapetol){
            return 1;
        }

        /*
         * Check for tee configuration
         */
        T detaTeeA = etaA - etaCD;
        T dphiTeeA = phiA - phiCD;
        if(dphiTeeA > M_PI){
            dphiTeeA -= 2*M_PI;
        }
        if(dphiTeeA < -M_PI){
            dphiTeeA += 2*M_PI;
        }

        const T dR2TeeA = detaTeeA*detaTeeA + dphiTeeA*dphiTeeA;

        if (dR2TeeA < shapetol){
            return 2;
        }

        T detaTeeB = etaB - etaCD;
        T dphiTeeB = phiB - phiCD;
        if(dphiTeeB > M_PI){
            dphiTeeB -= 2*M_PI;
        }
        if(dphiTeeB < -M_PI){
            dphiTeeB += 2*M_PI;
        }

        const T dR2TeeB = detaTeeB*detaTeeB + dphiTeeB*dphiTeeB;

        if (dR2TeeB < shapetol){
            return 2;
        }

        return 0;
    }

    template <typename T>
    void getrandangle(const T etaA, const T etaB, const T etaC, const T etaD,
                      const T phiA, const T phiB, const T phiC, const T phiD,
                      const T RL, const T RS,
                      T& r, T& theta){
        //we've already taken pains to ensure that the phi distances 
        //are less than pi
        //so we don't need to do the proper deltaphi calculation

        T etaAB = etaB - etaA;
        T phiAB = phiB - phiA;

        T etaCD = etaD - etaC;
        T phiCD = phiD - phiC;

        T dot = etaAB*etaCD + phiAB*phiCD;

        T cosTheta = dot/(RL*RS);
        
        r = RS/RL;
        theta = std::acos(cosTheta);
        /*if(std::isnan(theta)){
            printf("------------- NAN THETA!!!!! ---------\n");
            printf("\n");
            printf("(%g, %g) -> (%g, %g)\n", etaA, phiA, etaB, phiB);
            printf("(%g, %g) -> (%g, %g)\n", etaC, phiC, etaD, phiD);
            printf("RL = %g\n", RL);
            printf("RS = %g\n", RS);
            printf("etaAB = %g\n", etaAB);
            printf("phiAB = %g\n", phiAB);
            printf("etaCD = %g\n", etaCD);
            printf("phiCD = %g\n", phiCD);
            printf("dot = %g\n", dot);
            printf("cosTheta = %g\n", cosTheta);
            printf("\n");
            printf("--------------------------------------\n");
        }*/
        if (theta > M_PI/2){
            theta = M_PI - theta;
        }

        //printf("r = %g\n", r);
        //printf("theta = %g\n", theta);
    }

    template <typename T>
    bool lookForShapes(unsigned A,
                       unsigned B,
                       unsigned C,
                       unsigned D,

                       const jetDetails_t<T>& jetDetails,
                       const res4shapesAxes_t& res4ax,

                       prev_t<T, 5>& next){

        T RAB = jetDetails.floatDRs[A][B];
        T RCD = jetDetails.floatDRs[C][D];

        if (RCD > RAB){
            std::swap(RAB, RCD);
            std::swap(B, C);
            std::swap(A, D);
        }

        //printf("RAB = %g\n", RAB);
        //printf("RCD = %g\n", RCD);

        T etaA = jetDetails.etas[A];
        T etaB = jetDetails.etas[B];
        T etaC = jetDetails.etas[C];
        T etaD = jetDetails.etas[D];

        T phiA = jetDetails.phis[A];
        T phiB = jetDetails.phis[B];
        T phiC = jetDetails.phis[C];
        T phiD = jetDetails.phis[D];


        if (phiA > phiB){
            std::swap(phiA, phiB);
            std::swap(etaA, etaB);
        }
        if (phiC > phiD){
            std::swap(phiC, phiD);
            std::swap(etaC, etaD);
        }

        if (phiB - phiA > M_PI){
            phiB -= 2*M_PI;
            std::swap(phiA, phiB);
            std::swap(etaA, etaB);
        }

        if (phiD - phiC > M_PI){
            phiD -= 2*M_PI;
            std::swap(phiC, phiD);
            std::swap(etaC, etaD);
        }
        
        //now we know for certain that phiA < phiB and phiC < phiD
        //and that phiB - phiA <= M_PI and phiD - phiC <= M_PI

        //in principle we need to worry now about whether C-D would cross
        //if we were in the right +-2pi window
        //
        //the only solution I can think of is just to explicitly check
        //all three options
        //
        //Luckily I can't think of any way that it would be off by 
        //more than 2pi
        next.shape_res4_idx = linesCross(res4ax.shapetol,
                                        etaA, etaB, etaC, etaD,
                                        phiA, phiB, phiC, phiD);

        if(next.shape_res4_idx != 0){
            T r, theta;
            getrandangle(etaA, etaB, etaC, etaD,
                         phiA, phiB, phiC, phiD,
                         RAB, RCD,
                         r, theta);
            if(next.shape_res4_idx == 1){
                //printf("DIPOLE:\n");
                next.r_res4_idx = getIndex(r, res4ax.r_dipole);
                next.ct_res4_idx = getIndex(theta, res4ax.ct_dipole);
            } else {
                //printf("TEE:\n");
                next.r_res4_idx = getIndex(r, res4ax.r_tee);
                next.ct_res4_idx = getIndex(theta, res4ax.ct_tee);
            }
            //printf("\tr = %g\n", r);
            //printf("\trIdx = %u\n", next.r_res4_idx);
            //printf("\ttheta = %g\n", theta);
            //printf("\tthetaIdx = %u\n", next.ct_res4_idx);
            return true;
        } else {
            next.r_res4_idx = 0;
            next.ct_res4_idx = 0;
            return false;
        }
    }

    template <typename T>
    void runRes4(const jetDetails_t<T>& jetDetails,
                 const res4shapesAxes_t& res4ax,

                 prev_t<T, 5>& next){
        std::array<T, 6> dRs = {{
            jetDetails.floatDRs[next.is[0]][next.is[1]],
            jetDetails.floatDRs[next.is[0]][next.is[2]],
            jetDetails.floatDRs[next.is[0]][next.is[3]],
            jetDetails.floatDRs[next.is[1]][next.is[2]],
            jetDetails.floatDRs[next.is[1]][next.is[3]],
            jetDetails.floatDRs[next.is[2]][next.is[3]]
        }};
        auto maxel = max_element(dRs.begin(), dRs.end());
        T RL = *maxel;

        next.RL_res4_idx = getIndex(RL, res4ax.RL);

        auto minel = min_element(dRs.begin(), dRs.end());
        if (*minel < 1e-8){
            next.shape_res4_idx = 0;
            next.r_res4_idx = 0;
            next.ct_res4_idx = 0;
            return;
        }

        constexpr std::array<unsigned, 3> As = {{0, 0, 0}};
        constexpr std::array<unsigned, 3> Bs = {{1, 2, 3}};
        constexpr std::array<unsigned, 3> Cs = {{2, 3, 1}};
        constexpr std::array<unsigned, 3> Ds = {{3, 1, 2}};

        for (unsigned i=0; i<3; ++i){
            //printf("Trying %u, %u, %u, %u\n", As[i], Bs[i], Cs[i], Ds[i]);
            if (lookForShapes(next.is[As[i]], next.is[Bs[i]], next.is[Cs[i]], next.is[Ds[i]],
                              jetDetails, res4ax, next)){
                //printf("Found shape %u\n", next.shape_res4_idx);
                //printf("\tline 1: (%0.2f, %0.2f) -> (%0.2f, %0.2f)\n", 
                //        jetDetails.etas[next.is[As[i]]], jetDetails.phis[next.is[As[i]]],
                //        jetDetails.etas[next.is[Bs[i]]], jetDetails.phis[next.is[Bs[i]]]);
                //printf("\tline 2: (%0.2f, %0.2f) -> (%0.2f, %0.2f)\n", 
                //        jetDetails.etas[next.is[Cs[i]]], jetDetails.phis[next.is[Cs[i]]],
                //        jetDetails.etas[next.is[Ds[i]]], jetDetails.phis[next.is[Ds[i]]]);
                //printf("\t");
                //printVec(next.is);
                //printf("\n");
                return;
            }
        }
    } 
}

#endif
