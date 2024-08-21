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

        if (dR2Centroid < shapetol*shapetol){
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

        if (dR2TeeA < shapetol*shapetol){
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

        if (dR2TeeB < shapetol*shapetol){
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
        if (cosTheta > 1){
            theta = 0;
        } else {
            theta = std::acos(cosTheta);
}
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
        //printf("\tdot, %g\n", dot);
        //printf("\tcosTheta = %g\n", cosTheta);
        //printf("\tRL = %g\n", RL);
        //printf("\tRS = %g\n", RS);
        //printf("\tRS*RL = %g\n", RS*RL);
    }



    template <typename T>
    unsigned lookForShapes(const jetDetails_t<T>& jetDetails,
                           const res4shapesAxes_t& res4ax,
                           prev_t<T, 5>& next){

        T etaA, etaB, etaC, etaD;
        T phiA, phiB, phiC, phiD;
        T RAB, RCD;

        getEtasPhis(jetDetails, next, 
                etaA, etaB, etaC, etaD, 
                phiA, phiB, phiC, phiD, 
                RAB, RCD);

        if(RCD < 1e-8){
            //printf("\tRCD < 1e-8\n");
            next.shape_res4_idx[next.ires4] = 0;
            next.RL_res4_idx[next.ires4] = 0;
            next.r_res4_idx[next.ires4] = 0;
next.ct_res4_idx[next.ires4] = 0;

            return 0;
        }

        unsigned shape = linesCross(res4ax.shapetol,
                                    etaA, etaB, etaC, etaD,
                                    phiA, phiB, phiC, phiD);

        if(shape != 0){
            T r, theta;
            getrandangle(etaA, etaB, etaC, etaD,
                         phiA, phiB, phiC, phiD,
                         RAB, RCD,
                         r, theta);

            next.RL_res4_idx[next.ires4] = getIndex(RAB, res4ax.RL);

            if(shape == 1){
                //printf("\tDIPOLE:\n");
                next.shape_res4_idx[next.ires4] = 1;
                next.r_res4_idx[next.ires4] = getIndex(r, res4ax.r_dipole);
                next.ct_res4_idx[next.ires4] = getIndex(theta, res4ax.ct_dipole);
            } else {
                //printf("\tTEE:\n");
                next.shape_res4_idx[next.ires4] = 2;
                next.r_res4_idx[next.ires4] = getIndex(r, res4ax.r_tee);
                next.ct_res4_idx[next.ires4] = getIndex(theta, res4ax.ct_tee);
            }
            /*if ((RAB < 0.2) & (RAB > 0.1)){
                printf("%g, %g, %u, %u\n", 
                        theta, r,
                        next.ct_res4_idx[next.ires4], 
                        next.r_res4_idx[next.ires4]);
            }*/
            //printf("\tr = %g\n", r);
            //printf("\trIdx = %u\n", next.r_res4_idx);
            //printf("\ttheta = %g\n", theta);
            //printf("\tthetaIdx = %u\n", next.ct_res4_idx);
            return shape;
        } else {
            //printf("\tnoshape\n");
            next.shape_res4_idx[next.ires4] = 0;
            return 0;
        }
    }

    template <typename T>
    void getEtasPhis(const jetDetails_t<T>& jetDetails,
                     const prev_t<T, 5>& prev,

                     T& etaA, T& etaB, T& etaC, T& etaD,
                     T& phiA, T& phiB, T& phiC, T& phiD,

                     T& RAB, T& RCD){
        static constexpr std::array<unsigned, 3> As = {{0, 0, 0}};
        static constexpr std::array<unsigned, 3> Bs = {{1, 2, 3}};
        static constexpr std::array<unsigned, 3> Cs = {{2, 3, 1}};
        static constexpr std::array<unsigned, 3> Ds = {{3, 1, 2}};

        unsigned A = prev.is[As[prev.ires4]];
        unsigned B = prev.is[Bs[prev.ires4]];
        unsigned C = prev.is[Cs[prev.ires4]];
        unsigned D = prev.is[Ds[prev.ires4]];

        RAB = jetDetails.floatDRs[A][B];
        RCD = jetDetails.floatDRs[C][D];

        if (RCD > RAB){
            std::swap(RAB, RCD);
            std::swap(B, C);
            std::swap(A, D);
        }

        etaA = jetDetails.etas[A];
        etaB = jetDetails.etas[B];
        etaC = jetDetails.etas[C];
        etaD = jetDetails.etas[D];

        phiA = jetDetails.phis[A];
        phiB = jetDetails.phis[B];
        phiC = jetDetails.phis[C];
        phiD = jetDetails.phis[D];

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
    }


    template <typename T, int genOrder>
    void runRes4_nocheck(const jetDetails_t<T>& jetDetails,
                         const res4shapesAxes_t& res4ax,

                         const prev_t<T, genOrder>& prevGen,

                         prev_t<T, 5>& next){

        

        for(next.ires4=0; next.ires4<3; ++next.ires4){
            unsigned shapeGen = prevGen.shape_res4_idx[next.ires4];
            if(prevGen.shape_res4_idx[next.ires4] == 0){
                continue;
            }

            T etaA, etaB, etaC, etaD;
            T phiA, phiB, phiC, phiD;
            T RAB, RCD;

            getEtasPhis(jetDetails, next, 
                    etaA, etaB, etaC, etaD, 
                    phiA, phiB, phiC, phiD, 
                    RAB, RCD);

            T r, theta;
            getrandangle(etaA, etaB, etaC, etaD,
                         phiA, phiB, phiC, phiD,
                         RAB, RCD,
                         r, theta);

            next.RL_res4_idx[next.ires4] = getIndex(RAB, res4ax.RL);

            if (shapeGen == 1){
                next.shape_res4_idx[next.ires4] = 1;
                next.r_res4_idx[next.ires4] = getIndex(r, res4ax.r_dipole);
                next.ct_res4_idx[next.ires4] = getIndex(theta, res4ax.ct_dipole);
            } else {
                next.shape_res4_idx[next.ires4] = 2;
                next.r_res4_idx[next.ires4] = getIndex(r, res4ax.r_tee);
                next.ct_res4_idx[next.ires4] = getIndex(theta, res4ax.ct_tee);
            }
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

        auto minel = min_element(dRs.begin(), dRs.end());
        if (*minel < 1e-16){
            return;
        }

        int nfound = 0;
        for(next.ires4=0; next.ires4<3; ++next.ires4){
            //printf("Trying %u, %u, %u, %u\n", As[i], Bs[i], Cs[i], Ds[i]);
            //printf("%u %u %u %u:\n", As[i], Bs[i], Cs[i], Ds[i]);
            if (lookForShapes(jetDetails, res4ax, next)){
                ++nfound;
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
                //return;
            }
        }
        if(nfound > 2){
            /*printf("Found more than two shapes\n");
            printf("\t");
            printVec(next.is);
            printf("\n");
            printf("\t(%g, %g)\n", jetDetails.etas[next.is[0]], jetDetails.phis[next.is[0]]);
            printf("\t(%g, %g)\n", jetDetails.etas[next.is[1]], jetDetails.phis[next.is[1]]);
            printf("\t(%g, %g)\n", jetDetails.etas[next.is[2]], jetDetails.phis[next.is[2]]);
            printf("\t(%g, %g)\n", jetDetails.etas[next.is[3]], jetDetails.phis[next.is[3]]);
            printf("\tRL = %g\n", RL);*/
            /*printf("\tABCD: %u\n", lookForShapes(next.is[0],
                                                 next.is[1],
                                                 next.is[2],
                                                 next.is[3],
                                                 jetDetails,
                                                 res4ax,
                                                 next));
            printf("\tACBD: %u\n", lookForShapes(next.is[0],
                                                 next.is[2],
                                                 next.is[1],
                                                 next.is[3],
                                                 jetDetails,
                                                 res4ax,
                                                 next));
            printf("\tADBC: %u\n", lookForShapes(next.is[0],
                                                 next.is[3],
                                                 next.is[1],
                                                 next.is[2],
                                                 jetDetails,
                                                 res4ax,
                                                 next));*/
        }
    } 
}

#endif
