#ifndef EECS_RESOLVED4_H
#define EECS_RESOLVED4_H

constexpr double invsqrt2 = 0.7071067811865475;

namespace fastEEC{
    template <typename T>
    void fixedshape4(const unsigned i0,
                     const unsigned i1,
                     const unsigned i2,
                     const unsigned i3,
                     const resolvedInputs<T>& rin,
                     unsigned& shape_idx){
        printf("fixedshape4 %u %u %u %u\n", i0, i1, i2, i3);
        fflush(stdout);
        T R1 = rin.floatDRs[i0][i1];
        T R2 = rin.floatDRs[i0][i2];
        T R3 = rin.floatDRs[i0][i3];
        T R4 = rin.floatDRs[i1][i2];
        T R5 = rin.floatDRs[i1][i3];
        T R6 = rin.floatDRs[i2][i3];
        //printf("\tfloatDRs size is (%lu,%lu)\n", rin.floatDRs.shape()[0],
        //                                     rin.floatDRs.shape()[1]);
        //printf("\tfloatDRs ndim = %lu\n", rin.floatDRs.num_dimensions());
        printf("\tgot floatDRs\n");
        printf("\t\tR1 = [%u,%u]=%f\n", i0, i1, R1);
        printf("\t\tR2 = [%u,%u]=%f\n", i0, i2, R2);
        printf("\t\tR3 = [%u,%u]=%f\n", i0, i3, R3);
        printf("\t\tR4 = [%u,%u]=%f\n", i1, i2, R4);
        printf("\t\tR5 = [%u,%u]=%f\n", i1, i3, R5);
        printf("\t\tR6 = [%u,%u]=%f\n", i2, i3, R6);
        fflush(stdout);

        std::array<T, 6> subdRs({{R1, R2, R3, R4, R5, R6}});
        std::sort(subdRs.begin(), subdRs.end(), std::greater<T>());
        printf("\tsorted\n");
        printf("\t\t%g\n", subdRs[0]);
        printf("\t\t%g\n", subdRs[1]);
        printf("\t\t%g\n", subdRs[2]);
        printf("\t\t%g\n", subdRs[3]);
        printf("\t\t%g\n", subdRs[4]);
        printf("\t\t%g\n", subdRs[5]);
        fflush(stdout);

        T RL = subdRs[0];
        printf("\tRL = %g\n", RL);
        fflush(stdout);

        if ((RL - subdRs[1]) < rin.shapetol * RL &&
            (RL - subdRs[2]) < rin.shapetol * RL){

            if ((RL - subdRs[3]) < rin.shapetol * RL){
                T Rsquareedge = RL * invsqrt2;
                printf("\tRsquareedge = %g\n", Rsquareedge);
                fflush(stdout);
                if (std::abs(Rsquareedge - subdRs[4]) < rin.shapetol * RL &&
                    std::abs(Rsquareedge - subdRs[5]) < rin.shapetol * RL){

                    shape_idx = 1; //square
                    printf("\tshape_idx = 1\n");
                    fflush(stdout);
                    return;
                }
            } else {
                T Rinnertri = 0.6666666666666666 * RL;
                printf("\tRinnertri = %g\n", Rinnertri);
                fflush(stdout);

                if (std::abs(Rinnertri - subdRs[3]) < rin.shapetol * RL &&
                    std::abs(Rinnertri - subdRs[4]) < rin.shapetol * RL &&
                    std::abs(Rinnertri - subdRs[5]) < rin.shapetol * RL){

                    shape_idx = 2; //triangle 
                    printf("\tshape_idx = 2\n");
                    fflush(stdout);
                    return;
                }
            }
        }
        shape_idx = 0;
        printf("\tshape_idx = 0\n");
        fflush(stdout);
        return;
    }

    template <typename T>
    bool linesCross(const T etaA,
                    const T phiA,
                    const T etaB,
                    const T phiB,
                    const T etaC,
                    const T phiC,
                    const T etaD,
                    const T phiD,
                    const T shapetol){
        //bounding box check in phi
        if(phiD < phiA || phiB < phiC){
            return false;
        }

        //bounding box check in eta
        T sortedEtaAB[2] = {etaA, etaB};
        T sortedEtaCD[2] = {etaC, etaD};
        std::sort(sortedEtaAB, sortedEtaAB+2);
        std::sort(sortedEtaCD, sortedEtaCD+2);

        if(sortedEtaAB[1] < sortedEtaCD[0] || sortedEtaCD[1] < sortedEtaAB[0]){
            return false;
        }

        T ABeta = etaB - etaA;
        T ABphi = phiB - phiA;

        T CDeta = etaD - etaC;
        T CDphi = phiD - phiC;

        T slopecross = ABeta*CDphi - ABphi*CDeta;
        if (slopecross == 0){
            printf("colinear\n");
            return false;
        }

        T ACeta = etaC - etaA;
        T ACphi = phiC - phiA;
        T AC_ABcross = ACeta*ABphi - ACphi*ABeta;

        T invslopecross = 1/slopecross;
        T c = AC_ABcross*invslopecross;

        if (c < -shapetol || c > 1 + shapetol){
            printf("c = %g out of range\n", c);
            return false;
        }


        T AC_CDcross = ACeta*CDphi - ACphi*CDeta;
        T a = AC_CDcross*invslopecross;

        if (a < -shapetol || a > 1 + shapetol){
            printf("a = %g out of range\n", a);
            return false;
        }

        //now we know they do intersect,
        //the question is whether it's interesting

        if (std::abs(c-0.5) < shapetol){
            printf("c = %g close enough to 0.5\n", c);

            if(std::abs(a-0.5) < shapetol){
                printf("a = %g close enough to 0.5\n", a);
            } else if (std::abs(a-0.0) < shapetol ||
                       std::abs(a-1.0) < shapetol){
                printf("a = %g close enough to 0 or 1\n", a);
            }
        }

        return true;
    }

    template <typename T>
    void shapesInPlane(const T RL,
                       const T R2,
                       const T R3,
                       const T R4,
                       const T R5,
                       const T R6,
                       const unsigned qi0,
                       const unsigned qi1,
                       const unsigned qi2,
                       const unsigned qi3,
                       const resolvedInputs<T>& rin,
                       unsigned& shape_idx,
                       unsigned& r_idx,
                       unsigned& ct_idx){
        /*
         * A-B is RL
         * so C-D is the crossing line
         */
        T etaA = rin.etas[qi0]; 
        T etaB = rin.etas[qi1];
        T etaC = rin.etas[qi2];
        T etaD = rin.etas[qi3];
        T phiA = rin.phis[qi0];
        T phiB = rin.phis[qi1];
        T phiC = rin.phis[qi2];
        T phiD = rin.phis[qi3];

        //sort so that A is left of B
        if (phiA > phiB){
            std::swap(phiA, phiB);
            std::swap(etaA, etaB);
        }

        //sort so that C is left of D
        if (phiC > phiD){
            std::swap(phiC, phiD);
            std::swap(etaC, etaD);
        }

        //if the distance between A and B crosses the +-pi boundary, 
        //then we need to wrap around
        bool ABcrosses = phiB - phiA > M_PI;
        bool CDcrosses = phiD - phiC > M_PI;

        if(ABcrosses){
            phiB -= 2*M_PI;
            std::swap(phiA, phiB);
            std::swap(etaA, etaB);
        }

        if(CDcrosses){
            phiD -= 2*M_PI;
            std::swap(phiC, phiD);
            std::swap(etaC, etaD);
        }

        //we also need to worry whether C-D would cross, 
        //if only it were in the right phi +-2pi range
        //the only solution I can think of is the check three times...
        //good news is that if one crosses, then the others will not,
        //so we can terminate early

        std::array<T, 3> offsets = {0, 2*M_PI, -2*M_PI};
        for(const T& offset : offsets){
            if(linesCross(etaA, phiA+offset,
                          etaB, phiB+offset,
                          etaC, phiC,
                          etaD, phiD,
                          rin.shapetol)){
                printf("intersection\n");
                break;
            } else {
                printf("no interesection\n");
            }
        }
    }

    template <typename T>
    bool equilateralTriangle(const T RL,
                             const T R2,
                             const T R3,
                             const T R4,
                             const T R5,
                             const T R6,
                             const unsigned qi0,
                             const unsigned qi1,
                             const unsigned qi2,
                             const unsigned qi3,
                             const resolvedInputs<T>& rin,
                             unsigned& shape_idx,
                             unsigned& r_idx,
                             unsigned& ct_idx){
        double Rtol = rin.shapetol * RL;

        printf("RL = %g\n", RL);
        printf("R2 = %g\n", R2);
        printf("R3 = %g\n", R3);
        printf("R4 = %g\n", R4);
        printf("R5 = %g\n", R5);
        printf("R6 = %g\n", R6);
        //situation 1: triangle = (0, 1, 2)
        if (RL - R2 < Rtol && RL - R4 < Rtol){
            printf("triangle case 1:\n");

            shape_idx = 3; //equilateral triangle

            std::array<std::pair<int, T>, 3> Rs;
            Rs[0] = {0, R3};
            Rs[1] = {1, R5}; 
            Rs[2] = {2, R6};
            
            std::pair<int, T> minDist = *std::min_element(Rs.begin(), Rs.end(),
                    [](const std::pair<int, T>& a, const std::pair<int, T>& b){
                    return a.second < b.second;});

            T r = minDist.second / RL;
            
            T start_eta, start_phi, point_eta, point_phi, opp_eta, opp_phi;
            if (minDist.first == 0){
                start_eta = rin.etas[qi0];
                start_phi = rin.phis[qi0];
                point_eta = rin.etas[qi3];
                point_phi = rin.phis[qi3];
                opp_eta = 0.5 * (rin.etas[qi1] + rin.etas[qi2]);
                opp_phi = 0.5 * (rin.phis[qi1] + rin.phis[qi2]);
            } else if (minDist.first == 1){
                start_eta = rin.etas[qi1];
                start_phi = rin.phis[qi1];
                point_eta = rin.etas[qi3];
                point_phi = rin.phis[qi3];
                opp_eta = 0.5 * (rin.etas[qi0] + rin.etas[qi2]);
                opp_phi = 0.5 * (rin.phis[qi0] + rin.phis[qi2]);
            } else { //minDist.first == 2
                start_eta = rin.etas[qi2];
                start_phi = rin.phis[qi2];
                point_eta = rin.etas[qi3];
                point_phi = rin.phis[qi3];
                opp_eta = 0.5 * (rin.etas[qi0] + rin.etas[qi1]);
                opp_phi = 0.5 * (rin.phis[qi0] + rin.phis[qi1]);
            }

            T dphi_opp = opp_phi - start_phi;
            T deta_opp = opp_eta - start_eta;
            T dist_opp = std::sqrt(dphi_opp * dphi_opp + deta_opp * deta_opp);

            T dot = deta_opp * (point_eta - start_eta) +
                    dphi_opp * (point_phi - start_phi);
            T ct = dot / (minDist.second * dist_opp);

            if(ct > 1.0){
                ct = 1.0;
            }

            T theta = std::acos(ct);

            printf("r = %g\n", r);
            printf("ct = %g\n", ct);
            printf("theta = %g\n", std::acos(ct));

            r_idx = static_cast<unsigned>(rin.r_triangle->index(r) + 1);
            ct_idx = static_cast<unsigned>(rin.ct_triangle->index(theta) + 1);

            return true;
        }

        //situation 2: triangle = (0, 1, 3)
        if (RL - R3 < Rtol && RL - R5 < Rtol){
            printf("triangle case 2:\n");
            shape_idx = 3; //equilateral triangle

            std::array<std::pair<int, T>, 3> Rs;
            Rs[0] = {0, R2};
            Rs[1] = {1, R4}; 
            Rs[2] = {2, R6};
            
            std::pair<int, T> minDist = *std::min_element(Rs.begin(), Rs.end(),
                    [](const std::pair<int, T>& a, const std::pair<int, T>& b){
                    return a.second < b.second;});

            T r = minDist.second / RL;
            
            T start_eta, start_phi, point_eta, point_phi, opp_eta, opp_phi;
            if (minDist.first == 0){
                start_eta = rin.etas[qi0];
                start_phi = rin.phis[qi0];
                point_eta = rin.etas[qi2];
                point_phi = rin.phis[qi2];
                opp_eta = 0.5 * (rin.etas[qi1] + rin.etas[qi3]);
                opp_phi = 0.5 * (rin.phis[qi1] + rin.phis[qi3]);
            } else if (minDist.first == 1){
                start_eta = rin.etas[qi1];
                start_phi = rin.phis[qi1];
                point_eta = rin.etas[qi2];
                point_phi = rin.phis[qi2];
                opp_eta = 0.5 * (rin.etas[qi0] + rin.etas[qi3]);
                opp_phi = 0.5 * (rin.phis[qi0] + rin.phis[qi3]);
            } else { //minDist.first == 2
                start_eta = rin.etas[qi3];
                start_phi = rin.phis[qi3];
                point_eta = rin.etas[qi2];
                point_phi = rin.phis[qi2];
                opp_eta = 0.5 * (rin.etas[qi0] + rin.etas[qi1]);
                opp_phi = 0.5 * (rin.phis[qi0] + rin.phis[qi1]);
            }

            T dphi_opp = opp_phi - start_phi;
            T deta_opp = opp_eta - start_eta;
            T dist_opp = std::sqrt(dphi_opp * dphi_opp + deta_opp * deta_opp);

            T dot = deta_opp * (point_eta - start_eta) +
                    dphi_opp * (point_phi - start_phi);
            T ct = dot / (minDist.second * dist_opp);

            if(ct > 1.0){
                ct = 1.0;
            }

            T theta = std::acos(ct);

            printf("r = %g\n", r);
            printf("ct = %g\n", ct);
            printf("theta = %g\n", std::acos(ct));

            r_idx = static_cast<unsigned>(rin.r_triangle->index(r) + 1);
            ct_idx = static_cast<unsigned>(rin.ct_triangle->index(theta) + 1);

            return true;
        }

        return false;
    }

    template <typename T>
    void resolved4(const unsigned qi0, 
                   const unsigned qi1,
                   const unsigned qi2,
                   const unsigned qi3,
                   const resolvedInputs<T>& rin,
                   unsigned& shape_idx,
                   unsigned& RL_idx,
                   unsigned& r_idx,
                   unsigned& ct_idx){

        printf("resolved4 %u %u %u %u\n", qi0, qi1, qi2, qi3);
        T RL = rin.floatDRs[qi0][qi1];
        RL_idx = static_cast<unsigned>(rin.coarseRL->index(RL) + 1);

        T R2 = rin.floatDRs[qi0][qi2];
        T R3 = rin.floatDRs[qi0][qi3];
        T R4 = rin.floatDRs[qi1][qi2];
        T R5 = rin.floatDRs[qi1][qi3];
        T R6 = rin.floatDRs[qi2][qi3];

        //check if all four points are distinct
        if(R6 ==0 || R5 == 0 || R4 == 0 || R3==0 || R2==0 || RL == 0){
            shape_idx = 0;//no useful shape if any distances are zero
            r_idx = 0;
            ct_idx = 0;

            printf("\tone of the distances is zero\n");
            printf("\tshape_idx = 0\n");

            return;
        } 

        //check if we make an equilateral triangle
        
        if(equilateralTriangle(RL, R2, R3, R4, R5, R6, 
                               qi0, qi1, qi2, qi3,
                               rin, 
                               shape_idx, r_idx, ct_idx)){
            printf("\tshape_idx = 3\n");
            return;
        }
        printf("\tit's not a triangle...\n");

        shapesInPlane(RL, R2, R3, R4, R5, R6, 
                      qi0, qi1, qi2, qi3,
                      rin, 
                      shape_idx, r_idx, ct_idx);
    }
};

#endif
