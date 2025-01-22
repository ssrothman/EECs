#ifndef EECS_FAST_RES4_H
#define EECS_FAST_RES4_H

#include "util.h"
#include "fastStructs.h"
#include "SRothman/SimonTools/src/deltaR.h"

static constexpr double triangleL = 5;
static constexpr double triangleM = 4;
static constexpr double triangleS = 3;

static constexpr double triangleLM2 = (triangleL/triangleM)*(triangleL/triangleM);
static constexpr double triangleLS2 = (triangleL/triangleS)*(triangleL/triangleS);

namespace fastEEC{
    /*
     * Given points (1, 2, 3, 4) in the eta/phi plane:
     *      1. Check if (1, 2, 3) form a 3,4,5 right triangle
     *      2. If they do NOT return false
     *          In this case R, r, theta are untouched
     *      3. If they do, return true, and
     *          set R = the length of the hypotenuse 
     *          set r = the distance between point (4) 
     *                  and the right-angle vertex
     *                  divided by R
     *          set theta = the angle between segment r
     *                      and the short leg of the triangle
     * 
     *  The domains of the results are:
     *      R on [0, inf]
     *      r on [0, inf]
     *      theta in [0, 2pi]
     */
    template <typename T>
    bool checkTriangle(const float shapetol,
            const T eta1,
            const T eta2,
            const T eta3,
            const T eta4,
            const T phi1,
            const T phi2,
            const T phi3,
            const T phi4,
            T& R,
            T& r, 
            T& theta) noexcept {

        T dR2_12 = dR2(eta1, phi1, eta2, phi2);
        T dR2_13 = dR2(eta1, phi1, eta3, phi3);
        T dR2_23 = dR2(eta2, phi2, eta3, phi3);

        std::array<std::pair<char, T>, 3> dR2s = {{
            {0, dR2_12},
            {1, dR2_13},
            {2, dR2_23}
        }};
        std::sort(dR2s.begin(), dR2s.end(), [](std::pair<char, T> a, std::pair<char, T> b){
                return a.second > b.second;
        });

        T RL2 = dR2s[0].second;
        T RM2 = dR2s[1].second;
        T RS2 = dR2s[2].second;

        if (std::abs(RL2/RM2 - triangleLM2) > shapetol){
            return false;
        } else if(std::abs(RL2/RS2 - triangleLS2) > shapetol){
            return false;
        } else {
            /*
             * We pass the tringle selection
             *
             * Now we find A, B, C such that:
             *      AB is the "3" leg
             *      AC is the "4" leg
             *      BC is the "5" leg
             */
            T etaA=0, etaB=0, etaC=0;
            T phiA=0, phiB=0, phiC=0;
            
            switch(dR2s[0].first){
                case 0:
                    if(dR2s[1].first == 1){ 
                        //RL = (1,2), RM = (1,3)
                        //So A = 3, B = 2 C = 1
                        etaA = eta3;
                        phiA = phi3;

                        etaB = eta2;
                        phiB = phi2;

                        etaC = eta1;
                        phiC = phi1;
                    } else{
                        //RL = (1,2), RM = (2,3)
                        //so A = 3, B = 1, C = 2
                        etaA = eta3;
                        phiA = phi3;

                        etaB = eta1;
                        phiB = phi1;

                        etaC = eta2;
                        phiC = phi2;
                    }
                    break;
                case 1:
                    if(dR2s[1].first == 0){
                        //RL = (1,3), RM = (1,2)
                        //So A = 2, B = 3, C = 1
                        etaA = eta2;
                        phiA = phi2;

                        etaB = eta3;
                        phiB = phi3;

                        etaC = eta1;
                        phiC = phi1;
                    } else{
                        //RL = (1,3), RM = (2,3)
                        //So A = 2, B = 1, C = 3
                        etaA = eta2;
                        phiA = phi2;

                        etaB = eta1;
                        phiB = phi1;

                        etaC = eta3;
                        phiC = phi3;
                    }
                    break;
                case 2:
                    if (dR2s[1].first == 0){
                        //RL = (2,3), RM = (1,2)
                        //So A = 1, B = 3, C = 2
                        etaA = eta1;
                        phiA = phi1;

                        etaB = eta3;
                        phiB = phi3;

                        etaC = eta2;
                        phiC = phi2;
                    } else {
                        //RL = (2, 3), RM = (1,3)
                        //So A = 1, B = 2, C = 3
                        etaA = eta1;
                        phiA = phi1;

                        etaB = eta2;
                        phiB = phi2;

                        etaC = eta3;
                        phiC = phi3;
                    }
                    break;
                default:
                    assert(false);
            }

            T deta_A4 = eta4 - etaA;
            T dphi_A4 = deltaphi(phiA, phi4);

            const T dR2_A4 = deta_A4*deta_A4 + dphi_A4*dphi_A4;

            R = std::sqrt(RL2);
            T rA4 = std::sqrt(dR2_A4);

            r = rA4/R;

            T deta_AB = etaB-etaA;
            T dphi_AB = deltaphi(phiA, phiB);

            T RAB = std::sqrt(RS2);
            //T RAC = std::sqrt(RM2);

            T AB_dot_A4 = deta_AB * deta_A4 + dphi_AB * dphi_A4;
            T costheta = AB_dot_A4/(rA4 * RAB);

            //gives the acos in the range [0, pi]
            T acos_costheta = std::acos(costheta);

            /*
             * To check the quadrant, we check the sign of A4 dot AC
             * If it's potive we're in the upper half of the plane
             * If it's negative we're in the lower half of the plane
             */
            T deta_AC = etaC-etaA;
            T dphi_AC = deltaphi(phiA, phiC);
            T AC_dot_A4 = deta_AC * deta_A4 + dphi_AC * dphi_A4;

            if (AC_dot_A4 < 0){
                theta = 2*M_PI-acos_costheta;
            } else {
                theta = acos_costheta;
            }
                
            return true;
        }
    }

    /*
     * Given four points (A, B, C, D)
     * Check whether the pairs (A, B), (C, D) form crossing segments
     *  in either the dipole or tee configurations
     *
     * PRECONDITION:
     *  deltaphi(C, D) < pi
     *  deltaphi(A, BdetaTeeB*detaTeeB + dphiTeeB*dphiTeeB;) < pi
     *  deltaR(C, D) < deltaR(A, B)
     *
     * RESULT:
     *  if dR(centroidAB, centroidCD) < shapetol:
     *      we have a "dipole"
     *      return 1
     *  else if dR(centroidCD, A) < shapetol or dR(centroidCD, B) < shapetol:
     *      we have a "tee"
     *      return 2
     *  else:
     *      noshape
     *      return 0
     */
    template <typename T>
    unsigned linesCross(const float shapetol,
                    const T etaA,
                    const T etaB,
                    const T etaC,
                    const T etaD,
                    const T phiA,
                    const T phiB,
                    const T phiC,
                    const T phiD) noexcept {
        //centroid of CD
        T etaCD = 0.5 * (etaC + etaD);
        T phiCD = 0.5 * (phiC + phiD);

        //centroid of AB
        T etaAB = 0.5 * (etaA + etaB);
        T phiAB = 0.5 * (phiA + phiB);

        //check for dipole configuration
        T dR2Centroid = dR2(etaAB, phiAB, etaCD, phiCD);

        if (dR2Centroid < shapetol*shapetol){
            return 1;
        }

        //check for tee configuration
        T dR2TeeA = dR2(etaA, phiA, etaCD, phiCD);

        if (dR2TeeA < shapetol*shapetol){
            return 2;
        }

        const T dR2TeeB = dR2(etaB, phiB, etaCD, phiCD);

        if (dR2TeeB < shapetol*shapetol){
            return 2;
        }

        return 0;
    }

    /*
     * Given that (A, B) + (C, D) is an interesting shape
     *  compute the binning variables r and theta
     *
     * PRECONDITION:
     *  (A, B) + (C, D) is either a tee or a dipole
     *  RL = deltaR(A, B)
     *  RS = deltaR(C, D)
     *  RL > RS
     *
     * RESULT:
     *  set r = deltaR(C, D) / deltaR(A, B)
     *  set theta = angle between AB and CD
     *
     * The domain of the results is:
     *  r on [0, 1] due to the precondition
     *  theta on [0, pi/2] due to symmetry
     */
    template <typename T>
    void getrandangle(const T etaA, const T etaB, const T etaC, const T etaD,
                      const T phiA, const T phiB, const T phiC, const T phiD,
                      const T RL, const T RS,
                      T& r, T& theta) noexcept {

        T etaAB = etaB - etaA;
        T phiAB = deltaphi(phiA, phiB);

        T etaCD = etaD - etaC;
        T phiCD = deltaphi(phiC, phiD);

        T dot = etaAB*etaCD + phiAB*phiCD;

        T cosTheta = dot/(RL*RS);
        
        r = RS/RL;
        if (cosTheta > 1){
            if (cosTheta > 1.000001){
                printf("I hope this doesn't actually ever happen...\n");
                printf("\tetaAB = %f, phiAB = %f\n", etaAB, phiAB);
                printf("\tetaCD = %f, phiCD = %f\n", etaCD, phiCD);
                printf("\tRL = %f, RS = %f\n", RL, RS);
                printf("\tdot = %f, cosTheta = %f\n", dot, cosTheta);
            }
            theta = 0;
        } else {
            theta = std::acos(cosTheta);
        }

        if (theta > M_PI/2){
            theta = M_PI - theta;
        }
    }

    /*
     * Checks the set of four points pointed to by next
     *  with permutation set by next.ires4
     *  and looks for tees and dipoles
     *
     * If we find a dipole or tee we:
     *      set shapeindex, RL, r, theta indices correctly inside next
     *      return the shapeindex
     * else:
     *      set shapeindex=0 inside next
     *      return 0
     */
    template <typename T>
    unsigned lookForShapes(const jetDetails_t<T>& jetDetails,
                           const res4shapesAxes_t& res4ax,
                           prev_t<T, 5>& next){

        T etaA, etaB, etaC, etaD;
        T phiA, phiB, phiC, phiD;
        T RAB, RCD;

        //set AB, CD RAB, RCD such that
        //RAB > RCD
        //The coice of AB, CD from the four points
        //is determined by next.ires4
        getEtasPhis(jetDetails, next, 
                etaA, etaB, etaC, etaD, 
                phiA, phiB, phiC, phiD, 
                RAB, RCD);

        //reject configurations with too small shorter legs
        if(RCD < 1e-8){
            next.shape_res4_idx[next.ires4] = 0;
            next.RL_res4_idx[next.ires4] = 0;
            next.r_res4_idx[next.ires4] = 0;
            next.ct_res4_idx[next.ires4] = 0;

            return 0;
        }

        //get the shape index
        unsigned shape = linesCross(res4ax.shapetol,
                                    etaA, etaB, etaC, etaD,
                                    phiA, phiB, phiC, phiD);

        if(shape != 0){
            //found a tee or a dipole

            //look up r, theta
            T r, theta;
            getrandangle(etaA, etaB, etaC, etaD,
                         phiA, phiB, phiC, phiD,
                         RAB, RCD,
                         r, theta);

            //get RL index, which is same for all
            next.RL_res4_idx[next.ires4] = getIndex(RAB, res4ax.RL);

            if(shape == 1){
                //it's a dipole
                next.shape_res4_idx[next.ires4] = 1;
                next.r_res4_idx[next.ires4] = getIndex(r, res4ax.r_dipole);
                next.ct_res4_idx[next.ires4] = getIndex(theta, res4ax.ct_dipole);
            } else {
                //it's a tee
                next.shape_res4_idx[next.ires4] = 2;
                next.r_res4_idx[next.ires4] = getIndex(r, res4ax.r_tee);
                next.ct_res4_idx[next.ires4] = getIndex(theta, res4ax.ct_tee);
            }
            return shape;
        } else {
            //its not a tee or a dipole
            next.shape_res4_idx[next.ires4] = 0;
            return 0;
        }
    }

    /*
     * Look up etas phis, and distances pointed to by next
     * with premutation pointed to by next.ires4
     *
     * Result gaurentees that:
     *      RAB > RCD
     *      phiC < phiD
     *      phiA < phiB
     *      deltaphi(C, D) < pi
     *      deltaphi(A, B) < pi
     */
    template <typename T>
    void getEtasPhis(const jetDetails_t<T>& jetDetails,
                     const prev_t<T, 5>& prev,

                     T& etaA, T& etaB, T& etaC, T& etaD,
                     T& phiA, T& phiB, T& phiC, T& phiD,

                     T& RAB, T& RCD) noexcept {
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

    /*
     * Run res4 without checking the shape - for transfer calculations
     * This implementation should probably not be trusted
     * And hasn't been looked at in ages
     */
    template <typename T, int genOrder>
    void runRes4_nocheck(const jetDetails_t<T>& jetDetails,
                         const res4shapesAxes_t& res4ax,

                         const prev_t<T, genOrder>& prevGen,

                         prev_t<T, 5>& next) noexcept {

        

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

    /*
     * Entry point for res4 calculations
     * Checks for tees, dipoles and triangles
     * And stores the results in next
     */
    template <typename T>
    void runRes4(const jetDetails_t<T>& jetDetails,
                 const res4shapesAxes_t& res4ax,

                 prev_t<T, 5>& next) noexcept {
        //Look up deltaRs
        std::array<T, 6> dRs = {{
            jetDetails.floatDRs[next.is[0]][next.is[1]],
            jetDetails.floatDRs[next.is[0]][next.is[2]],
            jetDetails.floatDRs[next.is[0]][next.is[3]],
            jetDetails.floatDRs[next.is[1]][next.is[2]],
            jetDetails.floatDRs[next.is[1]][next.is[3]],
            jetDetails.floatDRs[next.is[2]][next.is[3]]
        }};

        //reject configurations where two points are the same
        auto minel = min_element(dRs.begin(), dRs.end());
        if (*minel < 1e-16){
            return;
        }

        //try all three tee/dipole permutations
        for(next.ires4=0; next.ires4<3; ++next.ires4){
            lookForShapes(jetDetails, res4ax, next);
        }

        //for the triangle there are actually four valid permutations
        //these arrays set that up
        std::array<T, 4> etas = {{
            jetDetails.etas[next.is[0]],
            jetDetails.etas[next.is[1]],
            jetDetails.etas[next.is[2]],
            jetDetails.etas[next.is[3]]
        }};
        std::array<T, 4> phis = {{
            jetDetails.phis[next.is[0]],
            jetDetails.phis[next.is[1]],
            jetDetails.phis[next.is[2]],
            jetDetails.phis[next.is[3]]
        }};

        static constexpr std::array<unsigned, 4> As = {{0, 1, 2, 3}};
        static constexpr std::array<unsigned, 4> Bs = {{1, 2, 3, 0}};
        static constexpr std::array<unsigned, 4> Cs = {{2, 3, 0, 1}};
        static constexpr std::array<unsigned, 4> Ds = {{3, 0, 1, 2}};

        //loop over all the triangle permutations
        for(unsigned iTri=0; iTri<4; ++iTri){
            T etaA = etas[As[iTri]];
            T etaB = etas[Bs[iTri]];
            T etaC = etas[Cs[iTri]];
            T etaD = etas[Ds[iTri]];

            T phiA = phis[As[iTri]];
            T phiB = phis[Bs[iTri]];
            T phiC = phis[Cs[iTri]];
            T phiD = phis[Ds[iTri]];

            T R, r, theta;
            //check for triangle
            if (checkTriangle(res4ax.shapetol,
                        etaA, etaB, etaC, etaD,
                        phiA, phiB, phiC, phiD,
                        R, r, theta)){

                next.istri_res4[iTri] = true;
                next.RL_res4tri_idx[iTri] = getIndex(R, res4ax.RL);
                next.r_res4tri_idx[iTri] = getIndex(r, res4ax.r_triangle);
                next.ct_res4tri_idx[iTri] = getIndex(theta, res4ax.ct_triangle);
            } else {
                next.istri_res4[iTri] = false;
            }
        }
    } 
}

#endif
