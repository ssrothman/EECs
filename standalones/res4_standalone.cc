#include "res4_standalone.h"
#include "adjacency.h"

#include "SRothman/SimonTools/src/deltaR.h"
#include "SRothman/SimonTools/src/histutil.h"
#include "SRothman/SimonTools/src/util.h"

#include <array>

//this is in principle configurable
//but it is important that the target ratios not be 1
constexpr double triangle_ratio_LM = 5.0/4.0;
constexpr double triangle_ratio_MS = 4.0/3.0;

/*
 * Compute angle between vectors (12) and (34)
 * the precomputed detas, dphis, and dRs are passed
 * 
 * the result is stored into c
 *      and is in the range [0, pi/2]
 */
static void compute_c(
        const double deta12, const double dphi12,
        const double deta34, const double dphi34,
        const double dR12, const double dR34,
        double& c) noexcept {

    double dot = deta12*deta34 + dphi12*dphi34;
    double ct = dot/(dR12*dR34);

    if (ct > 1){
        ct = 1;
    } else if (ct < -1){
        ct = -1;
    }
    c = std::acos(ct);
    /*
     * std::acos yields a result in [0, pi]
     * the symmetry of the problem actually puts c in [0, pi/2]
     * so we can compress the range
     */
    if(c > M_PI/2.0){
        c = M_PI - c;
    }
}

/*
 * Given two distances dR12 and dR34, compute R and r such that:
 *      R = max(dR12, dR34)
 *      r = min(dR12, dR34)/R
 */
static void compute_r_R(
        const double dR12,
        const double dR34,
        double& r,
        double& R) noexcept {

    if(dR12 > dR34){
        R = dR12;
        r = dR34/dR12;
    } else {
        R = dR34;
        r = dR12/dR34;
    }
}

/*
 * Compute the midpoint of two polar angles,
 * wrapping around as necessary
 * This is useful for the dipole and tee axes
 */
static double angle_midpoint(
        const double phi1, 
        const double phi2) noexcept {
    double result;

    /*
     * I wonder if this could be optimized slightly...
     * but I am confident that it works
     */

    if (phi2-phi1 > M_PI){
        result = (phi1 + phi2 - 2*M_PI)/2.0;
    } else if(phi2-phi1 < -M_PI){
        result = (phi1 + phi2 + 2*M_PI)/2.0;
    } else {
        return (phi1 + phi2)/2.0;
    }
    if (result > M_PI){
        result -= 2*M_PI;
    } else if (result < -M_PI){
        result += 2*M_PI;
    }
    return result;
}

struct res4_entry{
    unsigned idx_R, idx_r, idx_c;
    bool isShape=false;
};

/*
 * Check a given permutation of four particles 
 * for tee and dipole configurations
 *
 * If they are found, the result is stored into ans
 *
 * NOTE: particle order matters. 
 *     We check the configuration (12)(34)
 *
 * A boolean template parameter tracks whether the 
 * distances passed are squared or not
 */
template <class ResultType, bool distances_squared, bool actually_fill>
static void check_tee_dipole(
        ResultType& ans,

        res4_entry& dipole_entry,
        res4_entry& tee_entry,

        const double eta1, 
        const double phi1,
        const double eta2, 
        const double phi2,
        const double eta3, 
        const double phi3,
        const double eta4, 
        const double phi4,

        const double dR12_2,
        const double dR34_2,

        const double deta12,
        const double dphi12,
        const double deta34,
        const double dphi34,

        const double wt,

        const double tolerance2,

        const standaloneEEC::axis& ax_R,
        const standaloneEEC::axis& ax_r_dipole,
        const standaloneEEC::axis& ax_c_dipole,
        const standaloneEEC::axis& ax_r_tee,
        const standaloneEEC::axis& ax_c_tee) noexcept {

    /*
     * Check the configuration (12)(23) for tees and dipoles
     */

    bool isDipole, isTee;

    double mideta12 = (eta1 + eta2)/2.0;
    double midphi12 = angle_midpoint(phi1, phi2);

    double mideta34 = (eta3 + eta4)/2.0;
    double midphi34 = angle_midpoint(phi3, phi4);

    double dR2_midpoints = deltaR2(mideta12, midphi12,
                               mideta34, midphi34);
    isDipole = dR2_midpoints < tolerance2;

    if (dR12_2 < dR34_2){
        double dR2_mid12_3 = deltaR2(mideta12, midphi12,
                                eta3, phi3);
        double dR2_mid12_4 = deltaR2(mideta12, midphi12,
                                eta4, phi4);

        isTee = (dR2_mid12_3 < tolerance2) 
             || (dR2_mid12_4 < tolerance2);
    } else {
        double dR2_mid34_1 = deltaR2(mideta34, midphi34,
                                eta1, phi1);
        double dR2_mid34_2 = deltaR2(mideta34, midphi34,
                                eta2, phi2);

        isTee = (dR2_mid34_1 < tolerance2) 
             || (dR2_mid34_2 < tolerance2);
    }

    if(isTee || isDipole){
        double dR12, dR34;
        if constexpr(distances_squared){
            dR12 = std::sqrt(dR12_2);
            dR34 = std::sqrt(dR34_2);
        } else {
            dR12 = dR12_2;
            dR34 = dR34_2;
        }

        double R,r,c;
        compute_r_R(dR12, dR34, r, R);
        compute_c(deta12, dphi12,
                  deta34, dphi34,
                  dR12, dR34,
                  c);

        unsigned idx_R = getIndex(R, ax_R);
        if(isDipole){
            unsigned idx_r = getIndex(r, ax_r_dipole);
            unsigned idx_c = getIndex(c, ax_c_dipole);
            if constexpr (actually_fill){
                ans.fill_dipole(idx_R, idx_r, idx_c, wt);
            }
            dipole_entry.idx_R = idx_R;
            dipole_entry.idx_r = idx_r;
            dipole_entry.idx_c = idx_c;
            dipole_entry.isShape = true;
        }
        if (isTee){
            unsigned idx_r = getIndex(r, ax_r_tee);
            unsigned idx_c = getIndex(c, ax_c_tee);
            if constexpr (actually_fill){
                ans.fill_tee(idx_R, idx_r, idx_c, wt);
            }
            tee_entry.idx_R = idx_R;
            tee_entry.idx_r = idx_r;
            tee_entry.idx_c = idx_c;
            tee_entry.isShape = true;
        }
    }
}

/*
 * Check a given permutation of four particles
 * for triangle configurations
 *
 * If they are found, the result is stored into ans
 *
 * NOTE: particle order matters
 *    We check the configuration (123)(4)
 *    But order in the (123) does NOT matter
 *    And the function checks all 6 permutations
 *
 * A boolean template parameter tracks whether the 
 * distances passed are squared or not
 */
template <class ResultType, bool distances_squared, bool actually_fill>
static void check_triangle(
        ResultType& ans,

        res4_entry& triangle_entry,

        const double dR12_2, const double dR13_2,
        const double dR14_2, const double dR23_2,
        const double dR24_2, const double dR34_2,

        [[maybe_unused]] const double deta12, 
        [[maybe_unused]] const double dphi12,
        [[maybe_unused]] const double deta13, 
        [[maybe_unused]] const double dphi13,
        [[maybe_unused]] const double deta14, 
        [[maybe_unused]] const double dphi14,
        [[maybe_unused]] const double deta23, 
        [[maybe_unused]] const double dphi23,
        [[maybe_unused]] const double deta24, 
        [[maybe_unused]] const double dphi24,
        [[maybe_unused]] const double deta34, 
        [[maybe_unused]] const double dphi34,

        const double wt,

        const double tri_tolerance,

        const standaloneEEC::axis& ax_R,
        const standaloneEEC::axis& ax_r_triangle,
        const standaloneEEC::axis& ax_c_triangle) noexcept {

    /*
     * Check the configuration (123)(4) for a triangle
     */

    /*
     * First step is to find out the order of triangle sides
     * RL > RM > RS
     * 
     * We need to also keep track of which side is which
     * so that we can place the correct distance and angles
     * of the fourth particle w.r.t. the right angle
     *
     * There are ultimately 3! = 6 permutations of RL,RM,RS
     * which is annoying lol
     *
     * I can't think of a better way to do this than to 
     * just have a messy if-else block
     * but I'm sure a better way exists :(
     *
     * IF-ELSE BLOCK PLANS;
     * want to find A, B, C such that BC > AC > AB
     * ie BC is the hypotenuse and AB is the smallest side
     * this also means that A is the right-angle corner
     */

    double RAB_2, RAC_2, RBC_2;

    double RA4_2; //for r
    //for c dot products
    double detaA4, dphiA4; 
    double detaAC, dphiAC;
    double detaAB, dphiAB;

    if(dR12_2 > dR13_2){
        if(dR13_2 > dR23_2){
            // (12) > (13) > (23)
            // so A = 3, B = 2, C = 1
            RAB_2 = dR23_2;
            RAC_2 = dR13_2;
            RBC_2 = dR12_2;

            RA4_2 = dR34_2;
            
            detaA4 = deta34;
            dphiA4 = dphi34;

            detaAC = deta13;
            dphiAC = dphi13;

            detaAB = deta23;
            dphiAB = dphi23;

        } else if (dR12_2 > dR23_2){
            // (12) > (23) > (13)
            // so A = 3, B = 1, C = 2
            RAB_2 = dR13_2;
            RAC_2 = dR23_2;
            RBC_2 = dR12_2;

            RA4_2 = dR34_2;

            detaA4 = deta34;
            dphiA4 = dphi34;

            detaAC = deta23;
            dphiAC = dphi23;

            detaAB = deta13;
            dphiAB = dphi13;

        } else {
            // (23) > (12) > (13)
            // so A = 1, B = 3, C = 2
            RAB_2 = dR13_2;
            RAC_2 = dR12_2;
            RBC_2 = dR23_2;

            RA4_2 = dR14_2;

            detaA4 = deta14;
            dphiA4 = dphi14;

            detaAC = deta12;
            dphiAC = dphi12;

            detaAB = deta23;
            dphiAB = dphi23;

        }
    } else{
        if(dR12_2 > dR23_2){ 
            // (13) > (12) > (23)
            // so A = 2, B = 3, C = 1
            RAB_2 = dR23_2;
            RAC_2 = dR12_2;
            RBC_2 = dR13_2;

            RA4_2 = dR24_2;

            detaA4 = deta24;
            dphiA4 = dphi24;

            detaAC = deta12;
            dphiAC = dphi12;

            detaAB = deta23;
            dphiAB = dphi23;

        } else if(dR13_2 > dR23_2){ 
            // (13) > (23) > (12)
            // so A = 2, B = 1, C = 3
            RAB_2 = dR12_2;
            RAC_2 = dR23_2;
            RBC_2 = dR13_2;

            RA4_2 = dR24_2;

            detaA4 = deta24;
            dphiA4 = dphi24;

            detaAC = deta23;
            dphiAC = dphi23;

            detaAB = deta12;
            dphiAB = dphi12;

        } else {
            // (23) > (13) > (12)
            // so A = 1, B = 2, C = 3
            RAB_2 = dR12_2;
            RAC_2 = dR13_2;
            RBC_2 = dR23_2;

            RA4_2 = dR14_2;

            detaA4 = deta14;
            dphiA4 = dphi14;

            detaAC = deta13;
            dphiAC = dphi13;

            detaAB = deta12;
            dphiAB = dphi12;
        }
    }

    /*
     * Now we have the order of the triangle sides
     * we can compute the angles and distances
     */
    bool passLM; 
    bool passMS; 
    
    if constexpr(distances_squared){
        passLM = std::abs(std::sqrt(RBC_2/RAC_2) - triangle_ratio_LM) < tri_tolerance;
        passMS = std::abs(std::sqrt(RAC_2/RAB_2) - triangle_ratio_MS) < tri_tolerance;
    } else {
        passLM = std::abs(RBC_2/RAC_2 - triangle_ratio_LM) < tri_tolerance;
        passMS = std::abs(RAC_2/RAB_2 - triangle_ratio_MS) < tri_tolerance;
    }

    if(passLM && passMS){
        double R, r, c;

        double RAC, RA4;
        if constexpr(distances_squared){
            RAC = std::sqrt(RAC_2);
            RA4 = std::sqrt(RA4_2);
        } else{
            RAC = RAC_2;
            RA4 = RA4_2;
        }

        R = RAC;
        r = RA4/RAC;

        double dotA4AC = detaA4*detaAC + dphiA4*dphiAC;

        double cosc = dotA4AC/(RA4*RAC);
        if (cosc > 1){
            cosc = 1;
        } else if (cosc < -1){
            cosc = -1;
        }
        c = std::acos(cosc);
        
        /*
         * std::acos() yields a result in [0, pi]
         * in this case there is no symmetry
         * and c properly lives in [0, 2pi]
         *
         * we therefore need to check the other dot product 
         * to place c in the correct quadrant
         */

        double dotA4AB = detaA4*detaAB + dphiA4*dphiAB;

        if (dotA4AB < 0){
            c = -c;
        }

        unsigned idx_R = getIndex(R, ax_R);
        unsigned idx_r = getIndex(r, ax_r_triangle);
        unsigned idx_c = getIndex(c, ax_c_triangle);
        if constexpr(actually_fill){
            ans.fill_triangle(idx_R, idx_r, idx_c, wt);
        }
        triangle_entry.idx_R = idx_R;
        triangle_entry.idx_r = idx_r;
        triangle_entry.idx_c = idx_c;
        triangle_entry.isShape = true;
    }
}

template <class ResultType, bool distances_squared, bool actually_fill>
static void innermost_level(
        ResultType& ans,

        std::array<res4_entry, 3>& dipole_entries,
        std::array<res4_entry, 3>& tee_entries,
        std::array<res4_entry, 4>& triangle_entries,

        const double eta1, const double phi1,
        const double eta2, const double phi2,
        const double eta3, const double phi3, 
        const double eta4, const double phi4,

        const double deta12, const double dphi12,
        const double deta13, const double dphi13,
        const double deta14, const double dphi14,
        const double deta23, const double dphi23,
        const double deta24, const double dphi24,
        const double deta34, const double dphi34,

        const double dR12_2, const double dR13_2, 
        const double dR14_2, const double dR23_2, 
        const double dR24_2, const double dR34_2,

        const double wt,

        const standaloneEEC::axis& ax_R,
        const standaloneEEC::axis& ax_r_dipole,
        const standaloneEEC::axis& ax_c_dipole,
        const standaloneEEC::axis& ax_r_tee,
        const standaloneEEC::axis& ax_c_tee,
        const standaloneEEC::axis& ax_r_triangle,
        const standaloneEEC::axis& ax_c_triangle,

        const double tolerance2,
        const double tri_tolerance) noexcept {
    /*
     * Now we check the three permutations
     * of pairs of four particles
     * for tees and dipoles
     *
     * NB we do these checks together 
     * because they share the 
     * r, c, R calculations
     *
     * The permuatations are:
     * (12)(34)
     * (13)(24)
     * (14)(23)
     */

    //(12)(34)
    check_tee_dipole<ResultType, distances_squared, actually_fill>(
            ans,
            dipole_entries[0], 
            tee_entries[0],
            eta1, phi1,
            eta2, phi2,
            eta3, phi3,
            eta4, phi4,
            dR12_2, dR34_2,
            deta12, dphi12,
            deta34, dphi34,
            wt,
            tolerance2,
            ax_R,
            ax_r_dipole,
            ax_c_dipole,
            ax_r_tee,
            ax_c_tee);

    //(13)(24)
    check_tee_dipole<ResultType, distances_squared, actually_fill>(
            ans,
            dipole_entries[1],
            tee_entries[1],
            eta1, phi1,
            eta3, phi3,
            eta2, phi2,
            eta4, phi4,
            dR13_2, dR24_2,
            deta13, dphi13,
            deta24, dphi24,
            wt,
            tolerance2,
            ax_R,
            ax_r_dipole,
            ax_c_dipole,
            ax_r_tee,
            ax_c_tee);

    //(14)(23)
    check_tee_dipole<ResultType, distances_squared, actually_fill>(
            ans,
            dipole_entries[2],
            tee_entries[2],
            eta1, phi1,
            eta4, phi4,
            eta2, phi2,
            eta3, phi3,
            dR14_2, dR23_2,
            deta14, dphi14,
            deta23, dphi23,
            wt,
            tolerance2,
            ax_R,
            ax_r_dipole,
            ax_c_dipole,
            ax_r_tee,
            ax_c_tee);

    /*
     * Now we need to check for triangles
     * There are four triangle configurations
     * (ie 4 choose 3)
     */

    //(123)(4)
    check_triangle<ResultType, distances_squared, actually_fill>(
            ans,
            triangle_entries[0],

            dR12_2, dR13_2, 
            dR14_2, dR23_2,
            dR24_2, dR34_2,

            deta12, dphi12,
            deta13, dphi13,
            deta14, dphi14,
            deta23, dphi23,
            deta24, dphi24,
            deta34, dphi34,

            wt,

            tri_tolerance,

            ax_R,
            ax_r_triangle,
            ax_c_triangle);

    //(124)(3)
    check_triangle<ResultType, distances_squared, actually_fill>(
            ans,
            triangle_entries[1],

            dR12_2, dR14_2, 
            dR13_2, dR24_2,
            dR23_2, dR34_2,

            deta12, dphi12,
            deta14, dphi14,
            deta13, dphi13,
            deta24, dphi24,
            deta23, dphi23,
            deta34, dphi34,

            wt,

            tri_tolerance,

            ax_R,
            ax_r_triangle,
            ax_c_triangle);

    //(143)(2)
    check_triangle<ResultType, distances_squared, actually_fill>(
            ans,
            triangle_entries[2],

            dR14_2, dR13_2, 
            dR12_2, dR34_2,
            dR24_2, dR23_2,

            deta14, dphi14,
            deta13, dphi13,
            deta12, dphi12,
            deta34, dphi34,
            deta24, dphi24,
            deta23, dphi23,

            wt,

            tri_tolerance,

            ax_R,
            ax_r_triangle,
            ax_c_triangle);

    //(423)(1)
    check_triangle<ResultType, distances_squared, actually_fill>(
            ans,
            triangle_entries[3],

            dR24_2, dR34_2, 
            dR14_2, dR23_2,
            dR12_2, dR13_2,

            deta24, dphi24,
            deta34, dphi34,
            deta14, dphi14,
            deta23, dphi23,
            deta12, dphi12,
            deta13, dphi13,

            wt,

            tri_tolerance,

            ax_R,
            ax_r_triangle,
            ax_c_triangle);
}

/*
 * This is the main loop over four particles
 *
 * NB in the specific case of res4 we are only interested
 * in the case where the four particles are distinct
 * This allows us to save some loops
 * This also means we don't have to worry about 
 * symmetry factors
 *
 * Probably an unnecessary mount of effort was made
 * to avoid ANY repeated computations
 *
 * Thus all of the particle coordinates and pairwise 
 * quantities are computed at the top of each loop level
 * and passed down to inner loops
 *
 * We opt to only precompute the /squared/ distances
 * and leave the square root for later.
 * In principle this can lead to repeated square roots.
 * However, this is more than made up for by the fact
 * that the vast majority of cases are not an interesting
 * configuration, so we can avoid the square root altogether
 *
 * We template on the container type to allow for two 
 * versions of this method to share the same code
 * without any costly virtual function calls.
 * The two versions are:
 *     Filling a histogram directly via a multi_array
 *     Growing a vector of entries with bin indices and weights
 */
template <class BasicContainer>
static void res4_mainloop(
        standaloneEEC::res4_result<BasicContainer>& ans,

        const standaloneEEC::EECjet& thisjet,
        const standaloneEEC::axis& ax_R,
        const standaloneEEC::axis& ax_r_dipole,
        const standaloneEEC::axis& ax_c_dipole,
        const standaloneEEC::axis& ax_r_tee,
        const standaloneEEC::axis& ax_c_tee,
        const standaloneEEC::axis& ax_r_triangle,
        const standaloneEEC::axis& ax_c_triangle,

        const double tolerance2,
        const double tri_tolerance) noexcept {

    for (unsigned i1=0; i1<thisjet.N; ++i1){
        const auto&[E1, eta1, phi1] = thisjet.singles[i1];

        for (unsigned i2=i1+1; i2<thisjet.N; ++i2){
            const auto&[E2, eta2, phi2] = thisjet.singles[i2];

            const double E12 = E1 * E2;

            const double deta12 = deltaEta(eta1, eta2);
            const double dphi12 = deltaPhi(phi1, phi2);
            const double dR12_2 = square(deta12) + square(dphi12);

            for(unsigned i3=i2+1; i3<thisjet.N; ++i3){
                const auto&[E3, eta3, phi3] = thisjet.singles[i3];

                const double E123 = E12 * E3;

                const double deta13 = deltaEta(eta1, eta3);
                const double dphi13 = deltaPhi(phi1, phi3);
                const double dR13_2 = square(deta13) + square(dphi13);

                const double deta23 = deltaEta(eta2, eta3);
                const double dphi23 = deltaPhi(phi2, phi3);
                const double dR23_2 = square(deta23) + square(dphi23);

                for(unsigned i4=i3+1; i4<thisjet.N; ++i4){
                    const auto&[E4, eta4, phi4] = thisjet.singles[i4];

                    const double wt = E123 * E4;

                    const double deta14 = deltaEta(eta1, eta4);
                    const double dphi14 = deltaPhi(phi1, phi4);
                    const double dR14_2 = square(deta14) 
                        + square(dphi14);

                    const double deta24 = deltaEta(eta2, eta4);
                    const double dphi24 = deltaPhi(phi2, phi4);
                    const double dR24_2 = square(deta24) 
                        + square(dphi24);

                    const double deta34 = deltaEta(eta3, eta4);
                    const double dphi34 = deltaPhi(phi3, phi4);
                    const double dR34_2 = square(deta34) 
                        + square(dphi34);

                    std::array<res4_entry, 3> dipole_entries;
                    std::array<res4_entry, 3> tee_entries;
                    std::array<res4_entry, 4> triangle_entries;

                    innermost_level<standaloneEEC::res4_result<BasicContainer>, true, true>(
                            ans,

                            dipole_entries,
                            tee_entries,
                            triangle_entries,

                            eta1, phi1,
                            eta2, phi2,
                            eta3, phi3,
                            eta4, phi4,

                            deta12, dphi12,
                            deta13, dphi13,
                            deta14, dphi14,
                            deta23, dphi23,
                            deta24, dphi24,
                            deta34, dphi34,

                            dR12_2, dR13_2,
                            dR14_2, dR23_2,
                            dR24_2, dR34_2,

                            wt,

                            ax_R,
                            ax_r_dipole,
                            ax_c_dipole,
                            ax_r_tee,
                            ax_c_tee,
                            ax_r_triangle,
                            ax_c_triangle,

                            tolerance2,
                            tri_tolerance);
                }//end loop over i4
            }//end loop over i3
        }//end loop over i2
    }//end loop over i1
}//end res4_standalone()
 
template <class BasicContainer>
static void res4_mainloop_precomputed(
        standaloneEEC::res4_result<BasicContainer>& ans,

        const standaloneEEC::EECjet_precomputed& thisjet,
        const standaloneEEC::axis& ax_R,
        const standaloneEEC::axis& ax_r_dipole,
        const standaloneEEC::axis& ax_c_dipole,
        const standaloneEEC::axis& ax_r_tee,
        const standaloneEEC::axis& ax_c_tee,
        const standaloneEEC::axis& ax_r_triangle,
        const standaloneEEC::axis& ax_c_triangle,

        const double tolerance2,
        const double tri_tolerance) noexcept {

    for (unsigned i1=0; i1<thisjet.N; ++i1){
        const auto&[E1, eta1, phi1] = thisjet.singles[i1];

        for (unsigned i2=i1+1; i2<thisjet.N; ++i2){
            const auto&[E2, eta2, phi2] = thisjet.singles[i2];

            const double E12 = E1 * E2;

            const auto&[deta12, dphi12, dR12] = thisjet.pairs[i1][i2];

            for(unsigned i3=i2+1; i3<thisjet.N; ++i3){
                const auto&[E3, eta3, phi3] = thisjet.singles[i3];

                const double E123 = E12 * E3;

                const auto&[deta13, dphi13, dR13] = thisjet.pairs[i1][i3];
                const auto&[deta23, dphi23, dR23] = thisjet.pairs[i2][i3];

                for(unsigned i4=i3+1; i4<thisjet.N; ++i4){
                    const auto&[E4, eta4, phi4] = thisjet.singles[i4];

                    const double wt = E123 * E4;

                    const auto&[deta14, dphi14, dR14] = thisjet.pairs[i1][i4];
                    const auto&[deta24, dphi24, dR24] = thisjet.pairs[i2][i4];
                    const auto&[deta34, dphi34, dR34] = thisjet.pairs[i3][i4];

                    std::array<res4_entry, 3> dipole_entries;
                    std::array<res4_entry, 3> tee_entries;
                    std::array<res4_entry, 4> triangle_entries;

                    innermost_level<standaloneEEC::res4_result<BasicContainer>, false, true>(
                            ans,

                            dipole_entries,
                            tee_entries,
                            triangle_entries,

                            eta1, phi1,
                            eta2, phi2,
                            eta3, phi3,
                            eta4, phi4,

                            deta12, dphi12,
                            deta13, dphi13,
                            deta14, dphi14,
                            deta23, dphi23,
                            deta24, dphi24,
                            deta34, dphi34,

                            dR12, dR13,
                            dR14, dR23,
                            dR24, dR34,

                            wt,

                            ax_R,
                            ax_r_dipole,
                            ax_c_dipole,
                            ax_r_tee,
                            ax_c_tee,
                            ax_r_triangle,
                            ax_c_triangle,

                            tolerance2,
                            tri_tolerance);
                }//end loop over i4
            }//end loop over i3
        }//end loop over i2
    }//end loop over i1
}//end res4_standalone()

/*
 * This is the entry point for the res4 calculation 
 * with the multi_array container backend
 *
 * In this mode the result is accumulated directly into a
 * multi_array object representing the histogram. 
 *
 * Pros: dense representation is easier to work with
 *       fixed maximum size in memory
 *       no need for dynamic allocation to grow vectors
 *
 * Cons: less flexible
 *       can't benefit from any sparsity in memory usage
 *       requires seeking through large arrays to fill bins
 * 
 * I intent to test explicitely which backend is better,
 * but for now I expose both
 */
void standaloneEEC::res4_standalone_multi_array(
        standaloneEEC::res4_result_multi_array& ans,
        
        const simon_jet& J,
        const normType& nt,

        const standaloneEEC::axis& ax_R,
        const standaloneEEC::axis& ax_r_dipole,
        const standaloneEEC::axis& ax_c_dipole,
        const standaloneEEC::axis& ax_r_tee,
        const standaloneEEC::axis& ax_c_tee,
        const standaloneEEC::axis& ax_r_triangle,
        const standaloneEEC::axis& ax_c_triangle,

        const double tolerance,
        const double tri_tolerance) noexcept {

    standaloneEEC::EECjet thisjet(J, nt);

    res4_mainloop<res4_multi_array_container>(ans, thisjet, 
            ax_R, 
            ax_r_dipole, ax_c_dipole, 
            ax_r_tee, ax_c_tee, 
            ax_r_triangle, ax_c_triangle,
            tolerance*tolerance,
            tri_tolerance);
}

void standaloneEEC::res4_standalone_multi_array_precomputed(
        standaloneEEC::res4_result_multi_array& ans,
        
        const simon_jet& J,
        const normType& nt,

        const standaloneEEC::axis& ax_R,
        const standaloneEEC::axis& ax_r_dipole,
        const standaloneEEC::axis& ax_c_dipole,
        const standaloneEEC::axis& ax_r_tee,
        const standaloneEEC::axis& ax_c_tee,
        const standaloneEEC::axis& ax_r_triangle,
        const standaloneEEC::axis& ax_c_triangle,

        const double tolerance,
        const double tri_tolerance) noexcept {

    standaloneEEC::EECjet_precomputed thisjet(J, nt);

    res4_mainloop_precomputed<res4_multi_array_container>(ans, thisjet, 
            ax_R, 
            ax_r_dipole, ax_c_dipole, 
            ax_r_tee, ax_c_tee, 
            ax_r_triangle, ax_c_triangle,
            tolerance*tolerance,
            tri_tolerance);
}
 
/*
 * This is the entry point for the res4 calculation 
 * with the vector container backend
 *
 * In this mode the result is accumulated directly into a
 * multi_array object representing the histogram. 
 *
 * Cons: dense representation is easier to work with
 *       no maximum size in memory
 *       need for dynamic allocation to grow vectors
 *
 * Pros: more flexible
 *       benefits from sparsity in memory usage
 *       no seeking through large arrays to fill bins
 */
void standaloneEEC::res4_standalone_vector(
        standaloneEEC::res4_result_vector& ans,
        
        const simon_jet& J,
        const normType& nt,

        const standaloneEEC::axis& ax_R,
        const standaloneEEC::axis& ax_r_dipole,
        const standaloneEEC::axis& ax_c_dipole,
        const standaloneEEC::axis& ax_r_tee,
        const standaloneEEC::axis& ax_c_tee,
        const standaloneEEC::axis& ax_r_triangle,
        const standaloneEEC::axis& ax_c_triangle,

        const double tolerance,
        const double tri_tolerance) noexcept {

    standaloneEEC::EECjet thisjet(J, nt);

    res4_mainloop<res4_vector_container>(ans, thisjet, 
            ax_R, 
            ax_r_dipole, ax_c_dipole, 
            ax_r_tee, ax_c_tee, 
            ax_r_triangle, ax_c_triangle,
            tolerance*tolerance,
            tri_tolerance);
}

void standaloneEEC::res4_standalone_vector_precomputed(
        standaloneEEC::res4_result_vector& ans,
        
        const simon_jet& J,
        const normType& nt,

        const standaloneEEC::axis& ax_R,
        const standaloneEEC::axis& ax_r_dipole,
        const standaloneEEC::axis& ax_c_dipole,
        const standaloneEEC::axis& ax_r_tee,
        const standaloneEEC::axis& ax_c_tee,
        const standaloneEEC::axis& ax_r_triangle,
        const standaloneEEC::axis& ax_c_triangle,

        const double tolerance,
        const double tri_tolerance) noexcept {

    standaloneEEC::EECjet_precomputed thisjet(J, nt);

    res4_mainloop_precomputed<res4_vector_container>(ans, thisjet, 
            ax_R, 
            ax_r_dipole, ax_c_dipole, 
            ax_r_tee, ax_c_tee, 
            ax_r_triangle, ax_c_triangle,
            tolerance*tolerance,
            tri_tolerance);
}

template <class TransferContainer, class BasicContainer>
static void res4_transferloop(
        standaloneEEC::res4_transfer_result<TransferContainer, BasicContainer>& ans,

        std::array<res4_entry, 3>& dipole_entries_1,
        std::array<res4_entry, 3>& tee_entries_1,
        std::array<res4_entry, 4>& triangle_entries_1,

        const standaloneEEC::EECjet_precomputed& otherjet,

        const standaloneEEC::neighborhood* const n1,
        const standaloneEEC::neighborhood* const n2,
        const standaloneEEC::neighborhood* const n3,
        const standaloneEEC::neighborhood* const n4,

        const standaloneEEC::axis& ax_R,
        const standaloneEEC::axis& ax_r_dipole,
        const standaloneEEC::axis& ax_c_dipole,
        const standaloneEEC::axis& ax_r_tee,
        const standaloneEEC::axis& ax_c_tee,
        const standaloneEEC::axis& ax_r_triangle,
        const standaloneEEC::axis& ax_c_triangle,

        const double tolerance2,
        const double tri_tolerance,

        const double wt) noexcept {

    for(const standaloneEEC::neighbor& j1: *n1){
        const double twt1 = wt * j1.wt;
        const auto&[E1, eta1, phi1] = otherjet.singles[j1.idx];

        for(const standaloneEEC::neighbor& j2: *n2){
            const double twt2 = twt1 * j2.wt;
            const auto&[E2, eta2, phi2] = otherjet.singles[j2.idx];
            
            const auto&[deta12, dphi12, dR12] = otherjet.pairs[j1.idx][j2.idx];

            for(const standaloneEEC::neighbor& j3 : *n3){
                const double twt3 = twt2 * j3.wt;
                const auto&[E3, eta3, phi3] = otherjet.singles[j3.idx];

                const auto&[deta13, dphi13, dR13] = otherjet.pairs[j1.idx][j3.idx];
                const auto&[deta23, dphi23, dR23] = otherjet.pairs[j2.idx][j3.idx];

                for(const standaloneEEC::neighbor& j4 : *n4){
                    const double twt4 = twt3 * j4.wt;
                    const auto&[E4, eta4, phi4] = otherjet.singles[j4.idx];

                    const auto&[deta14, dphi14, dR14] = otherjet.pairs[j1.idx][j4.idx];
                    const auto&[deta24, dphi24, dR24] = otherjet.pairs[j2.idx][j4.idx];
                    const auto&[deta34, dphi34, dR34] = otherjet.pairs[j3.idx][j4.idx];

                    std::array<res4_entry, 3> dipole_entries_2;
                    std::array<res4_entry, 3> tee_entries_2;
                    std::array<res4_entry, 4> triangle_entries_2;

                    innermost_level<standaloneEEC::res4_transfer_result<TransferContainer, BasicContainer>, false, false>(
                            ans,

                            dipole_entries_2,
                            tee_entries_2,
                            triangle_entries_2,

                            eta1, phi1,
                            eta2, phi2,
                            eta3, phi3,
                            eta4, phi4,

                            deta12, dphi12,
                            deta13, dphi13,
                            deta14, dphi14,
                            deta23, dphi23,
                            deta24, dphi24,
                            deta34, dphi34,

                            dR12, dR13,
                            dR14, dR23,
                            dR24, dR34,

                            twt4,

                            ax_R,
                            ax_r_dipole,
                            ax_c_dipole,
                            ax_r_tee,
                            ax_c_tee,
                            ax_r_triangle,
                            ax_c_triangle,

                            tolerance2,
                            tri_tolerance);

                    for(unsigned i=0; i<3; ++i){
                        if(dipole_entries_1[i].isShape && dipole_entries_2[i].isShape){
                            ans.fill_dipole(dipole_entries_1[i].idx_R, 
                                            dipole_entries_1[i].idx_r, 
                                            dipole_entries_1[i].idx_c,
                                            dipole_entries_2[i].idx_R, 
                                            dipole_entries_2[i].idx_r, 
                                            dipole_entries_2[i].idx_c,
                                            twt4);
                        } else if(dipole_entries_1[i].isShape){
                            ans.unmatched1.fill_dipole(
                                    dipole_entries_1[i].idx_R, 
                                    dipole_entries_1[i].idx_r, 
                                    dipole_entries_1[i].idx_c,
                                    twt4);
                        } else if(dipole_entries_2[i].isShape){
                            ans.unmatched2.fill_dipole(
                                    dipole_entries_2[i].idx_R, 
                                    dipole_entries_2[i].idx_r, 
                                    dipole_entries_2[i].idx_c,
                                    twt4);
                        }//end if dipole entries match

                        if(tee_entries_1[i].isShape && tee_entries_2[i].isShape){
                            ans.fill_tee(tee_entries_1[i].idx_R, 
                                         tee_entries_1[i].idx_r, 
                                         tee_entries_1[i].idx_c,
                                         tee_entries_2[i].idx_R, 
                                         tee_entries_2[i].idx_r, 
                                         tee_entries_2[i].idx_c,
                                         twt4);
                        } else if(tee_entries_1[i].isShape){
                            ans.unmatched1.fill_tee(
                                    tee_entries_1[i].idx_R, 
                                    tee_entries_1[i].idx_r, 
                                    tee_entries_1[i].idx_c,
                                    twt4);
                        } else if(tee_entries_2[i].isShape){
                            ans.unmatched2.fill_tee(
                                    tee_entries_2[i].idx_R, 
                                    tee_entries_2[i].idx_r, 
                                    tee_entries_2[i].idx_c,
                                    twt4);
                        }//end if tee entries match
                    }//end loop over tee/dipole entries

                    for(unsigned i=0; i<4; ++i){
                        if(triangle_entries_1[i].isShape && triangle_entries_2[i].isShape){
                            ans.fill_triangle(triangle_entries_1[i].idx_R, 
                                              triangle_entries_1[i].idx_r, 
                                              triangle_entries_1[i].idx_c,
                                              triangle_entries_2[i].idx_R, 
                                              triangle_entries_2[i].idx_r, 
                                              triangle_entries_2[i].idx_c,
                                              twt4);
                        } else if(triangle_entries_1[i].isShape){
                            ans.unmatched1.fill_triangle(
                                    triangle_entries_1[i].idx_R, 
                                    triangle_entries_1[i].idx_r, 
                                    triangle_entries_1[i].idx_c,
                                    twt4);
                        } else if(triangle_entries_2[i].isShape){
                            ans.unmatched2.fill_triangle(
                                    triangle_entries_2[i].idx_R, 
                                    triangle_entries_2[i].idx_r, 
                                    triangle_entries_2[i].idx_c,
                                    twt4);
                        }//end if triangle entries match
                    }//end loop over triangle entries
                }//end loop over j4
            }//end loop over j3
        }//end loop over j2
    }//end loop over j1
}//end res4_transferloop()

template <class BasicContainer, class TransferContainer>
static void res4_mainloop_transfer_precomputed(
        standaloneEEC::res4_result<BasicContainer>& ans,
        standaloneEEC::res4_transfer_result<TransferContainer, BasicContainer>& transfer_ans,

        const standaloneEEC::EECjet_precomputed& thisjet,
        const standaloneEEC::EECjet_precomputed& otherjet,

        const standaloneEEC::adjacency& adj,

        const standaloneEEC::axis& ax_R,
        const standaloneEEC::axis& ax_r_dipole,
        const standaloneEEC::axis& ax_c_dipole,
        const standaloneEEC::axis& ax_r_tee,
        const standaloneEEC::axis& ax_c_tee,
        const standaloneEEC::axis& ax_r_triangle,
        const standaloneEEC::axis& ax_c_triangle,

        const double tolerance2,
        const double tri_tolerance) noexcept {

    for (unsigned i1=0; i1<thisjet.N; ++i1){
        const auto&[E1, eta1, phi1] = thisjet.singles[i1];
        const bool hasMatch1 = adj.hasMatch[i1];
        const standaloneEEC::neighborhood* const n1 = hasMatch1 ? &adj.adj[i1] : nullptr;

        for (unsigned i2=i1+1; i2<thisjet.N; ++i2){
            const auto&[E2, eta2, phi2] = thisjet.singles[i2];
            const bool hasMatch2 = hasMatch1 && adj.hasMatch[i2];
            const standaloneEEC::neighborhood* const n2 = hasMatch2 ? &adj.adj[i2] : nullptr;

            const double E12 = E1 * E2;

            const auto&[deta12, dphi12, dR12] = thisjet.pairs[i1][i2];

            for(unsigned i3=i2+1; i3<thisjet.N; ++i3){
                const auto&[E3, eta3, phi3] = thisjet.singles[i3];
                const bool hasMatch3 = hasMatch2 && adj.hasMatch[i3];
                const standaloneEEC::neighborhood* const n3 = hasMatch3 ? &adj.adj[i3] : nullptr;

                const double E123 = E12 * E3;

                const auto&[deta13, dphi13, dR13] = thisjet.pairs[i1][i3];
                const auto&[deta23, dphi23, dR23] = thisjet.pairs[i2][i3];

                for(unsigned i4=i3+1; i4<thisjet.N; ++i4){
                    const auto&[E4, eta4, phi4] = thisjet.singles[i4];
                    const bool hasMatch4 = hasMatch3 && adj.hasMatch[i4];
                    const standaloneEEC::neighborhood* const n4 = hasMatch4 ? &adj.adj[i4] : nullptr;

                    const double wt = E123 * E4;

                    const auto&[deta14, dphi14, dR14] = thisjet.pairs[i1][i4];
                    const auto&[deta24, dphi24, dR24] = thisjet.pairs[i2][i4];
                    const auto&[deta34, dphi34, dR34] = thisjet.pairs[i3][i4];

                    std::array<res4_entry, 3> dipole_entries;
                    std::array<res4_entry, 3> tee_entries;
                    std::array<res4_entry, 4> triangle_entries;

                    innermost_level<standaloneEEC::res4_result<BasicContainer>, false, true>(
                            ans,

                            dipole_entries,
                            tee_entries,
                            triangle_entries,

                            eta1, phi1,
                            eta2, phi2,
                            eta3, phi3,
                            eta4, phi4,

                            deta12, dphi12,
                            deta13, dphi13,
                            deta14, dphi14,
                            deta23, dphi23,
                            deta24, dphi24,
                            deta34, dphi34,

                            dR12, dR13,
                            dR14, dR23,
                            dR24, dR34,

                            wt,

                            ax_R,
                            ax_r_dipole,
                            ax_c_dipole,
                            ax_r_tee,
                            ax_c_tee,
                            ax_r_triangle,
                            ax_c_triangle,

                            tolerance2,
                            tri_tolerance);

                    if(hasMatch4){
                        res4_transferloop<TransferContainer, BasicContainer>(
                            transfer_ans,
                            dipole_entries,
                            tee_entries,
                            triangle_entries,
                            otherjet,
                            n1, n2, n3, n4,
                            ax_R,
                            ax_r_dipole,
                            ax_c_dipole,
                            ax_r_tee,
                            ax_c_tee,
                            ax_r_triangle,
                            ax_c_triangle,
                            tolerance2,
                            tri_tolerance,
                            wt);
                    }
                }//end loop over i4
            }//end loop over i3
        }//end loop over i2
    }//end loop over i1
}//end res4_standalone()

void standaloneEEC::res4_standalone_transfer_multi_array(
        standaloneEEC::res4_result_multi_array& ans,
        standaloneEEC::res4_transfer_result_multi_array& transfer_ans,

        const simon_jet& J1,
        const simon_jet& J2,
        const normType& nt,

        const Eigen::MatrixXd& adjmat,

        const standaloneEEC::axis& ax_R,
        const standaloneEEC::axis& ax_r_dipole,
        const standaloneEEC::axis& ax_c_dipole,
        const standaloneEEC::axis& ax_r_tee,
        const standaloneEEC::axis& ax_c_tee,
        const standaloneEEC::axis& ax_r_triangle,
        const standaloneEEC::axis& ax_c_triangle,

        const double tolerance,
        const double tri_tolerance) noexcept {

    standaloneEEC::EECjet_precomputed thisjet(J1, nt);
    standaloneEEC::EECjet_precomputed otherjet(J2, nt);

    standaloneEEC::adjacency adj(adjmat);

    res4_mainloop_transfer_precomputed<res4_multi_array_container, res4_transfer_multi_array_container>(
            ans, transfer_ans,
            thisjet, otherjet,
            adj,
            ax_R,
            ax_r_dipole, ax_c_dipole,
            ax_r_tee, ax_c_tee,
            ax_r_triangle, ax_c_triangle,
            tolerance*tolerance,
            tri_tolerance);
}

void standaloneEEC::res4_standalone_transfer_vector(
        standaloneEEC::res4_result_multi_array& ans,
        standaloneEEC::res4_transfer_result_vector& transfer_ans,

        const simon_jet& J1,
        const simon_jet& J2,
        const normType& nt,

        const Eigen::MatrixXd& adjmat,

        const standaloneEEC::axis& ax_R,
        const standaloneEEC::axis& ax_r_dipole,
        const standaloneEEC::axis& ax_c_dipole,
        const standaloneEEC::axis& ax_r_tee,
        const standaloneEEC::axis& ax_c_tee,
        const standaloneEEC::axis& ax_r_triangle,
        const standaloneEEC::axis& ax_c_triangle,

        const double tolerance,
        const double tri_tolerance) noexcept {

    standaloneEEC::EECjet_precomputed thisjet(J1, nt);
    standaloneEEC::EECjet_precomputed otherjet(J2, nt);

    standaloneEEC::adjacency adj(adjmat);

    res4_mainloop_transfer_precomputed<res4_multi_array_container, res4_transfer_vector_container>(
            ans, transfer_ans,
            thisjet, otherjet,
            adj,
            ax_R,
            ax_r_dipole, ax_c_dipole,
            ax_r_tee, ax_c_tee,
            ax_r_triangle, ax_c_triangle,
            tolerance*tolerance,
            tri_tolerance);
}
