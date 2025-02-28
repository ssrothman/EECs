#ifndef SROTHMAN_EEC_V2_RES4_BACKEND_H
#define SROTHMAN_EEC_V2_RES4_BACKEND_H

#include "Adjacency.h"

#include "SRothman/SimonTools/src/deltaR.h"
#include "SRothman/SimonTools/src/histutil.h"
#include "SRothman/SimonTools/src/util.h"

#include <memory>

#include <array>

/*
 * Compute angle between vectors (12) and (34)
 * the precomputed detas, dphis, and dRs are passed
 * 
 * the result is stored into c
 *      and is in the range [0, pi/2]
 */
inline void compute_c(
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
inline void compute_r_R(
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
inline double angle_midpoint(
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
    bool isShape;
    res4_entry() noexcept : idx_R(0), idx_r(0), idx_c(0), isShape(false) {}
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
inline void check_tee_dipole(
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

        const EEC::Res4Axes axes) noexcept {

    /*
     * Check the configuration (12)(23) for tees and dipoles
     */

    bool isDipole, isTee;

    double mideta12 = (eta1 + eta2)/2.0;
    double midphi12 = angle_midpoint(phi1, phi2);

    double mideta34 = (eta3 + eta4)/2.0;
    double midphi34 = angle_midpoint(phi3, phi4);

    double dR2_midpoints = simon::deltaR2(mideta12, midphi12,
                               mideta34, midphi34);
    isDipole = dR2_midpoints < tolerance2;

    if (dR12_2 < dR34_2){
        double dR2_mid12_3 = simon::deltaR2(mideta12, midphi12,
                                eta3, phi3);
        double dR2_mid12_4 = simon::deltaR2(mideta12, midphi12,
                                eta4, phi4);

        isTee = (dR2_mid12_3 < tolerance2) 
             || (dR2_mid12_4 < tolerance2);
    } else {
        double dR2_mid34_1 = simon::deltaR2(mideta34, midphi34,
                                eta1, phi1);
        double dR2_mid34_2 = simon::deltaR2(mideta34, midphi34,
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

        unsigned idx_R = simon::getIndex(R, axes.R);
        if(isDipole){
            unsigned idx_r = simon::getIndex(r, axes.r_dipole);
            unsigned idx_c = simon::getIndex(c, axes.c_dipole);
            if constexpr (actually_fill){
                ans.fill_dipole(idx_R, idx_r, idx_c, wt);
            }
            dipole_entry.idx_R = idx_R;
            dipole_entry.idx_r = idx_r;
            dipole_entry.idx_c = idx_c;
            dipole_entry.isShape = true;
        } else{
            dipole_entry.isShape = false;
        }
        if (isTee){
            unsigned idx_r = simon::getIndex(r, axes.r_tee);
            unsigned idx_c = simon::getIndex(c, axes.c_tee);
            if constexpr (actually_fill){
                ans.fill_tee(idx_R, idx_r, idx_c, wt);
            }
            tee_entry.idx_R = idx_R;
            tee_entry.idx_r = idx_r;
            tee_entry.idx_c = idx_c;
            tee_entry.isShape = true;
        } else {
            tee_entry.isShape = false;
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
inline void check_triangle(
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

        const EEC::Res4Axes& axes) noexcept {

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
        passLM = std::abs(std::sqrt(RBC_2/RAC_2) - EEC::TRIANGLE_RATIO_LM) < tri_tolerance;
        passMS = std::abs(std::sqrt(RAC_2/RAB_2) - EEC::TRIANGLE_RATIO_MS) < tri_tolerance;
    } else {
        passLM = std::abs(RBC_2/RAC_2 - EEC::TRIANGLE_RATIO_LM) < tri_tolerance;
        passMS = std::abs(RAC_2/RAB_2 - EEC::TRIANGLE_RATIO_MS) < tri_tolerance;
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

        unsigned idx_R = simon::getIndex(R, axes.R);
        unsigned idx_r = simon::getIndex(r, axes.r_triangle);
        unsigned idx_c = simon::getIndex(c, axes.c_triangle);
        if constexpr(actually_fill){
            ans.fill_triangle(idx_R, idx_r, idx_c, wt);
        }
        triangle_entry.idx_R = idx_R;
        triangle_entry.idx_r = idx_r;
        triangle_entry.idx_c = idx_c;
        triangle_entry.isShape = true;
    } else{
        triangle_entry.isShape = false;
    }
}

template <class ResultType, bool distances_squared, bool actually_fill>
inline void innermost_level(
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

        const EEC::Res4Axes& axes,

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
            axes);

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
            axes);

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
            axes);

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

            axes);

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

            axes);

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

            axes);

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

            axes);
}

template <class TransferResultType, class ResultType, class JetType>
inline void res4_transferloop(
        TransferResultType& ans,
        ResultType& untransfered_reco,
        ResultType& untransfered_gen,

        const std::array<res4_entry, 3>& dipole_entries_gen,
        const std::array<res4_entry, 3>& tee_entries_gen,
        const std::array<res4_entry, 4>& triangle_entries_gen,

        const std::shared_ptr<const JetType> thisjet_reco,

        const EEC::neighborhood& n1,
        const EEC::neighborhood& n2,
        const EEC::neighborhood& n3,
        const EEC::neighborhood& n4,

        const EEC::Res4Axes& axes_reco,

        const double tolerance2,
        const double tri_tolerance,

        const double wt_gen) noexcept {

    for(const EEC::neighbor& j1: n1){
        const double twt1 = wt_gen * j1.wt;
        //unused E1 throws compiler warnings
        //this seems to be an unavoidable language limitation
        //we have to wait until c++26 for a solution 
        //(https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2020/p2169r0.pdf)
        const auto&[E1, eta1, phi1] = thisjet_reco->singles.get(j1.idx);

        for(const EEC::neighbor& j2: n2){
            const double twt2 = twt1 * j2.wt;
            const auto&[E2, eta2, phi2] = thisjet_reco->singles.get(j2.idx);
            
            const auto&[deta12, dphi12, dR12] = thisjet_reco->pairs.get(j1.idx, j2.idx);

            for(const EEC::neighbor& j3 : n3){
                const double twt3 = twt2 * j3.wt;
                const auto&[E3, eta3, phi3] = thisjet_reco->singles.get(j3.idx);

                const auto&[deta13, dphi13, dR13] = thisjet_reco->pairs.get(j1.idx, j3.idx);
                const auto&[deta23, dphi23, dR23] = thisjet_reco->pairs.get(j2.idx, j3.idx);

                for(const EEC::neighbor& j4 : n4){
                    const double twt4 = twt3 * j4.wt;
                    const auto&[E4, eta4, phi4] = thisjet_reco->singles.get(j4.idx);

                    const auto&[deta14, dphi14, dR14] = thisjet_reco->pairs.get(j1.idx, j4.idx);
                    const auto&[deta24, dphi24, dR24] = thisjet_reco->pairs.get(j2.idx, j4.idx);
                    const auto&[deta34, dphi34, dR34] = thisjet_reco->pairs.get(j3.idx, j4.idx);

                    std::array<res4_entry, 3> dipole_entries_reco;
                    std::array<res4_entry, 3> tee_entries_reco;
                    std::array<res4_entry, 4> triangle_entries_reco;

                    innermost_level<TransferResultType, JetType::pairType::distances_squared, false>(
                            ans,

                            dipole_entries_reco,
                            tee_entries_reco,
                            triangle_entries_reco,

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

                            axes_reco,

                            tolerance2,
                            tri_tolerance);

                    for(unsigned i=0; i<3; ++i){
                        if(dipole_entries_gen[i].isShape && dipole_entries_reco[i].isShape){
                            ans.fill_dipole(dipole_entries_reco[i].idx_R, 
                                            dipole_entries_reco[i].idx_r, 
                                            dipole_entries_reco[i].idx_c,
                                            dipole_entries_gen[i].idx_R, 
                                            dipole_entries_gen[i].idx_r, 
                                            dipole_entries_gen[i].idx_c,
                                            twt4, wt_gen);
                        } else if(dipole_entries_reco[i].isShape){
                            untransfered_reco.fill_dipole(
                                    dipole_entries_reco[i].idx_R, 
                                    dipole_entries_reco[i].idx_r, 
                                    dipole_entries_reco[i].idx_c,
                                    twt4);
                        } else if(dipole_entries_gen[i].isShape){
                            untransfered_gen.fill_dipole(
                                    dipole_entries_gen[i].idx_R, 
                                    dipole_entries_gen[i].idx_r, 
                                    dipole_entries_gen[i].idx_c,
                                    twt4);
                        }//end if dipole entries match

                        if(tee_entries_gen[i].isShape && tee_entries_reco[i].isShape){
                            ans.fill_tee(tee_entries_reco[i].idx_R, 
                                         tee_entries_reco[i].idx_r, 
                                         tee_entries_reco[i].idx_c,
                                         tee_entries_gen[i].idx_R, 
                                         tee_entries_gen[i].idx_r, 
                                         tee_entries_gen[i].idx_c,
                                         twt4, wt_gen);
                        } else if(tee_entries_reco[i].isShape){
                            untransfered_reco.fill_tee(
                                    tee_entries_reco[i].idx_R, 
                                    tee_entries_reco[i].idx_r, 
                                    tee_entries_reco[i].idx_c,
                                    twt4);
                        } else if(tee_entries_gen[i].isShape){
                            untransfered_gen.fill_tee(
                                    tee_entries_gen[i].idx_R, 
                                    tee_entries_gen[i].idx_r, 
                                    tee_entries_gen[i].idx_c,
                                    twt4);
                        }//end if tee entries match
                    }//end loop over tee/dipole entries

                    for(unsigned i=0; i<4; ++i){
                        if(triangle_entries_gen[i].isShape && triangle_entries_reco[i].isShape){
                            ans.fill_triangle(triangle_entries_reco[i].idx_R, 
                                              triangle_entries_reco[i].idx_r, 
                                              triangle_entries_reco[i].idx_c,
                                              triangle_entries_gen[i].idx_R, 
                                              triangle_entries_gen[i].idx_r, 
                                              triangle_entries_gen[i].idx_c,
                                              twt4, wt_gen);
                        } else if(triangle_entries_reco[i].isShape){
                            untransfered_reco.fill_triangle(
                                    triangle_entries_reco[i].idx_R, 
                                    triangle_entries_reco[i].idx_r, 
                                    triangle_entries_reco[i].idx_c,
                                    twt4);
                        } else if(triangle_entries_gen[i].isShape){
                            untransfered_gen.fill_triangle(
                                    triangle_entries_gen[i].idx_R, 
                                    triangle_entries_gen[i].idx_r, 
                                    triangle_entries_gen[i].idx_c,
                                    twt4);
                        }//end if triangle entries match
                    }//end loop over triangle entries
                }//end loop over j4
            }//end loop over j3
        }//end loop over j2
    }//end loop over j1
}//end res4_transferloop()

template <class ResultType, class JetType, bool doUnmatched, class TransferResultType, bool doTransfer>
inline void res4_mainloop(
        ResultType& ans,
        ResultType* unmatched_gen,
        TransferResultType* transfer_ans,
        ResultType* untransfered_reco,
        ResultType* untransfered_gen,

        const std::shared_ptr<const JetType> thisjet_reco,
        const std::shared_ptr<const JetType> thisjet_gen,
        const std::vector<bool>* matched,

        const std::shared_ptr<const EEC::Adjacency> adj,

        const EEC::Res4Axes * const axes_reco,
        const EEC::Res4Axes& axes_gen,

        const double tolerance2,
        const double tri_tolerance) noexcept {

    for (unsigned i1=0; i1<thisjet_gen->N; ++i1){
        const auto&[E1, eta1, phi1] = thisjet_gen->singles.get(i1);
        [[maybe_unused]] bool matched1;
        if constexpr(doUnmatched){
            matched1 = matched->at(i1);
        }
        [[maybe_unused]] EEC::neighborhood const * n1;
        if constexpr(doTransfer){
            n1 = &(adj->get_neighborhood(i1));
        }

        for (unsigned i2=i1+1; i2<thisjet_gen->N; ++i2){
            const auto&[E2, eta2, phi2] = thisjet_gen->singles.get(i2);
            [[maybe_unused]] bool matched2;
            if constexpr(doUnmatched){
                matched2 = matched1 && matched->at(i2);
            }
            [[maybe_unused]] EEC::neighborhood const * n2;
            if constexpr(doTransfer){
                n2 = &(adj->get_neighborhood(i2));
            }

            const double E12 = E1 * E2;

            const auto&[deta12, dphi12, dR12] = thisjet_gen->pairs.get(i1, i2);

            for(unsigned i3=i2+1; i3<thisjet_gen->N; ++i3){
                const auto&[E3, eta3, phi3] = thisjet_gen->singles.get(i3);
                [[maybe_unused]] bool matched3;
                if constexpr(doUnmatched){
                    matched3 = matched2 && matched->at(i3);
                }
                [[maybe_unused]] EEC::neighborhood const * n3;
                if constexpr(doTransfer){
                    n3 = &(adj->get_neighborhood(i3));
                }

                const double E123 = E12 * E3;

                const auto&[deta13, dphi13, dR13] = thisjet_gen->pairs.get(i1, i3);
                const auto&[deta23, dphi23, dR23] = thisjet_gen->pairs.get(i2, i3);

                for(unsigned i4=i3+1; i4<thisjet_gen->N; ++i4){
                    const auto&[E4, eta4, phi4] = thisjet_gen->singles.get(i4);
                    [[maybe_unused]] bool matched4;
                    if constexpr(doUnmatched){
                        matched4 = matched3 && matched->at(i4);
                    }
                    [[maybe_unused]] EEC::neighborhood const * n4;
                    if constexpr(doTransfer){
                        n4 = &(adj->get_neighborhood(i4));
                    }

                    const double wt = E123 * E4;

                    const auto&[deta14, dphi14, dR14] = thisjet_gen->pairs.get(i1, i4);
                    const auto&[deta24, dphi24, dR24] = thisjet_gen->pairs.get(i2, i4);
                    const auto&[deta34, dphi34, dR34] = thisjet_gen->pairs.get(i3, i4);

                    std::array<res4_entry, 3> dipole_entries;
                    std::array<res4_entry, 3> tee_entries;
                    std::array<res4_entry, 4> triangle_entries;

                    innermost_level<ResultType, JetType::pairType::distances_squared, true>(
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

                            axes_gen,

                            tolerance2,
                            tri_tolerance);

                    if constexpr(doUnmatched){
                        if(!matched4){
                            for(unsigned i=0; i<3; ++i){
                                if(dipole_entries[i].isShape){
                                    unmatched_gen->fill_dipole(
                                            dipole_entries[i].idx_R,
                                            dipole_entries[i].idx_r,
                                            dipole_entries[i].idx_c,
                                            wt);
                                }
                                if(tee_entries[i].isShape){
                                    unmatched_gen->fill_tee(
                                            tee_entries[i].idx_R,
                                            tee_entries[i].idx_r,
                                            tee_entries[i].idx_c,
                                            wt);
                                }
                            }
                            for(unsigned i=0; i<4; ++i){
                                if(triangle_entries[i].isShape){
                                    unmatched_gen->fill_triangle(
                                            triangle_entries[i].idx_R,
                                            triangle_entries[i].idx_r,
                                            triangle_entries[i].idx_c,
                                            wt);
                                }
                            }
                        }
                    }//end if constexpr (doUnmatched)

                    if constexpr(doTransfer){
                        res4_transferloop<TransferResultType, ResultType>(
                            *transfer_ans,
                            *untransfered_reco,
                            *untransfered_gen,
                            dipole_entries,
                            tee_entries,
                            triangle_entries,
                            thisjet_reco,
                            *n1, *n2, *n3, *n4,
                            *axes_reco,
                            tolerance2,
                            tri_tolerance,
                            wt);
                    }//end if constexpr (doTransfer) 
                }//end loop over i4
            }//end loop over i3
        }//end loop over i2
    }//end loop over i1
}//end res4_standalone()

#endif
