#include "res4_standalone.h"
#include "SRothman/SimonTools/src/deltaR.h"
#include "SRothman/SimonTools/src/histutil.h"
#include "SRothman/SimonTools/src/util.h"

//this is in principle configurable
//but it is important that the target ratios not be 1
constexpr double triangle_ratio_LM = 5.0/4.0;
constexpr double triangle_ratio_MS = 4.0/3.0;

static void compute_c(
        const double deta12, const double dphi12,
        const double deta34, const double dphi34,
        const double dR12, const double dR34,
        double& c){

    double dot = deta12*deta34 + dphi12*dphi34;
    double ct = dot/(dR12*dR34);

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

static void compute_r_R(
        const double dR12,
        const double dR34,
        double& r,
        double& R){

    if(dR12 > dR34){
        R = dR12;
        r = dR34/dR12;
    } else {
        R = dR34;
        r = dR12/dR34;
    }
}

static double angle_midpoint(
        const double phi1, 
        const double phi2){
    /*
     * Compute the midpoint of two polar angles,
     * wrapping around as necessary
     * This is useful for the dipole and tee axes
     */
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

template <class Container>
static void check_tee_dipole(
        standaloneEEC::res4_result<Container>& ans,

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
        const standaloneEEC::axis& ax_c_tee){

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
        double dR12 = std::sqrt(dR12_2);
        double dR34 = std::sqrt(dR34_2);

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
            ans.fill_dipole(idx_R, idx_r, idx_c, wt);
        }
        if (isTee){
            unsigned idx_r = getIndex(r, ax_r_tee);
            unsigned idx_c = getIndex(c, ax_c_tee);
            ans.fill_tee(idx_R, idx_r, idx_c, wt);
        }
    }
}

template <class Container>
static void check_triangle(
        standaloneEEC::res4_result<Container>& ans,

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

        const double tolerance2,

        const standaloneEEC::axis& ax_R,
        const standaloneEEC::axis& ax_r_triangle,
        const standaloneEEC::axis& ax_c_triangle){

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
    bool passLM = std::abs(std::sqrt(RBC_2/RAC_2) - triangle_ratio_LM) < tolerance2;
    bool passMS = std::abs(std::sqrt(RAC_2/RAB_2) - triangle_ratio_MS) < tolerance2;

    if(passLM && passMS){
        double R, r, c;

        double RAC = std::sqrt(RAC_2);
        double RA4 = std::sqrt(RA4_2);

        R = RAC;
        r = RA4/RAC;

        double dotA4AC = detaA4*detaAC + dphiA4*dphiAC;

        double cosc = dotA4AC/(RA4*RAC);
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
        ans.fill_triangle(idx_R, idx_r, idx_c, wt);
    }
}

template <class Container>
static void res4_mainloop(
        standaloneEEC::res4_result<Container>& ans,

        const standaloneEEC::EECjet& thejet,
        const standaloneEEC::axis& ax_R,
        const standaloneEEC::axis& ax_r_dipole,
        const standaloneEEC::axis& ax_c_dipole,
        const standaloneEEC::axis& ax_r_tee,
        const standaloneEEC::axis& ax_c_tee,
        const standaloneEEC::axis& ax_r_triangle,
        const standaloneEEC::axis& ax_c_triangle,

        double tolerance2){
    /*
     * This is the main loop over four particles
     * NB in the specific case of res4 we are only interested
     * in the case where the four particles are distinct
     *
     * This allows us to save some loops
     * This also means we don't have to worry about 
     * symmetry factors
     */

    for (unsigned i1=0; i1<thejet.N; ++i1){
        double E1 = thejet.Es[i1];
        double eta1 = thejet.etas[i1];
        double phi1 = thejet.phis[i1];

        for (unsigned i2=i1+1; i2<thejet.N; ++i2){
            double E12 = E1 * thejet.Es[i2];
            double eta2 = thejet.etas[i2];
            double phi2 = thejet.phis[i2];

            double deta12 = deltaEta(eta1, eta2);
            double dphi12 = deltaPhi(phi1, phi2);
            double dR12_2 = square(deta12) + square(dphi12);

            for(unsigned i3=i2+1; i3<thejet.N; ++i3){
                double E123 = E12 * thejet.Es[i3];
                double eta3 = thejet.etas[i3];
                double phi3 = thejet.phis[i3];

                double deta13 = deltaEta(eta1, eta3);
                double dphi13 = deltaPhi(phi1, phi3);
                double dR13_2 = square(deta13) + square(dphi13);

                double deta23 = deltaEta(eta2, eta3);
                double dphi23 = deltaPhi(phi2, phi3);
                double dR23_2 = square(deta23) + square(dphi23);

                for(unsigned i4=i3+1; i4<thejet.N; ++i4){
                    double wt = E123 * thejet.Es[i4];
                    double eta4 = thejet.etas[i4];
                    double phi4 = thejet.phis[i4];

                    double deta14 = deltaEta(eta1, eta4);
                    double dphi14 = deltaPhi(phi1, phi4);
                    double dR14_2 = square(deta14) 
                        + square(dphi14);

                    double deta24 = deltaEta(eta2, eta4);
                    double dphi24 = deltaPhi(phi2, phi4);
                    double dR24_2 = square(deta24) 
                        + square(dphi24);

                    double deta34 = deltaEta(eta3, eta4);
                    double dphi34 = deltaPhi(phi3, phi4);
                    double dR34_2 = square(deta34) 
                        + square(dphi34);

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
                    check_tee_dipole(
                            ans,
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
                    check_tee_dipole(
                            ans,
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
                    check_tee_dipole(
                            ans,
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
                    check_triangle(
                            ans,

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

                            tolerance2,

                            ax_R,
                            ax_r_triangle,
                            ax_c_triangle);

                    //(124)(3)
                    check_triangle(
                            ans,

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

                            tolerance2,

                            ax_R,
                            ax_r_triangle,
                            ax_c_triangle);

                    //(143)(2)
                    check_triangle(
                            ans,

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

                            tolerance2,

                            ax_R,
                            ax_r_triangle,
                            ax_c_triangle);

                    //(423)(1)
                    check_triangle(
                            ans,

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

                            tolerance2,

                            ax_R,
                            ax_r_triangle,
                            ax_c_triangle);

                }//end loop over i4
            }//end loop over i3
        }//end loop over i2
    }//end loop over i1
}//end res4_standalone()
 
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

        double tolerance){

    standaloneEEC::EECjet thejet(J, nt);

    res4_mainloop<res4_multi_array_container>(ans, thejet, 
            ax_R, 
            ax_r_dipole, ax_c_dipole, 
            ax_r_tee, ax_c_tee, 
            ax_r_triangle, ax_c_triangle,
            tolerance*tolerance);
}
 
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

        double tolerance){

    standaloneEEC::EECjet thejet(J, nt);

    res4_mainloop<res4_vector_container>(ans, thejet, 
            ax_R, 
            ax_r_dipole, ax_c_dipole, 
            ax_r_tee, ax_c_tee, 
            ax_r_triangle, ax_c_triangle,
            tolerance*tolerance);
}
