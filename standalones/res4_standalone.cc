#include "res4_standalone.h"
#include "SRothman/SimonTools/src/deltaR.h"
#include "SRothman/SimonTools/src/histutil.h"

//this is in principle configurable
//but it is important that the target ratios not be 1
constexpr double triangle_ratio_LM = 5.0/4.0;
constexpr double triangle_ratio_MS = 4.0/3.0;

static void compute_c(
        const double eta1, const double phi1,
        const double eta2, const double phi2,
        const double eta3, const double phi3,
        const double eta4, const double phi4,
        const double dR12, const double dR34,
        double& c){

    double deta12 = deltaEta(eta1, eta2);
    double dphi12 = deltaPhi(phi1, phi2);
    double deta34 = deltaEta(eta3, eta4);
    double dphi34 = deltaPhi(phi3, phi4);

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

        const double dR12,
        const double dR34,

        const double wt,

        const double tolerance2,

        const standaloneEEC::axis& ax_R,
        const standaloneEEC::axis& ax_r_dipole,
        const standaloneEEC::axis& ax_c_dipole,
        const standaloneEEC::axis& ax_r_tee,
        const standaloneEEC::axis& ax_c_tee){

    /*
     * Check the configuration (12)(23) for a triangle
     */

    bool isDipole, isTee;

    double mideta12 = (eta1 + eta2)/2.0;
    double midphi12 = angle_midpoint(phi1, phi2);

    double mideta34 = (eta3 + eta4)/2.0;
    double midphi34 = angle_midpoint(phi3, phi4);

    double dR2_midpoints = deltaR2(mideta12, midphi12,
                               mideta34, midphi34);
    isDipole = dR2_midpoints < tolerance2;

    if (dR12 < dR34){
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
        double R,r,c;
        compute_r_R(dR12, dR34, r, R);
        compute_c(eta1, phi1,
                  eta2, phi2,
                  eta3, phi3,
                  eta4, phi4,
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

        const double eta1, const double phi1,
        const double eta2, const double phi2,
        const double eta3, const double phi3,
        const double eta4, const double phi4,

        const double dR12, const double dR13,
        const double dR14, const double dR23,
        const double dR24, const double dR34,

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

    double etaA, phiA;
    double etaB, phiB;
    double etaC, phiC;
    double RAB, RAC, RBC;
    double RA4; //distance of 4th particle from right angle

    if(dR12 > dR13){
        if(dR13 > dR23){
            // (12) > (13) > (23)
            // so A = 3, B = 2, C = 1
            etaA = eta3;
            phiA = phi3;
            
            etaB = eta2;
            phiB = phi2;
            
            etaC = eta1;
            phiC = phi1;

            RAB = dR23;
            RAC = dR13;
            RBC = dR12;

            RA4 = dR34;
            
        } else if (dR12 > dR23){
            // (12) > (23) > (13)
            // so A = 3, B = 1, C = 2
            etaA = eta3;
            phiA = phi3;

            etaB = eta1;
            phiB = phi1;

            etaC = eta2;
            phiC = phi2;

            RAB = dR13;
            RAC = dR23;
            RBC = dR12;

            RA4 = dR34;

        } else {
            // (23) > (12) > (13)
            // so A = 1, B = 3, C = 2
            etaA = eta1;
            phiA = phi1;

            etaB = eta3;
            phiB = phi3;

            etaC = eta2;
            phiC = phi2;

            RAB = dR13;
            RAC = dR12;
            RBC = dR23;

            RA4 = dR14;

        }
    } else{
        if(dR12 > dR23){ 
            // (13) > (12) > (23)
            // so A = 2, B = 3, C = 1
            etaA = eta2;
            phiA = phi2;

            etaB = eta3;
            phiB = phi3;

            etaC = eta1;
            phiC = phi1;

            RAB = dR23;
            RAC = dR12;
            RBC = dR13;

            RA4 = dR24;

        } else if(dR13 > dR23){ 
            // (13) > (23) > (12)
            // so A = 2, B = 1, C = 3
            etaA = eta2;
            phiA = phi2;

            etaB = eta1;
            phiB = phi1;

            etaC = eta3;
            phiC = phi3;

            RAB = dR12;
            RAC = dR23;
            RBC = dR13;

            RA4 = dR24;

        } else {
            // (23) > (13) > (12)
            // so A = 1, B = 2, C = 3
            etaA = eta1;
            phiA = phi1;

            etaB = eta2;
            phiB = phi2;

            etaC = eta3;
            phiC = phi3;

            RAB = dR12;
            RAC = dR13;
            RBC = dR23;

            RA4 = dR14;
        }
    }

    /*
     * Now we have the order of the triangle sides
     * we can compute the angles and distances
     */
    bool passLM = std::abs(RBC/RAC - triangle_ratio_LM) < tolerance2;
    bool passMS = std::abs(RAC/RAB - triangle_ratio_MS) < tolerance2;
    if(passLM && passMS){
        double R, r, c;

        R = RAC;
        r = RA4/RAC;

        double detaA4 = deltaEta(etaA, eta4);
        double dphiA4 = deltaPhi(phiA, phi4);

        double detaAC = deltaEta(etaA, etaC);
        double dphiAC = deltaPhi(phiA, phiC);

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

        double detaAB = deltaEta(etaA, etaB);
        double dphiAB = deltaPhi(phiA, phiB);

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
            double E2 = thejet.Es[i2];
            double eta2 = thejet.etas[i2];
            double phi2 = thejet.phis[i2];
            double dR12 = deltaR(eta1, phi1, eta2, phi2);

            for(unsigned i3=i2+1; i3<thejet.N; ++i3){
                double E3 = thejet.Es[i3];
                double eta3 = thejet.etas[i3];
                double phi3 = thejet.phis[i3];
                double dR13 = deltaR(eta1, phi1, eta3, phi3);
                double dR23 = deltaR(eta2, phi2, eta3, phi3);

                for(unsigned i4=i3+1; i4<thejet.N; ++i4){
                    double E4 = thejet.Es[i4];
                    double eta4 = thejet.etas[i4];
                    double phi4 = thejet.phis[i4];
                    double dR14 = deltaR(eta1, phi1, eta4, phi4);
                    double dR24 = deltaR(eta2, phi2, eta4, phi4);
                    double dR34 = deltaR(eta3, phi3, eta4, phi4);

                    double wt = E1*E2*E3*E4;

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
                            dR12, dR34,
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
                            dR13, dR24,
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
                            dR14, dR23,
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
                            eta1, phi1,
                            eta2, phi2,
                            eta3, phi3,
                            eta4, phi4,

                            dR12, dR13, 
                            dR14, dR23,
                            dR24, dR34,

                            wt,

                            tolerance2,

                            ax_R,
                            ax_r_triangle,
                            ax_c_triangle);

                    //(124)(3)
                    check_triangle(
                            ans,
                            eta1, phi1,
                            eta2, phi2,
                            eta4, phi4,
                            eta3, phi3,

                            dR12, dR14, 
                            dR13, dR24,
                            dR23, dR34,

                            wt,

                            tolerance2,

                            ax_R,
                            ax_r_triangle,
                            ax_c_triangle);

                    //(143)(2)
                    check_triangle(
                            ans,
                            eta1, phi1,
                            eta4, phi4,
                            eta3, phi3,
                            eta2, phi2,

                            dR14, dR13, 
                            dR12, dR34,
                            dR24, dR23,

                            wt,

                            tolerance2,

                            ax_R,
                            ax_r_triangle,
                            ax_c_triangle);

                    //(423)(1)
                    check_triangle(
                            ans,
                            eta4, phi4,
                            eta2, phi2,
                            eta3, phi3,
                            eta1, phi1,

                            dR24, dR34, 
                            dR14, dR23,
                            dR12, dR13,

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
