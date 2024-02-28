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
        T R1 = rin.floatDRs[i0][i1];
        T R2 = rin.floatDRs[i0][i2];
        T R3 = rin.floatDRs[i0][i3];
        T R4 = rin.floatDRs[i1][i2];
        T R5 = rin.floatDRs[i1][i3];
        T R6 = rin.floatDRs[i2][i3];

        std::array<T, 6> subdRs({{R1, R2, R3, R4, R5, R6}});
        std::sort(subdRs.begin(), subdRs.end());

        T RL = subdRs[0];

        if ((RL - subdRs[1]) < rin.shapetol * RL &&
            (RL - subdRs[2]) < rin.shapetol * RL){

            if ((RL - subdRs[3]) < rin.shapetol * RL){
                T Rsquareedge = RL * invsqrt2;
                if (std::abs(Rsquareedge - subdRs[4]) < rin.shapetol * RL &&
                    std::abs(Rsquareedge - subdRs[5]) < rin.shapetol * RL){

                    shape_idx = 1; //square
                    return;
                }
            } else {
                T Rinnertri = 0.6666666666666666 * RL;

                if (std::abs(Rinnertri - subdRs[3]) < rin.shapetol * RL &&
                    std::abs(Rinnertri - subdRs[4]) < rin.shapetol * RL &&
                    std::abs(Rinnertri - subdRs[5]) < rin.shapetol * RL){

                    shape_idx = 2; //triangle 
                    return;
                }
            }
        }
        shape_idx = 0;
        return;
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

        T RL = rin.floatDRs[qi0][qi1];
        RL_idx = static_cast<unsigned>(rin.coarseRL->index(RL) + 1);

        T R2 = rin.floatDRs[qi0][qi2];
        T R3 = rin.floatDRs[qi0][qi3];
        T R4 = rin.floatDRs[qi1][qi2];
        T R5 = rin.floatDRs[qi1][qi3];
        T R6 = rin.floatDRs[qi2][qi3];

        std::array<std::pair<unsigned, T>, 4> subdRs;
        subdRs[0] = {0,R2};
        subdRs[1] = {1,R3};
        subdRs[2] = {2,R4};
        subdRs[3] = {3,R5}; 
        std::sort(subdRs.begin(), subdRs.end(), [](auto& left, auto& right) {
            return left.second < right.second;
        }); 

        if (subdRs[3].second ==0){
            shape_idx = 0;//no useful shape if any distances are zero
            r_idx = 0;
            ct_idx = 0;

            return;
        } else if((RL - subdRs[0].second)/RL < rin.shapetol &&
                  (RL - subdRs[1].second)/RL < rin.shapetol){

            shape_idx = 3;//equilateral triangle

            T r, ct;
            switch(subdRs[0].first){
                case 0: //the fourth point is qi3, and it's connected via the 102 vertex
                    {r = R3/RL;
                    T eta01 = rin.etas[qi1] - rin.etas[qi0];
                    T phi01 = rin.phis[qi1] - rin.phis[qi0];
                    T eta03 = rin.etas[qi3] - rin.etas[qi0];
                    T phi03 = rin.phis[qi3] - rin.phis[qi0];
                    ct = (eta01*eta03 + phi01*phi03)/(RL*R3);
                    break;}
                case 1: //the fourth point is qi2, and it's connected via the 103 vertex
                    {r = R2/RL;
                    T eta01 = rin.etas[qi1] - rin.etas[qi0];
                    T phi01 = rin.phis[qi1] - rin.phis[qi0];
                    T eta02 = rin.etas[qi2] - rin.etas[qi0];
                    T phi02 = rin.phis[qi2] - rin.phis[qi0];
                    ct = (eta01*eta02 + phi01*phi02)/(RL*R2);
                    break;}
                case 2: //the fourth point is qi3, and it's connected via the 210 vertex
                    {r = R5/RL;
                    T eta10 = rin.etas[qi0] - rin.etas[qi1];
                    T phi10 = rin.phis[qi0] - rin.phis[qi1];
                    T eta13 = rin.etas[qi3] - rin.etas[qi1];
                    T phi13 = rin.phis[qi3] - rin.phis[qi1];
                    ct = (eta10*eta13 + phi10*phi13)/(RL*R5);
                    break;}
                case 3: //the fourth point is q2, and it's connected via the 310 vertex
                    {r = R4/RL;
                    T eta10 = rin.etas[qi0] - rin.etas[qi1];
                    T phi10 = rin.phis[qi0] - rin.phis[qi1];
                    T eta12 = rin.etas[qi2] - rin.etas[qi1];
                    T phi12 = rin.phis[qi2] - rin.phis[qi1];
                    ct = (eta10*eta12 + phi10*phi12)/(RL*R4);
                    break;}
            }
            r_idx = static_cast<unsigned>(rin.r_triangle->index(r)+1);
            ct_idx = static_cast<unsigned>(rin.ct_triangle->index(ct)+1);
            return;
        }

        //if we have gotten this far, then it's either tee or dipole (or nothing)
        T etaA = rin.etas[qi0]; 
        T etaB = rin.etas[qi1];
        T etaC = rin.etas[qi2];
        T etaD = rin.etas[qi3];
        T phiA = rin.phis[qi0];
        T phiB = rin.phis[qi1];
        T phiC = rin.phis[qi2];
        T phiD = rin.phis[qi3];

        if (phiA > phiB){
            std::swap(phiA, phiB);
            std::swap(etaA, etaB);
        }

        if (phiC > phiD){
            std::swap(phiC, phiD);
            std::swap(etaC, etaD);
        }

        //need to handle crossing the +-pi boundary
        bool ABcrosses = phiB - phiA > M_PI;
        bool CDcrosses = phiD - phiD > M_PI;

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

        if(ABcrosses && !CDcrosses){
            //need to check which way to transform AB
            if(phiB < phiC){ 
                //if the rightmost piece of AB is to the left of the leftmost piece of CD
                //then transform AB to the right
                phiB += 2*M_PI; 
                phiA += 2*M_PI; 
            }
        } else if(CDcrosses){
            if(phiD < phiA){
                //if the rightmost piece of CD is to the left of the leftmost piece of AB
                //then transform CD to the right
                phiD += 2*M_PI;
                phiC += 2*M_PI;
            }
        }

        //now we are guarentees to have both lines in the same range

        //bounding box check
        if(phiD < phiA || phiB < phiC){
            shape_idx = 0;
            r_idx = 0;
            ct_idx = 0;
            return;
        } else {
            T minEtaAB = etaA;
            T minEtaCD = etaC;
            T maxEtaAB = etaB;
            T maxEtaCD = etaD;
            if(minEtaAB > maxEtaAB){
                std::swap(minEtaAB, maxEtaAB);
            }
            if(minEtaCD > maxEtaCD){
                std::swap(minEtaCD, maxEtaCD);
            }
            if(minEtaAB > maxEtaCD || minEtaCD > maxEtaAB){
                shape_idx = 0;
                r_idx = 0;
                ct_idx = 0;
                return;
            } else {//pass bounding box check
               
                T ABeta = etaB - etaA;
                T ABphi = phiB - phiA;

                T CDeta = etaD - etaC;
                T CDphi = phiD - phiC;

                T slopecross = ABeta*CDphi - ABphi*CDeta;
                if(slopecross == 0){//colinear
                    shape_idx = 0;
                    r_idx = 0;
                    ct_idx = 0;
                    return;
                } else {
                    T ACeta = etaC - etaA;
                    T ACphi = phiC - phiA;
                    T AC_ABcross = ACeta*ABphi - ACphi*ABeta;

                    T invslopecross = 1/slopecross;
                    T c = AC_ABcross*invslopecross;

                    if (std::abs(c - 0.5) > rin.shapetol){
                        shape_idx = 0;
                        r_idx = 0;
                        ct_idx = 0;
                        return;
                    } else {
                        T AC_CDcross = ACeta*CDphi - ACphi*CDeta;     
                        T a = AC_CDcross*invslopecross;
                        if(std::abs(a - 0.5) < rin.shapetol){
                            shape_idx = 1; //dipole configuration
                            r_idx = static_cast<unsigned>(rin.r_dipole->index(R6/RL) + 1);
                            T ct = (ABeta * CDeta + ABphi * CDphi)/(RL*R6);
                            ct_idx = static_cast<unsigned>(rin.ct_dipole->index(ct) + 1);
                            return;
                        } else if(std::abs(a-0) < rin.shapetol || std::abs(a-1) < rin.shapetol){
                            shape_idx = 2; //tee configuration
                            r_idx = static_cast<unsigned>(rin.r_tee->index(R6/RL) + 1);
                            T ct = (ABeta * CDeta + ABphi * CDphi)/(RL*R6);
                            ct_idx = static_cast<unsigned>(rin.ct_tee->index(ct) + 1);
                            return;
                        } else{
                            shape_idx = 0;
                            r_idx = 0;
                            ct_idx = 0;
                            return;
                        }
                    }
                }
            }
        }
        shape_idx = 0;
        r_idx = 0;
        ct_idx = 0;
        return;
    }
};

#endif
