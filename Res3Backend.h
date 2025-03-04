#ifndef SROTHMAN_EEC_RES3_BACKEND_H
#define SROTHMAN_EEC_RES3_BACKEND_H

#include "Adjacency.h"

#include "SRothman/SimonTools/src/deltaR.h"
#include "SRothman/SimonTools/src/histutil.h"
#include "SRothman/SimonTools/src/util.h"

#include <memory>

#include <array>

template <bool distances_squared>
inline void get_indices_CA(const double dR12, const double dR13, const double dR23,
                           const double E1, const double eta1, const double phi1,
                           const double E2, const double eta2, const double phi2,
                           const double E3, const double eta3, const double phi3,
                           const double deta12, const double dphi12,
                           const double deta13, const double dphi13,
                           const double deta23, const double dphi23,
                           unsigned& R_idx, unsigned& r_idx, unsigned& c_idx,
                           const EEC::Res3Axes& axes) noexcept {

    std::array<std::pair<double, unsigned>, 3> dRs = {
        {std::make_pair(dR12, 0), std::make_pair(dR13, 1), std::make_pair(dR23, 2)}
    };
    if constexpr(distances_squared){
        for(auto& dR : dRs){
            dR.first = std::sqrt(dR.first);
        }
    }
    std::sort(dRs.begin(), dRs.end(),
              [](const std::pair<double, unsigned>& a, const std::pair<double, unsigned>& b){
                return a.first > b.first;
            }
    );

    if (dRs[0].first == 0){
        R_idx = 0;
        r_idx = 0;
        c_idx = 0;
    } else if (dRs[2].first == 0){
        R_idx = simon::getIndex(dRs[0].first, axes.R);
        r_idx = 0;
        c_idx = 0;
    } else {
        double etaFar, phiFar;
        double etaMid, phiMid;
        double detaSmall, dphiSmall;
        if (dRs[2].second == 0){//shortest leg is dR12
            etaFar = eta3;
            phiFar = phi3;

            double Esum = E1 + E2;
            etaMid = (E1*eta1 + E2*eta2) / Esum;
            phiMid = (E1*phi1 + E2*phi2) / Esum;

            detaSmall = deta12;
            dphiSmall = dphi12;
        } else if(dRs[2].second == 1){//shortest leg is dR13
            etaFar = eta2;
            phiFar = phi2;

            double Esum = E1 + E3;
            etaMid = (E1*eta1 + E3*eta3) / Esum;
            phiMid = (E1*phi1 + E3*phi3) / Esum;

            detaSmall = deta13;
            dphiSmall = dphi13;
        } else {//shortest leg is dR23
            etaFar = eta1;
            phiFar = phi1;

            double Esum = E2 + E3;
            etaMid = (E2*eta2 + E3*eta3) / Esum;
            phiMid = (E2*phi2 + E3*phi3) / Esum;

            detaSmall = deta23;
            dphiSmall = dphi23;
        }

        double detaLong = simon::deltaEta(etaFar, etaMid);
        double dphiLong = simon::deltaPhi(phiFar, phiMid);

        double R = std::sqrt(simon::square(detaLong) + simon::square(dphiLong));
        double r = dRs[2].first / R;

        double dot = (detaLong*detaSmall + dphiLong*dphiSmall) / (R * dRs[2].first);
        if (dot > 1){
            dot = 1;
        } else if (dot < -1){
            dot = -1;
        }
        double c = std::acos(dot);
        if (c > M_PI/2){
            c = M_PI - c;
        }

        R_idx = simon::getIndex(R, axes.R);
        r_idx = simon::getIndex(r, axes.r);
        c_idx = simon::getIndex(c, axes.c);
    }
}

template <bool distances_squared>
inline void get_indices(const double dR12, const double dR13, const double dR23,
                        unsigned& R_idx, unsigned& r_idx, unsigned& c_idx,
                        const EEC::Res3Axes& axes) noexcept {

    std::array<double, 3> dRs = {dR12, dR13, dR23};
    if constexpr(distances_squared){
        for(auto& dR : dRs){
            dR = std::sqrt(dR);
        }
    }
    std::sort(dRs.begin(), dRs.end(),
             [](const double a, const double b){return a > b;});
    if (dRs[0] == 0){
        R_idx = 0;
        r_idx = 0;
        c_idx = 0;
    } else if (dRs[2] == 0){
        R_idx = simon::getIndex(dRs[0], axes.R);
        r_idx = 0;
        c_idx = 0;
    } else {
        double R = dRs[0];
        double r = dRs[2]/dRs[1];
        double c = std::asin(std::sqrt(1 - simon::square((dRs[0] - dRs[1])/dRs[2])));
        R_idx = simon::getIndex(R, axes.R);
        r_idx = simon::getIndex(r, axes.r);
        c_idx = simon::getIndex(c, axes.c);
    }
}

template <bool isCA, class JetType, class TransferContainer>
inline void res3_transferloop(
        EEC::Res3TransferResult<TransferContainer>& transfer,

        const std::shared_ptr<const JetType> thisjet_reco,

        const EEC::Res3Axes& axes_reco,

        const unsigned R_idx_gen,
        const unsigned r_idx_gen,
        const unsigned c_idx_gen,

        const EEC::neighborhood& n1,
        const EEC::neighborhood& n2,
        const EEC::neighborhood& n3,

        const double wt_gen) noexcept {

    for(const EEC::neighbor& j1 : n1){
        [[maybe_unused]] auto&[E1, eta1, phi1] = thisjet_reco->singles.get(j1.idx);
        const double twt1 = wt_gen * j1.wt;

        for(const EEC::neighbor& j2 : n2){
            [[maybe_unused]] auto&[E2, eta2, phi2] = thisjet_reco->singles.get(j2.idx);
            const double twt2 = twt1 * j2.wt;

            const auto&[deta12, dphi12, dR12] = thisjet_reco->pairs.get(j1.idx, j2.idx);

            for(const EEC::neighbor& j3 : n3){
                [[maybe_unused]] auto&[E3, eta3, phi3] = thisjet_reco->singles.get(j3.idx);
                const double twt3 = twt2 * j3.wt;

                const auto&[deta13, dphi13, dR13] = thisjet_reco->pairs.get(j1.idx, j3.idx);
                const auto&[deta23, dphi23, dR23] = thisjet_reco->pairs.get(j2.idx, j3.idx);

                unsigned R_idx_reco, r_idx_reco, c_idx_reco;
                if constexpr (isCA){
                    get_indices_CA<JetType::pairType::distances_squared>(
                        dR12, dR13, dR23,
                        E1, eta1, phi1,
                        E2, eta2, phi2,
                        E3, eta3, phi3,
                        deta12, dphi12,
                        deta13, dphi13,
                        deta23, dphi23,
                        R_idx_reco, r_idx_reco, c_idx_reco,
                        axes_reco);
                } else {
                    get_indices<JetType::pairType::distances_squared>(
                            dR12, dR13, dR23, 
                            R_idx_reco, r_idx_reco, c_idx_reco, 
                            axes_reco);
                }
#ifdef CHECK_BY_HAND
                printf("(%g, %g), (%g, %g), (%g, %g): [transfer]\n",
                        eta1, phi1, eta2, phi2, eta3, phi3);
                printf("\tE1*E2*E3 = %g\n", E1*E2*E3);
                printf("\tR_idx = %u\n", R_idx_reco);
                printf("\tr_idx = %u\n", r_idx_reco);

                printf("\tE1 = %g, E2 = %g, E3 = %g\n", E1, E2, E3);
#endif

                transfer.fill(R_idx_reco, r_idx_reco, c_idx_reco,
                              R_idx_gen, r_idx_gen, c_idx_gen, 
                              twt3, wt_gen);
            }
        }
    }
}

template <bool isCA, class BasicContainer, class JetType, bool doUnmatched, class TransferContainer, bool doTransfer>
inline void res3_mainloop(
        EEC::Res3Result<BasicContainer>& result,
        [[maybe_unused]] EEC::Res3Result<BasicContainer>* unmatched_gen,
        [[maybe_unused]] EEC::Res3TransferResult<TransferContainer>* transfer,

        [[maybe_unused]] const std::shared_ptr<const JetType> thisjet_reco,
        const std::shared_ptr<const JetType> thisjet_gen,
        [[maybe_unused]] const std::vector<bool>* const matched,

        [[maybe_unused]] const std::shared_ptr<const EEC::Adjacency> adj,

        [[maybe_unused]] const EEC::Res3Axes* const axes_reco,
        const EEC::Res3Axes& axes_gen) noexcept {//TODO: REMOVE [[maybe_unused]] here

    for(unsigned i1=0; i1<thisjet_gen->N; ++i1){
        const auto&[E1, eta1, phi1] = thisjet_gen->singles.get(i1);
        [[maybe_unused]] bool matched1;
        if constexpr(doUnmatched){
            matched1 = matched->at(i1);
        }
        EEC::neighborhood const * n1 = nullptr;
        if constexpr(doTransfer){
            n1 = &(adj->get_neighborhood(i1));
        }

        /*
         * One-particle part
         * This has symmetry factor 1
         * And coordinates (R, r, c) = (0, 0, 0)
         */
        double wt = E1*E1*E1;
#ifdef CHECK_BY_HAND
        printf("(%g, %g): [1-particle]\n", eta1, phi1);
        printf("\tE1^3 = %g\n", E1*E1*E1);
        printf("\tR_idx = 0\n");
        printf("\tr_idx = 0\n");
        printf("\tc_idx = 0\n");
        if constexpr(doUnmatched){
            printf("\tmatched = %d\n", matched1);
        }
#endif
        result.fill(0, 0, 0, wt);
        if constexpr(doUnmatched){
            if(!matched1){
                unmatched_gen->fill(0, 0, 0, wt);
            }
        }
        if constexpr(doTransfer){
            res3_transferloop<isCA>(
                    *transfer,
                    thisjet_reco,
                    *axes_reco,
                    0, 0, 0,
                    *n1,
                    *n1,
                    *n1,
                    wt);
        }


        for(unsigned i2=i1+1; i2<thisjet_gen->N; ++i2){
            const auto&[E2, eta2, phi2] = thisjet_gen->singles.get(i2);
            [[maybe_unused]] bool matched2;
            if constexpr(doUnmatched){
                matched2 = matched1 && matched->at(i2);
            }
            EEC::neighborhood const * n2 = nullptr;
            if constexpr(doTransfer){
                n2 = &(adj->get_neighborhood(i2));
            }

            const double E12 = E1*E2;

            const auto&[deta12, dphi12, dR12] = thisjet_gen->pairs.get(i1, i2);

            /*
             * Two-particle part
             * This has symmetry factor 3!/2! = 3
             * And RL = R12, r=c=0
             */
            unsigned R_idx;
            if constexpr(JetType::pairType::distances_squared){
                R_idx = simon::getIndex(std::sqrt(dR12), axes_gen.R);
            } else {
                R_idx = simon::getIndex(dR12, axes_gen.R);
            }
            wt = 3 * E12 * (E1 + E2);
            result.fill(R_idx, 0, 0, wt);
#ifdef CHECK_BY_HAND
            printf("(%g, %g), (%g, %g): [2-particle]\n", eta1, phi1, eta2, phi2);
            printf("\tE1*E2*(E1+E2) = %g\n", E12*(E1+E2));
            printf("\tR_idx = %u\n", R_idx);
            printf("\tr_idx = 0\n");
            printf("\tc_idx = 0\n");
            if constexpr(doUnmatched){
                printf("\tmatched = %d\n", matched2);
            }
#endif
            if constexpr(doUnmatched){
                if(!matched2){
                    unmatched_gen->fill(R_idx, 0, 0, wt);
                }
            }
            if constexpr(doTransfer){
                double wt1 = 3*E12*E1;
                res3_transferloop<isCA>(
                        *transfer,
                        thisjet_reco,
                        *axes_reco,
                        R_idx, 0, 0,
                        *n1,
                        *n1,
                        *n2,
                        wt1);
                double wt2 = 3*E12*E2;
                res3_transferloop<isCA>(
                        *transfer,
                        thisjet_reco,
                        *axes_reco,
                        R_idx, 0, 0,
                        *n1,
                        *n2,
                        *n2,
                        wt2);
            }

            for(unsigned i3=i2+1; i3<thisjet_gen->N; ++i3){
                const auto&[E3, eta3, phi3] = thisjet_gen->singles.get(i3);
                [[maybe_unused]] bool matched3;
                if constexpr(doUnmatched){
                    matched3 = matched2 && matched->at(i3);
                }
                EEC::neighborhood const * n3 = nullptr;
                if constexpr(doTransfer){
                    n3 = &(adj->get_neighborhood(i3));
                }

                const double E123 = E12*E3;

                const auto&[deta13, dphi13, dR13] = thisjet_gen->pairs.get(i1, i3);
                const auto&[deta23, dphi23, dR23] = thisjet_gen->pairs.get(i2, i3);

                /*
                 * Three-particle part
                 * This has symmetry factor 3! = 6
                 * and full coordinates
                 */
                unsigned r_idx, c_idx;
                if constexpr (isCA){
                    get_indices_CA<JetType::pairType::distances_squared>(
                        dR12, dR13, dR23,
                        E1, eta1, phi1,
                        E2, eta2, phi2,
                        E3, eta3, phi3,
                        deta12, dphi12,
                        deta13, dphi13,
                        deta23, dphi23,
                        R_idx, r_idx, c_idx,
                        axes_gen);
                } else {
                    get_indices<JetType::pairType::distances_squared>(
                            dR12, dR13, dR23, 
                            R_idx, r_idx, c_idx, 
                            axes_gen);
                }
#ifdef CHECK_BY_HAND
                printf("(%g, %g), (%g, %g), (%g, %g): [3-particle]\n",
                        eta1, phi1, eta2, phi2, eta3, phi3);
                printf("\tE1*E2*E3 = %g\n", E1*E2*E3);
                printf("\tR_idx = %u\n", R_idx);
                printf("\tr_idx = %u\n", r_idx);
                printf("\tc_idx = %u\n", c_idx);
                if constexpr(doUnmatched){
                    printf("\tmatched = %d\n", matched3);
                }
#endif
                
                wt = E123*6;
                result.fill(R_idx, r_idx, c_idx, wt);
                if constexpr(doUnmatched){
                    if(!matched3){
                        unmatched_gen->fill(R_idx, r_idx, c_idx, wt);
                    }
                }
                if constexpr(doTransfer){
                    res3_transferloop<isCA>(
                            *transfer,
                            thisjet_reco,
                            *axes_reco,
                            R_idx, r_idx, c_idx,
                            *n1,
                            *n2,
                            *n3,
                            wt);
                }
            }
        }
    }
}

#endif
