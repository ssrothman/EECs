#ifndef SROTHMAN_EEC_RES3_BACKEND_H
#define SROTHMAN_EEC_RES3_BACKEND_H

#include "Adjacency.h"

#include "SRothman/SimonTools/src/deltaR.h"
#include "SRothman/SimonTools/src/histutil.h"
#include "SRothman/SimonTools/src/util.h"

#include <memory>

#include <array>

template <class ResultType, bool distances_squared>
inline void get_indices_CA(const double dR12, const double dR13, const double dR23,
                           const double E1, const double eta1, const double phi1,
                           const double E2, const double eta2, const double phi2,
                           const double E3, const double eta3, const double phi3,
                           const double deta12, const double dphi12,
                           const double deta13, const double dphi13,
                           const double deta23, const double dphi23,
                           typename ResultType::T& R_idx, 
                           typename ResultType::T& r_idx, 
                           typename ResultType::T& c_idx,
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
        if constexpr(ResultType::SHOULD_BIN){
            R_idx = simon::getIndex(dRs[0].first, axes.R);
        } else {
            R_idx = dRs[0].first;
        }
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

        if constexpr(ResultType::SHOULD_BIN){
            R_idx = simon::getIndex(R, axes.R);
            r_idx = simon::getIndex(r, axes.r);
            c_idx = simon::getIndex(c, axes.c);
        } else {
            R_idx = R;
            r_idx = r;
            c_idx = c;
        }
    }
}

template <class ResultType, bool distances_squared>
inline void get_indices(const double dR12, const double dR13, const double dR23,
                        typename ResultType::T& R_idx,
                        typename ResultType::T& r_idx,
                        typename ResultType::T& c_idx,
                        const EEC::Res3Axes& axes) noexcept {

    std::array<double, 3> dRs = {{dR12, dR13, dR23}};
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
        if constexpr(ResultType::SHOULD_BIN){
            R_idx = simon::getIndex(dRs[0], axes.R);
        } else {
            R_idx = dRs[0];
        }
        r_idx = 0;
        c_idx = 0;
    } else {
        double R = dRs[0];
        double r = dRs[2]/dRs[1];
        double c = std::asin(std::sqrt(1 - simon::square((dRs[0] - dRs[1])/dRs[2])));
        if constexpr(ResultType::SHOULD_BIN){
            R_idx = simon::getIndex(R, axes.R);
            r_idx = simon::getIndex(r, axes.r);
            c_idx = simon::getIndex(c, axes.c);
        } else {
            R_idx = R;
            r_idx = r;
            c_idx = c;
        }
    }
}

template <bool isCA, class JetType, class TransferResultType>
inline void res3_transferloop(
        TransferResultType& transfer,

        const std::shared_ptr<const JetType> thisjet_reco,

        const EEC::Res3Axes& axes_reco,

        const typename TransferResultType::T R_idx_gen,
        const typename TransferResultType::T r_idx_gen,
        const typename TransferResultType::T c_idx_gen,

        const EEC::neighborhood& n1,
        const EEC::neighborhood& n2,
        const EEC::neighborhood& n3,

        const double wt_gen) noexcept {

    for(const EEC::neighbor& j1 : n1){
        const auto& p1 = thisjet_reco->singles.get(j1.idx);
        const double twt1 = wt_gen * j1.wt;

        for(const EEC::neighbor& j2 : n2){
            const auto& p2 = thisjet_reco->singles.get(j2.idx);
            const double twt2 = twt1 * j2.wt;

            const auto& pair12 = thisjet_reco->pairs.get(j1.idx, j2.idx);

            for(const EEC::neighbor& j3 : n3){
                const auto& p3 = thisjet_reco->singles.get(j3.idx);
                const double twt3 = twt2 * j3.wt;

                const auto& pair13 = thisjet_reco->pairs.get(j1.idx, j3.idx);
                const auto& pair23 = thisjet_reco->pairs.get(j2.idx, j3.idx);

                typename TransferResultType::T R_idx_reco, r_idx_reco, c_idx_reco;
                if constexpr (isCA){
                    get_indices_CA<TransferResultType, JetType::pairType::distances_squared>(
                        pair12.floatDR, pair13.floatDR, pair23.floatDR,
                        p1.pt, p1.eta, p1.phi,
                        p2.pt, p2.eta, p2.phi,
                        p3.pt, p3.eta, p3.phi,
                        pair12.deta, pair12.dphi,
                        pair13.deta, pair13.dphi,
                        pair23.deta, pair23.dphi,
                        R_idx_reco, r_idx_reco, c_idx_reco,
                        axes_reco);
                } else {
                    get_indices<TransferResultType, JetType::pairType::distances_squared>(
                            pair12.floatDR, pair13.floatDR, pair23.floatDR, 
                            R_idx_reco, r_idx_reco, c_idx_reco, 
                            axes_reco);
                }
#ifdef CHECK_BY_HAND
                printf("(%g, %g), (%g, %g), (%g, %g): [transfer]\n",
                        p1.eta, p1.phi, p2.eta, p2.phi, p3.eta, p3.phi);
                printf("\tE1*E2*E3 = %g\n", p1.pt*p2.pt*p3.pt);
                std::cout << "\tR_idx = " << R_idx_reco << std::endl;
                std::cout << "\tc_idx = " << c_idx_reco << std::endl;

                printf("\tE1 = %g, E2 = %g, E3 = %g\n", p1.pt, p2.pt, p3.pt);
#endif

                transfer.fill(R_idx_reco, r_idx_reco, c_idx_reco,
                              R_idx_gen, r_idx_gen, c_idx_gen, 
                              twt3, wt_gen);
            }
        }
    }
}

template <bool isCA, class ResultType, class JetType, bool doUnmatched, class TransferResultType, bool doTransfer>
inline void res3_mainloop(
        ResultType& result,
        [[maybe_unused]] ResultType* unmatched_gen,
        [[maybe_unused]] TransferResultType* transfer,

        [[maybe_unused]] const std::shared_ptr<const JetType> thisjet_reco,
        const std::shared_ptr<const JetType> thisjet_gen,
        [[maybe_unused]] const std::vector<bool>* const matched,

        [[maybe_unused]] const std::shared_ptr<const EEC::Adjacency> adj,

        [[maybe_unused]] const EEC::Res3Axes* const axes_reco,
        const EEC::Res3Axes& axes_gen) noexcept {

    for(unsigned i1=0; i1<thisjet_gen->N; ++i1){
        const auto& p1 = thisjet_gen->singles.get(i1);
        [[maybe_unused]] bool matched1;
        if constexpr(doUnmatched){
            matched1 = matched->at(i1);
        }
        [[maybe_unused]] EEC::neighborhood const * n1 = nullptr;
        if constexpr(doTransfer){
            n1 = &(adj->get_neighborhood(i1));
        }

        /*
         * One-particle part
         * This has symmetry factor 1
         * And coordinates (R, r, c) = (0, 0, 0)
         */
        double wt = p1.pt*p1.pt*p1.pt;
#ifdef CHECK_BY_HAND
        printf("(%g, %g): [1-particle]\n", p1.eta, p1.phi);
        printf("\tE1^3 = %g\n", p1.pt*p1.pt*p1.pt);
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
            const auto& p2 = thisjet_gen->singles.get(i2);
            [[maybe_unused]] bool matched2;
            if constexpr(doUnmatched){
                matched2 = matched1 && matched->at(i2);
            }
            [[maybe_unused]] EEC::neighborhood const * n2 = nullptr;
            if constexpr(doTransfer){
                n2 = &(adj->get_neighborhood(i2));
            }

            const double E12 = p1.pt*p2.pt;

            const auto& pair12 = thisjet_gen->pairs.get(i1, i2);

            /*
             * Two-particle part
             * This has symmetry factor 3!/2! = 3
             * And RL = R12, r=c=0
             */
            typename ResultType::T R_idx;
            if constexpr(JetType::pairType::distances_squared){
                R_idx = simon::getIndex(std::sqrt(pair12.floatDR), axes_gen.R);
            } else {
                R_idx = simon::getIndex(pair12.floatDR, axes_gen.R);
            }
            wt = 3 * E12 * (p1.pt + p2.pt);


            result.fill(R_idx, 0, 0, wt);
#ifdef CHECK_BY_HAND
            printf("(%g, %g), (%g, %g): [2-particle]\n", p1.eta, p1.phi, p2.eta, p2.phi);
            printf("\tE1*E2*(E1+E2) = %g\n", E12*(p1.pt+p2.pt));
            std::cout << "\tR = " << R_idx << std::endl;
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
                double wt1 = 3*E12*p1.pt;
                res3_transferloop<isCA>(
                        *transfer,
                        thisjet_reco,
                        *axes_reco,
                        R_idx, 0, 0,
                        *n1,
                        *n1,
                        *n2,
                        wt1);
                double wt2 = 3*E12*p2.pt;
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
                const auto& p3 = thisjet_gen->singles.get(i3);
                [[maybe_unused]] bool matched3;
                if constexpr(doUnmatched){
                    matched3 = matched2 && matched->at(i3);
                }
                [[maybe_unused]] EEC::neighborhood const * n3 = nullptr;
                if constexpr(doTransfer){
                    n3 = &(adj->get_neighborhood(i3));
                }

                const double E123 = E12*p3.pt;

                const auto& pair13 = thisjet_gen->pairs.get(i1, i3);
                const auto& pair23 = thisjet_gen->pairs.get(i2, i3);

                /*
                 * Three-particle part
                 * This has symmetry factor 3! = 6
                 * and full coordinates
                 */
                typename ResultType::T r_idx, c_idx;
                if constexpr (isCA){
                    get_indices_CA<ResultType, JetType::pairType::distances_squared>(
                        pair12.floatDR, pair13.floatDR, pair23.floatDR,
                        p1.pt, p1.eta, p1.phi,
                        p2.pt, p2.eta, p2.phi,
                        p3.pt, p3.eta, p3.phi,
                        pair12.deta, pair12.dphi,
                        pair13.deta, pair13.dphi,
                        pair23.deta, pair23.dphi,
                        R_idx, r_idx, c_idx,
                        axes_gen);
                } else {
                    get_indices<ResultType, JetType::pairType::distances_squared>(
                            pair12.floatDR, pair13.floatDR, pair23.floatDR, 
                            R_idx, r_idx, c_idx, 
                            axes_gen);
                }
#ifdef CHECK_BY_HAND
                printf("(%g, %g), (%g, %g), (%g, %g): [3-particle]\n",
                        p1.eta, p1.phi, p2.eta, p2.phi, p3.eta, p3.phi);
                printf("\tE1*E2*E3 = %g\n", p1.pt*p2.pt*p3.pt);
                std::cout << "\tR_idx = " << R_idx << std::endl;
                std::cout << "\tr_idx = " << r_idx << std::endl;
                std::cout << "\tc_idx = " << c_idx << std::endl;
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
