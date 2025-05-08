#ifndef SROTHMAN_EEC_CA_RES4_BACKEND_H
#define SROTHMAN_EEC_CA_RES4_BACKEND_H

#include "usings.h"

#include "SRothman/SimonTools/src/deltaR.h"
#include "SRothman/SimonTools/src/histutil.h"
#include "SRothman/SimonTools/src/util.h"
#include "SRothman/SimonTools/src/CA.h"

#include <memory>

#include <array>


/*
 * (1, 1, 1, 1) can be interpreted equally well
 * as a chain or symmetric topology
 * with R = r1 = r2 = 0
 *
 * we don't really care about the R=0 case that much,
 * so to simplify things we'll just call it a chain
 *
 * this prevents double-counting, and makes 
 * it easier to validate my implementation
 */
template <class ResultType>
inline void get_CAentry_onepart(
        EEC::CAres4_entry<typename ResultType::T>& entry,
        const EEC::CARes4Axes& axes) noexcept {
    entry.is_symmetric = false;
    entry.is_chain = true;

    if constexpr(ResultType::SHOULD_BIN){
        entry.template fill_chain<true>(simon::getIndex(0, axes.R),
                         0, 0, 0, axes);
    } else {
        entry.template fill_chain<false>(0, 0, 0, 0, axes);
    }
}

/*
 * (1, 1, 2, 2) and (1, 1, 1, 2) are logically different imo
 *
 * (1, 1, 2, 2) is definitely a valid symmetric topology
 * with r1 = r2 = 0
 * but it is NOT a valid chain topology
 * because p4 is never duplicated in the chain
 *
 * (1, 1, 1, 2), or equivalently (1, 2, 2, 2) 
 * is obviously a valid chain with r1=r2=0
 * but it is NOT a valid symmetric topology
 */
template <class ResultType, bool distances_squared>
inline void get_CAentry_twopart(
        EEC::CAres4_entry<typename ResultType::T>& entry22,
        EEC::CAres4_entry<typename ResultType::T>& entry13,
        const double dR12,
        const EEC::CARes4Axes& axes) noexcept {
    typename ResultType::T Ridx;
    if constexpr (ResultType::SHOULD_BIN){
        if constexpr(distances_squared){
            Ridx = simon::getIndex(std::sqrt(dR12), axes.R);
        } else {
            Ridx = simon::getIndex(dR12, axes.R);
        }
    } else {
        if constexpr(distances_squared){
            Ridx = std::sqrt(dR12);
        } else {
            Ridx = dR12;
        }
    }
    entry22.R_idx = Ridx;
    entry13.R_idx = Ridx;

    entry22.template fill_symmetric<ResultType::SHOULD_BIN>(Ridx, 0, 0, 0, 0, 0, 0, axes);
    entry13.template fill_chain<ResultType::SHOULD_BIN>(Ridx, 0, 0, 0, axes);
}

/*
 * For the three-particle part it 
 * actually matters which particle is p3
 *
 * if p3 is duplicated, then it is 
 * a symmetric topology with r2 = 0
 *
 * if p1 or p2 are duplicated then it is 
 * a chain with r2 = 0
 */
template <class ResultType, bool distances_squared>
inline void get_CAentry_threepart(
        EEC::CAres4_entry<typename ResultType::T>& entry1,
        EEC::CAres4_entry<typename ResultType::T>& entry2,
        EEC::CAres4_entry<typename ResultType::T>& entry3,
        const EEC::oneentry& p1,
        const EEC::oneentry& p2,
        const EEC::oneentry& p3,
        const simon::pairentry_resolved<distances_squared>& p1p2,
        const simon::pairentry_resolved<distances_squared>& p1p3,
        const simon::pairentry_resolved<distances_squared>& p2p3,
        const EEC::CARes4Axes& axes) noexcept {

#ifdef CHECK_BY_HAND
    printf("running CAentry_threepart\n");
#endif
    std::array<const EEC::oneentry*, 3> particles = {{
        &p1, &p2, &p3
    }};
#ifdef CHECK_BY_HAND
    printf("p1 [pt, eta, phi] = [%g, %g, %g]\n", p1.pt, p1.eta, p1.phi);
    printf("p2 [pt, eta, phi] = [%g, %g, %g]\n", p2.pt, p2.eta, p2.phi);
    printf("p3 [pt, eta, phi] = [%g, %g, %g]\n", p3.pt, p3.eta, p3.phi);
#endif

    EEC::oneentry r, R;
    simon::pairentry_resolved<distances_squared> rp3;
    unsigned bkpindex;
    
    simon::CA3(particles, 
               p1p2, p1p3, p2p3, 
               r, R, rp3, 
               bkpindex);

#ifdef CHECK_BY_HAND
    printf("r [pt, eta, phi] = [%g, %g, %g]\n", r.pt, r.eta, r.phi);
    printf("R [pt, eta, phi] = [%g, %g, %g]\n", R.pt, R.eta, R.phi);
#endif

    double RL = rp3.floatDR;
    if constexpr(distances_squared){
        RL = std::sqrt(RL);
    }
    typename ResultType::T Ridx;
    if constexpr (ResultType::SHOULD_BIN){
        Ridx = simon::getIndex(RL, axes.R);
    } else {
        Ridx = RL;
    }

    double RS1; 
    double cRr1; 
    if (bkpindex == 2){
        RS1 = p1p2.floatDR;
        cRr1 = simon::angle_between(p1p2, rp3);
    } else if (bkpindex == 1){
        RS1 = p1p3.floatDR;
        cRr1 = simon::angle_between(p1p3, rp3);
    } else {
        RS1 = p2p3.floatDR;
        cRr1 = simon::angle_between(p2p3, rp3);
    }
    if constexpr(distances_squared){
        RS1 = std::sqrt(RS1);
    }

    double RS2 = 0;
    double cRr2 = 0;

    double cr1r2 = 0;

#ifdef CHECK_BY_HAND
    printf("RL: %g\n", RL);
    printf("RS1: %g\n", RS1);
    printf("cRr1: %g\n", cRr1);
#endif

    if (bkpindex == 0){
        entry1.template fill_symmetric<ResultType::SHOULD_BIN>(
                Ridx, 
                RL, RS1, RS2, 
                cRr1, cRr2, cr1r2, 
                axes
        );
        entry2.template fill_chain<ResultType::SHOULD_BIN>(Ridx, RL, RS1, cRr1, axes);
        entry3.template fill_chain<ResultType::SHOULD_BIN>(Ridx, RL, RS1, cRr1, axes);
    } else if (bkpindex == 1){
        entry1.template fill_chain<ResultType::SHOULD_BIN>(Ridx, RL, RS1, cRr1, axes);
        entry2.template fill_symmetric<ResultType::SHOULD_BIN>(
                Ridx, 
                RL, RS1, RS2, 
                cRr1, cRr2, cr1r2, 
                axes
        );
        entry3.template fill_chain<ResultType::SHOULD_BIN>(Ridx, RL, RS1, cRr1, axes);
    } else {
        entry1.template fill_chain<ResultType::SHOULD_BIN>(Ridx, RL, RS1, cRr1, axes);
        entry2.template fill_chain<ResultType::SHOULD_BIN>(Ridx, RL, RS1, cRr1, axes);
        entry3.template fill_symmetric<ResultType::SHOULD_BIN>(
                Ridx, 
                RL, RS1, RS2, 
                cRr1, cRr2, cr1r2, 
                axes
        );
    }
}

/*
 * The four-particle part is unique, making my life easier :)
 */
template <class ResultType, bool distances_squared>
inline void get_CAentry_fourpart(
        EEC::CAres4_entry<typename ResultType::T>& entry,
        const EEC::oneentry& p1,
        const EEC::oneentry& p2,
        const EEC::oneentry& p3,
        const EEC::oneentry& p4,
        const simon::pairentry_resolved<distances_squared>& p1p2,
        const simon::pairentry_resolved<distances_squared>& p1p3,
        const simon::pairentry_resolved<distances_squared>& p1p4,
        const simon::pairentry_resolved<distances_squared>& p2p3,
        const simon::pairentry_resolved<distances_squared>& p2p4,
        const simon::pairentry_resolved<distances_squared>& p3p4,
        const EEC::CARes4Axes& axes) noexcept {

    std::array<const EEC::oneentry*, 4> particles = {{
        &p1, &p2, &p3, &p4
    }};
#ifdef CHECK_BY_HAND 
    printf("running CAentry_fourpart\n");
    printf("p1 [pt, eta, phi] = [%g, %g, %g]\n", p1.pt, p1.eta, p1.phi);
    printf("p2 [pt, eta, phi] = [%g, %g, %g]\n", p2.pt, p2.eta, p2.phi);
    printf("p3 [pt, eta, phi] = [%g, %g, %g]\n", p3.pt, p3.eta, p3.phi);
    printf("p4 [pt, eta, phi] = [%g, %g, %g]\n", p4.pt, p4.eta, p4.phi);
#endif
    EEC::oneentry r1, r2, R;
    std::vector<simon::pairentry_resolved<distances_squared>> return_pairs;
    int topologyflag;
    std::array<unsigned, 4> return_indices; 
    std::vector<const simon::pairentry_resolved<distances_squared> *> return_orig_pairs;

    simon::CA4(
            particles, 
            p1p2, p1p3, p1p4,
            p2p3, p2p4, p3p4, 
            r1, r2, R, 
            return_pairs, 
            topologyflag, 
            return_indices,
            return_orig_pairs
    );
#ifdef CHECK_BY_HAND
    printf("r1 [pt, eta, phi] = [%g, %g, %g]\n", r1.pt, r1.eta, r1.phi);
    printf("r2 [pt, eta, phi] = [%g, %g, %g]\n", r2.pt, r2.eta, r2.phi);
    printf("R [pt, eta, phi] = [%g, %g, %g]\n", R.pt, R.eta, R.phi);
#endif

    double RL = return_pairs[0].floatDR;
    if constexpr(distances_squared){
        RL = std::sqrt(RL);
    }
    typename ResultType::T Ridx;
    if constexpr(ResultType::SHOULD_BIN){
        Ridx = simon::getIndex(RL, axes.R);
    } else {
        Ridx = RL;
    }

#ifdef CHECK_BY_HAND
    printf("RL: %g\n", RL);
#endif

    if (RL==0){
#ifdef CHECK_BY_HAND
        printf("RL is zero - short-circuiting\n");
#endif
        entry.template fill_chain<ResultType::SHOULD_BIN>(Ridx, RL, 0, 0, axes);
        return;
    }

    if (topologyflag == simon::CHAIN){
        double RS1 = return_pairs[1].floatDR;
        if constexpr(distances_squared){
            RS1 = std::sqrt(RS1);
        }
        double c = simon::angle_between(
            return_pairs[0], return_pairs[1]
        );
        entry.template fill_chain<ResultType::SHOULD_BIN>(Ridx, RL, RS1, c, axes);
#ifdef CHECK_BY_HAND
        printf("chain topology\n");
        printf("RS1: %g\n", RS1);
        printf("c: %g\n", c);
#endif
    } else {
        double RS2 = return_orig_pairs[0]->floatDR;
        double RS1 = return_orig_pairs[1]->floatDR;
        if constexpr(distances_squared){
            RS1 = std::sqrt(RS1);
            RS2 = std::sqrt(RS2);
        }
        double cRr1 = simon::angle_between(
            *return_orig_pairs[1], return_pairs[0]
        );
        double cRr2 = simon::angle_between(
            *return_orig_pairs[0], return_pairs[0]
        );
        double cr1r2 = simon::angle_between(
            *return_orig_pairs[0], *return_orig_pairs[1]
        );
#ifdef CHECK_BY_HAND
        printf("symmetric topology\n");
        printf("RS1: %g\n", RS1);
        printf("RS2: %g\n", RS2);
        printf("cRr1: %g\n", cRr1);
        printf("cRr2: %g\n", cRr2);
        printf("cr1r2: %g\n", cr1r2);
#endif
        entry.template fill_symmetric<ResultType::SHOULD_BIN>(
            Ridx, 
            RL, RS1, RS2, 
            cRr1, cRr2, cr1r2, 
            axes
        );
    }
}

template <class TransferResultType, class ResultType, class JetType>
inline void CAres4_transferloop(
        TransferResultType& transfer,

        const EEC::CAres4_entry<typename ResultType::T>& entry_gen,

        const JetType& thisjet_reco,

        const EEC::neighborhood& n1,
        const EEC::neighborhood& n2,
        const EEC::neighborhood& n3,
        const EEC::neighborhood& n4,

        const EEC::CARes4Axes& axes_reco,

        double wt_gen) noexcept {

    for(const EEC::neighbor& j1:  n1){
        const double twt1 = wt_gen * j1.wt;
        const auto& p1 = thisjet_reco.singles.get(j1.idx);

        for(const EEC::neighbor& j2: n2){
            const double twt2 = twt1 * j2.wt;
            const auto& p2 = thisjet_reco.singles.get(j2.idx);

            const auto& p1p2 = thisjet_reco.pairs.get(j1.idx, j2.idx);

            for(const EEC::neighbor& j3: n3){
                const double twt3 = twt2 * j3.wt;
                const auto& p3 = thisjet_reco.singles.get(j3.idx);
    
                const auto& p1p3 = thisjet_reco.pairs.get(j1.idx, j3.idx);
                const auto& p2p3 = thisjet_reco.pairs.get(j2.idx, j3.idx);

                for(const EEC::neighbor& j4: n4){
                    const double twt4 = twt3 * j4.wt;
                    const auto& p4 = thisjet_reco.singles.get(j4.idx);

                    const auto& p1p4 = thisjet_reco.pairs.get(j1.idx, j4.idx);
                    const auto& p2p4 = thisjet_reco.pairs.get(j2.idx, j4.idx);
                    const auto& p3p4 = thisjet_reco.pairs.get(j3.idx, j4.idx);

                    EEC::CAres4_entry<typename ResultType::T> entry_reco;
                    get_CAentry_fourpart<ResultType>(
                        entry_reco,
                        p1, p2, p3, p4,
                        p1p2, p1p3, p1p4, 
                        p2p3, p2p4, p3p4,
                        axes_reco
                    ); 

                    transfer.fill(entry_reco, entry_gen, 
                                  twt4, wt_gen);
#ifdef CHECK_BY_HAND
                    printf("transferring %g -> %g\n",
                            wt_gen, twt4);
                    if (entry_gen.is_chain && entry_reco.is_chain){
                        printf("\tchain -> chain\n");
                    } else if (entry_gen.is_chain){
                        printf("\tchain -> symmetric\n");
                    } else if (entry_reco.is_chain){
                        printf("\tsymmetric -> chain\n");
                    } else {
                        printf("\tsymmetric -> symmetric\n");
                    }
#endif
                }
            }
        }
    }
}

template <class ResultType, class JetType, bool doUnmatched, class TransferResultType, bool doTransfer>
inline void CAres4_mainloop(
        ResultType& result,
        [[maybe_unused]] ResultType* unmatched_gen,
        [[maybe_unused]] TransferResultType* transfer,
        [[maybe_unused]] const JetType * const thisjet_reco,
        const JetType& thisjet_gen,
        [[maybe_unused]] const std::vector<bool> * const matched,
        [[maybe_unused]] const EEC::Adjacency* const adj,
        [[maybe_unused]] const EEC::CARes4Axes* const axes_reco,
        const EEC::CARes4Axes& axes_gen) noexcept {

    for(unsigned i1=0; i1<thisjet_gen.N; ++i1){
        const auto& p1 = thisjet_gen.singles.get(i1);
        [[maybe_unused]] bool matched1;
        if constexpr(doUnmatched){
            matched1 = matched->at(i1);
        }
        [[maybe_unused]] EEC::neighborhood const * n1 = nullptr;
        if constexpr(doTransfer){
            n1 = &adj->get_neighborhood(i1);
        }

        const double E1_p2 = p1.pt*p1.pt;
        const double E1_p3 = E1_p2*p1.pt;

        /*
         * One-particle piece
         * Symmetry factor is 4!/4! = 1
         */
        double wt = E1_p3*p1.pt;

        EEC::CAres4_entry<typename ResultType::T> entry;
        get_CAentry_onepart<ResultType>(entry, axes_gen);
#ifdef CHECK_BY_HAND
        printf("1-particle part\n");
        printf("(%g %g)\n", p1.eta, p1.phi);
        printf("E1: %g\n", p1.pt);
        //printf("\tR: %u\n", entry.R_idx);
        //printf("\tis_symmetric: %u\n", entry.is_symmetric);
        //printf("\t\tr_wrtr: %u\n", entry.symmetric_wrtr_r_idx);
        //printf("\t\tc_wrtr: %u\n", entry.symmetric_wrtr_c_idx);
        //printf("\t\tr1_wrtR: %u\n", entry.symmetric_wrtR1_r_idx);
        //printf("\t\tc1_wrtR: %u\n", entry.symmetric_wrtR1_c_idx);
        //printf("\t\tr2_wrtR: %u\n", entry.symmetric_wrtR2_r_idx);
        //printf("\t\tc2_wrtR: %u\n", entry.symmetric_wrtR2_c_idx);
        //printf("\tis_chain: %u\n", entry.is_chain);
        //printf("\t\tr_chain: %u\n", entry.chain_r_idx);
        //printf("\t\tc_chain: %u\n", entry.chain_c_idx);
        printf("\twt: %f\n", wt);
#endif
        result.fill(entry, wt);

        if constexpr(doUnmatched){
            if (!matched1){
                unmatched_gen->fill(entry, wt);
            }
        }

        if constexpr(doTransfer){
            CAres4_transferloop<TransferResultType, ResultType, JetType>(
                *transfer,
                entry,
                *thisjet_reco,
                *n1,
                *n1,
                *n1,
                *n1,
                *axes_reco,
                wt
            );
        }

        for(unsigned i2=i1+1; i2<thisjet_gen.N; ++i2){
            const auto& p2 = thisjet_gen.singles.get(i2);
            [[maybe_unused]] bool matched2;
            if constexpr(doUnmatched){
                matched2 = matched1 && matched->at(i2);
            }
            [[maybe_unused]] EEC::neighborhood const * n2 = nullptr;
            if constexpr(doTransfer){
                n2 = &adj->get_neighborhood(i2);
            }

            const auto& p1p2 = thisjet_gen.pairs.get(i1, i2);

            const double E12 = p1.pt*p2.pt;

            /*
             * Two-particle piece
             * Has three different components
             * with two different symmetry factors
             *
             * (1,1,1,2) and (1, 2, 2, 2) have 
             * symmetry factor 4!/3! = 4
             *
             * (1,1,2,2) has symmetry factor 4!/(2!x2!) = 6
             */

            double wt_1 = 4*E12 * E1_p2;
            double wt_2 = 4*E12 * p2.pt*p2.pt;
            
            double wt_12 = 6*E12 * E12;

            EEC::CAres4_entry<typename ResultType::T> entry22, entry13;
            double RL = p1p2.floatDR;
            get_CAentry_twopart<ResultType, JetType::pairType::distances_squared>(
                    entry22, entry13, 
                    RL, axes_gen
            );
            //we actually know that 
            //22 will always be symmetric
            //13 will always be chain
            result.template fill_chain<typename ResultType::T>(entry13, wt_1+wt_2);
            result.template fill_symmetric<typename ResultType::T>(entry22, wt_12);
#ifdef CHECK_BY_HAND
            printf("2-particle part\n");
            printf("(%g %g), (%g %g)\n", p1.eta, p1.phi, p2.eta, p2.phi);
            printf("E1: %g\n", p1.pt);
            printf("E2: %g\n", p2.pt);
            printf("\t1112 + 12222:\n");
            //printf("\t\tR: %u\n", entry13.R_idx);
            //printf("\t\tis_symmetric: %u\n", entry13.is_symmetric);
            //printf("\t\t\tr_wrtr: %u\n", entry13.symmetric_wrtr_r_idx);
            //printf("\t\t\tc_wrtr: %u\n", entry13.symmetric_wrtr_c_idx);
            //printf("\t\t\tr1_wrtR: %u\n", entry13.symmetric_wrtR1_r_idx);
            //printf("\t\t\tc1_wrtR: %u\n", entry13.symmetric_wrtR1_c_idx);
            //printf("\t\t\tr2_wrtR: %u\n", entry13.symmetric_wrtR2_r_idx);
            //printf("\t\t\tc2_wrtR: %u\n", entry13.symmetric_wrtR2_c_idx);
            //printf("\t\tis_chain: %u\n", entry13.is_chain);
            //printf("\t\t\tr_chain: %u\n", entry13.chain_r_idx);
            //printf("\t\t\tc_chain: %u\n", entry13.chain_c_idx);
            printf("\t\twt: %f\n", wt_1+wt_2);
            printf("\t1122:\n");
            //printf("\t\tR: %u\n", entry22.R_idx);
            //printf("\t\tis_symmetric: %u\n", entry22.is_symmetric);
            //printf("\t\t\tr_wrtr: %u\n", entry22.symmetric_wrtr_r_idx);
            //printf("\t\t\tc_wrtr: %u\n", entry22.symmetric_wrtr_c_idx);
            //printf("\t\t\tr1_wrtR: %u\n", entry22.symmetric_wrtR1_r_idx);
            //printf("\t\t\tc1_wrtR: %u\n", entry22.symmetric_wrtR1_c_idx);
            //printf("\t\t\tr2_wrtR: %u\n", entry22.symmetric_wrtR2_r_idx);
            //printf("\t\t\tc2_wrtR: %u\n", entry22.symmetric_wrtR2_c_idx);
            //printf("\t\tis_chain: %u\n", entry22.is_chain);
            //printf("\t\t\tr_chain: %u\n", entry22.chain_r_idx);
            //printf("\t\t\tc_chain: %u\n", entry22.chain_c_idx);
            printf("\t\twt: %f\n", wt_12);
#endif
            if constexpr(doUnmatched){
                if (!matched2){
                    unmatched_gen->template fill_chain<typename ResultType::T>(entry13, wt_1+wt_2);
                    unmatched_gen->template fill_symmetric<typename ResultType::T>(entry22, wt_12);
                }
            }

            if constexpr(doTransfer){
                CAres4_transferloop<TransferResultType, ResultType, JetType>(
                    *transfer,
                    entry13,
                    *thisjet_reco,
                    *n1,
                    *n1,
                    *n1,
                    *n2,
                    *axes_reco,
                    wt_1
                );
                CAres4_transferloop<TransferResultType, ResultType, JetType>(
                    *transfer,
                    entry13,
                    *thisjet_reco,
                    *n1,
                    *n2,
                    *n2,
                    *n2,
                    *axes_reco,
                    wt_2
                );
                CAres4_transferloop<TransferResultType, ResultType, JetType>(
                    *transfer,
                    entry22,
                    *thisjet_reco,
                    *n1,
                    *n1,
                    *n2,
                    *n2,
                    *axes_reco,
                    wt_12
                );
            }
        
            for(unsigned i3=i2+1; i3<thisjet_gen.N; ++i3){
                const auto& p3 = thisjet_gen.singles.get(i3);
                [[maybe_unused]] bool matched3;
                if constexpr(doUnmatched){
                    matched3 = matched2 && matched->at(i3);
                }
                [[maybe_unused]] EEC::neighborhood const * n3 = nullptr;
                if constexpr(doTransfer){
                    n3 = &adj->get_neighborhood(i3);
                }

                const auto& p1p3 = thisjet_gen.pairs.get(i1, i3);
                const auto& p2p3 = thisjet_gen.pairs.get(i2, i3);

                const double E123 = E12*p3.pt;

                /*
                 * Three-particle piece
                 *
                 * This has three different components
                 * all with the same symmetry factor
                 * of 4! / 2! = 12
                 */
                wt_1 = 12 * E123 * p1.pt;
                wt_2 = 12 * E123 * p2.pt;
                double wt_3 = 12 * E123 * p3.pt;

                EEC::CAres4_entry<typename ResultType::T> entry1, entry2, entry3;
                get_CAentry_threepart<ResultType, JetType::pairType::distances_squared>(
                        entry1, entry2, entry3,
                        p1, p2, p3,
                        p1p2, p1p3, p2p3,
                        axes_gen);

#ifdef CHECK_BY_HAND
                printf("3-particle part\n");
                printf("(%g %g), (%g %g), (%g %g)\n", p1.eta, p1.phi, p2.eta, p2.phi, p3.eta, p3.phi);
                printf("E1: %g\n", p1.pt);
                printf("E2: %g\n", p2.pt);
                printf("E3: %g\n", p3.pt);
                printf("\t1123:\n");
                //printf("\t\tR: %u\n", entry1.R_idx);
                //printf("\t\tis_symmetric: %u\n", entry1.is_symmetric);
                //printf("\t\t\tr_wrtr: %u\n", entry1.symmetric_wrtr_r_idx);
                //printf("\t\t\tc_wrtr: %u\n", entry1.symmetric_wrtr_c_idx);
                //printf("\t\t\tr1_wrtR: %u\n", entry1.symmetric_wrtR1_r_idx);
                //printf("\t\t\tc1_wrtR: %u\n", entry1.symmetric_wrtR1_c_idx);
                //printf("\t\t\tr2_wrtR: %u\n", entry1.symmetric_wrtR2_r_idx);
                //printf("\t\t\tc2_wrtR: %u\n", entry1.symmetric_wrtR2_c_idx);
                //printf("\t\tis_chain: %u\n", entry1.is_chain);
                //printf("\t\t\tr_chain: %u\n", entry1.chain_r_idx);
                //printf("\t\t\tc_chain: %u\n", entry1.chain_c_idx);
                printf("\t\twt: %f\n", wt_1);
                printf("\t1223:\n");
                //printf("\t\tR: %u\n", entry2.R_idx);
                //printf("\t\tis_symmetric: %u\n", entry2.is_symmetric);
                //printf("\t\t\tr_wrtr: %u\n", entry2.symmetric_wrtr_r_idx);
                //printf("\t\t\tc_wrtr: %u\n", entry2.symmetric_wrtr_c_idx);
                //printf("\t\t\tr1_wrtR: %u\n", entry2.symmetric_wrtR1_r_idx);
                //printf("\t\t\tc1_wrtR: %u\n", entry2.symmetric_wrtR1_c_idx);
                //printf("\t\t\tr2_wrtR: %u\n", entry2.symmetric_wrtR2_r_idx);
                //printf("\t\t\tc2_wrtR: %u\n", entry2.symmetric_wrtR2_c_idx);
                //printf("\t\tis_chain: %u\n", entry2.is_chain);
                //printf("\t\t\tr_chain: %u\n", entry2.chain_r_idx);
                //printf("\t\t\tc_chain: %u\n", entry2.chain_c_idx);
                printf("\t\twt: %f\n", wt_2);
                printf("\t1233:\n");
                //printf("\t\tR: %u\n", entry3.R_idx);
                //printf("\t\tis_symmetric: %u\n", entry3.is_symmetric);
                //printf("\t\t\tr_wrtr: %u\n", entry3.symmetric_wrtr_r_idx);
                //printf("\t\t\tc_wrtr: %u\n", entry3.symmetric_wrtr_c_idx);
                //printf("\t\t\tr1_wrtR: %u\n", entry3.symmetric_wrtR1_r_idx);
                //printf("\t\t\tc1_wrtR: %u\n", entry3.symmetric_wrtR1_c_idx);
                //printf("\t\t\tr2_wrtR: %u\n", entry3.symmetric_wrtR2_r_idx);
                //printf("\t\t\tc2_wrtR: %u\n", entry3.symmetric_wrtR2_c_idx);
                //printf("\t\tis_chain: %u\n", entry3.is_chain);
                //printf("\t\t\tr_chain: %u\n", entry3.chain_r_idx);
                //printf("\t\t\tc_chain: %u\n", entry3.chain_c_idx);
                printf("\t\twt: %f\n", wt_3);
#endif
                result.fill(entry1, wt_1);
                result.fill(entry2, wt_2);
                result.fill(entry3, wt_3);
                if constexpr(doUnmatched){
                    if (!matched3){
                        unmatched_gen->fill(entry1, wt_1);
                        unmatched_gen->fill(entry2, wt_2);
                        unmatched_gen->fill(entry3, wt_3);
                    }
                }

                if constexpr(doTransfer){
                    CAres4_transferloop<TransferResultType, ResultType, JetType>(
                        *transfer,
                        entry1,
                        *thisjet_reco,
                        *n1,
                        *n1,
                        *n2,
                        *n3,
                        *axes_reco,
                        wt_1
                    );
                    CAres4_transferloop<TransferResultType, ResultType, JetType>(
                        *transfer,
                        entry2,
                        *thisjet_reco,
                        *n1,
                        *n2,
                        *n2,
                        *n3,
                        *axes_reco,
                        wt_2
                    );  
                    CAres4_transferloop<TransferResultType, ResultType, JetType>(
                        *transfer,
                        entry3,
                        *thisjet_reco,
                        *n1,
                        *n2,
                        *n3,
                        *n3,
                        *axes_reco,
                        wt_3
                    );
                }

                for(unsigned i4=i3+1; i4<thisjet_gen.N; ++i4){
                    const auto& p4 = thisjet_gen.singles.get(i4);
                    [[maybe_unused]] bool matched4;
                    if constexpr(doUnmatched){
                        matched4 = matched3 && matched->at(i4);
                    }
                    [[maybe_unused]] EEC::neighborhood const * n4 = nullptr;
                    if constexpr(doTransfer){
                        n4 = &adj->get_neighborhood(i4);
                    }

                    const auto& p1p4 = thisjet_gen.pairs.get(i1, i4);
                    const auto& p2p4 = thisjet_gen.pairs.get(i2, i4);
                    const auto& p3p4 = thisjet_gen.pairs.get(i3, i4);

                    const double E1234 = E123*p4.pt;

                    /*
                     * Four-particle piece
                     *
                     * This has one component 
                     * with symmetry factor
                     * 4! = 24
                     */
                    wt = 24 * E1234;

                    EEC::CAres4_entry<typename ResultType::T> entry;
                    get_CAentry_fourpart<ResultType, JetType::pairType::distances_squared>(
                                    entry,
                                    p1, p2, p3, p4,
                                    p1p2, p1p3, p1p4,
                                    p2p3, p2p4, p3p4,
                                    axes_gen);
#ifdef CHECK_BY_HAND
                    printf("4-particle part\n");
                    printf("(%g %g), (%g %g), (%g %g), (%g %g)\n", p1.eta, p1.phi, p2.eta, p2.phi, p3.eta, p3.phi, p4.eta, p4.phi);
                    printf("E1: %g\n", p1.pt);
                    printf("E2: %g\n", p2.pt);
                    printf("E3: %g\n", p3.pt);
                    printf("E4: %g\n", p4.pt);
                    //printf("\tR: %u\n", entry.R_idx);
                    //printf("\tis_symmetric: %u\n", entry.is_symmetric);
                    //printf("\t\tr_wrtr: %u\n", entry.symmetric_wrtr_r_idx);
                    //printf("\t\tc_wrtr: %u\n", entry.symmetric_wrtr_c_idx);
                    //printf("\t\tr1_wrtR: %u\n", entry.symmetric_wrtR1_r_idx);
                    //printf("\t\tc1_wrtR: %u\n", entry.symmetric_wrtR1_c_idx);
                    //printf("\t\tr2_wrtR: %u\n", entry.symmetric_wrtR2_r_idx);
                    //printf("\t\tc2_wrtR: %u\n", entry.symmetric_wrtR2_c_idx);
                    //printf("\tis_chain: %u\n", entry.is_chain);
                    //printf("\t\tr_chain: %u\n", entry.chain_r_idx);
                    //printf("\t\tc_chain: %u\n", entry.chain_c_idx);
                    printf("\twt: %f\n", wt);
#endif
                    result.fill(entry, wt);
                    if constexpr(doUnmatched){
                        if (!matched4){
                            unmatched_gen->fill(entry, wt);
                        }
                    }
                    if constexpr(doTransfer){
                        CAres4_transferloop<TransferResultType, ResultType, JetType>(
                            *transfer,
                            entry,
                            *thisjet_reco,
                            *n1,
                            *n2,
                            *n3,
                            *n4,
                            *axes_reco,
                            wt
                        );
                    }
                }
            }
        }
    }
}

#endif
