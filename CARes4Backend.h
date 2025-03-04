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
inline void get_CAentry_onepart(
        EEC::CAres4_entry& entry,
        const EEC::CARes4Axes& axes) noexcept {
    entry.is_symmetric = false;
    entry.is_chain = true;

    entry.fill_chain(simon::getIndex(0, axes.R),
                     0, 0, 0, axes);
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
inline void get_CAentry_twopart(
        EEC::CAres4_entry& entry22,
        EEC::CAres4_entry& entry13,
        const double dR12,
        const EEC::CARes4Axes& axes) noexcept {
    unsigned Ridx = simon::getIndex(dR12, axes.R);
    entry22.R_idx = Ridx;
    entry13.R_idx = Ridx;

    entry22.fill_symmetric(Ridx, 0, 0, 0, 0, 0, 0, axes);
    entry13.fill_chain(Ridx, 0, 0, 0, axes);
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
template <bool distances_squared>
inline void get_CAentry_threepart(
        EEC::CAres4_entry& entry1,
        EEC::CAres4_entry& entry2,
        EEC::CAres4_entry& entry3,
        const EEC::oneentry& p1,
        const EEC::oneentry& p2,
        const EEC::oneentry& p3,
        const simon::pairentry_resolved<distances_squared>& p1p2,
        const simon::pairentry_resolved<distances_squared>& p1p3,
        const simon::pairentry_resolved<distances_squared>& p2p3,
        const EEC::CARes4Axes& axes) noexcept {

    std::array<const EEC::oneentry*, 3> particles = {
        &p1, &p2, &p3
    };

    EEC::oneentry r, R;
    simon::pairentry_resolved<distances_squared> rp3;
    unsigned bkpindex;
    
    simon::CA3(particles, p1p2, p1p3, p2p3, r, R, rp3, bkpindex);

    double RL = rp3.floatDR;
    if constexpr(distances_squared){
        RL = std::sqrt(RL);
    }
    unsigned Ridx = simon::getIndex(RL, axes.R);

    double RS1 = p1p2.floatDR;
    if constexpr(distances_squared){
        RS1 = std::sqrt(RS1);
    }
    double cRr1 = simon::angle_between(p1p2, rp3);

    double RS2 = 0;
    double cRr2 = 0;

    double cr1r2 = 0;

    if (bkpindex == 0){
        entry1.fill_symmetric(
                Ridx, 
                RL, RS1, RS2, 
                cRr1, cRr2, cr1r2, 
                axes
        );
        entry2.fill_chain(Ridx, RL, RS1, cRr1, axes);
        entry3.fill_chain(Ridx, RL, RS1, cRr1, axes);
    } else if (bkpindex == 1){
        entry1.fill_chain(Ridx, RL, RS1, cRr1, axes);
        entry2.fill_symmetric(
                Ridx, 
                RL, RS1, RS2, 
                cRr1, cRr2, cr1r2, 
                axes
        );
        entry3.fill_chain(Ridx, RL, RS1, cRr1, axes);
    } else {
        entry1.fill_chain(Ridx, RL, RS1, cRr1, axes);
        entry2.fill_chain(Ridx, RL, RS1, cRr1, axes);
        entry3.fill_symmetric(
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
template <bool distances_squared>
inline void get_CAentry_fourpart(
        EEC::CAres4_entry& entry,
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

    std::array<const EEC::oneentry*, 4> particles = {
        &p1, &p2, &p3, &p4
    };
    
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

    double RL = return_pairs[0].floatDR;
    if constexpr(distances_squared){
        RL = std::sqrt(RL);
    }
    unsigned Ridx = simon::getIndex(RL, axes.R);

    if (topologyflag == simon::CHAIN){
        double RS1 = return_pairs[1].floatDR;
        if constexpr(distances_squared){
            RS1 = std::sqrt(RS1);
        }
        double c = simon::angle_between(
            return_pairs[1], return_pairs[2]
        );
        entry.fill_chain(Ridx, RL, RS1, c, axes);
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
        entry.fill_symmetric(
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

        const EEC::CAres4_entry& entry_gen,

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

                    EEC::CAres4_entry entry_reco;
                    get_CAentry_fourpart(
                        entry_reco,
                        p1, p2, p3, p4,
                        p1p2, p1p3, p1p4, 
                        p2p3, p2p4, p3p4,
                        axes_reco
                    ); 

                    transfer.fill(entry_reco, entry_gen, 
                                  twt4, wt_gen);
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

        const double E1_p2 = p1.pt*p1.pt;
        const double E1_p3 = E1_p2*p1.pt;

        /*
         * One-particle piece
         * Symmetry factor is 4!/4! = 1
         */
        double wt = E1_p3*p1.pt;

        EEC::CAres4_entry entry;
        get_CAentry_onepart(entry, axes_gen);
        result.fill(entry, wt);

        if constexpr(doUnmatched){
            if (!matched1){
                unmatched_gen->fill(entry, wt);
            }
        }

        for(unsigned i2=i1+1; i2<thisjet_gen.N; ++i2){
            const auto& p2 = thisjet_gen.singles.get(i2);
            [[maybe_unused]] bool matched2;
            if constexpr(doUnmatched){
                matched2 = matched1 && matched->at(i2);
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

            EEC::CAres4_entry entry22, entry13;
            double RL = p1p2.floatDR;
            if constexpr(JetType::pairType::distances_squared){
                RL = std::sqrt(RL);
            }
            get_CAentry_twopart(entry22, entry13, 
                                RL,
                                axes_gen);
            //we actually know that 
            //22 will always be symmetric
            //13 will always be chain
            result.fill_chain(entry13, wt_1+wt_2);
            result.fill_symmetric(entry22, wt_12);
            if constexpr(doUnmatched){
                if (!matched2){
                    unmatched_gen->fill_chain(entry13, wt_1+wt_2);
                    unmatched_gen->fill_symmetric(entry22, wt_12);
                }
            }
        
            for(unsigned i3=i2+1; i3<thisjet_gen.N; ++i3){
                const auto& p3 = thisjet_gen.singles.get(i3);
                [[maybe_unused]] bool matched3;
                if constexpr(doUnmatched){
                    matched3 = matched2 && matched->at(i3);
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
                wt_1 = E123 * p1.pt;
                wt_2 = E123 * p2.pt;
                double wt_3 = E123 * p3.pt;

                EEC::CAres4_entry entry1, entry2, entry3;
                get_CAentry_threepart(entry1, entry2, entry3,
                                      p1, p2, p3,
                                      p1p2, p1p3, p2p3,
                                      axes_gen);

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

                for(unsigned i4=i3+1; i4<thisjet_gen.N; ++i4){
                    const auto& p4 = thisjet_gen.singles.get(i4);
                    [[maybe_unused]] bool matched4;
                    if constexpr(doUnmatched){
                        matched4 = matched3 && matched->at(i4);
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
                    wt = E1234;

                    EEC::CAres4_entry entry;
                    get_CAentry_fourpart(entry,
                                        p1, p2, p3, p4,
                                        p1p2, p1p3, p1p4,
                                        p2p3, p2p4, p3p4,
                                        axes_gen);
                    result.fill(entry, wt);
                    if constexpr(doUnmatched){
                        if (!matched4){
                            unmatched_gen->fill(entry, wt);
                        }
                    }
                }
            }
        }
    }
}

#endif
