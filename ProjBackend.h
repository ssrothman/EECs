#ifndef SROTHMN_EEC_PROJ_BACKEND_H
#define SROTHMN_EEC_PROJ_BACKEND_H

#include "Adjacency.h"

#include "SRothman/SimonTools/src/histutil.h"
#include "SRothman/SimonTools/src/util.h"

#include <memory>

constexpr unsigned NWT(unsigned nneighbor){
    if (nneighbor == 1) return 5;
    else return 7 - nneighbor;
}

template <unsigned N, class JetType, class TransferResultType>
inline void proj_transferloop(
        TransferResultType& transfer,

        const std::shared_ptr<const JetType> thisjet_reco,

        const std::array<const EEC::neighborhood *, N>& ns,

        const double dR_gen,
        const std::array<double, NWT(N)>& wt_gen) noexcept {

    for (unsigned i=0; i<N; ++i){
        if (ns[i]->empty()){
            return;
        } else if (ns[i]->size() != 1){
            printf("ERROR: only 1-to-1 matching is supported! Results will be incomplete\n");
        }
    }

    typename TransferResultType::T dR = 0;
    double basewt = ns[0]->at(0).wt;
    double wtsum = basewt;
    if constexpr(N > 1){
        for (unsigned i=1; i<N; ++i){
            basewt *= ns[i]->at(0).wt;
            wtsum += ns[i]->at(0).wt;
            for (unsigned j=0; j<i; ++j){
                dR = std::max(dR, thisjet_reco->pairs.getDR(
                            ns[j]->at(0).idx, ns[i]->at(0).idx
                ));
            }
        }
    }

    printf("Transfer %u\n", N);

    double partial = basewt;
    if constexpr (N == 1){
        partial *= wtsum;
    }
    printf("\tbasewt = %g\n", basewt);
    printf("\twtsum = %g\n", wtsum);

    //transfer.fill(2, dR, dR_gen, wt_gen[0], wt_gen[0]);

    transfer.fill(
        7-NWT(N),
        dR, dR_gen,
        partial*wt_gen[0], wt_gen[0]
    );
    printf("\tTransfer fill %u:\n", 7-NWT(N));
    std::cout << "\t\t" << dR_gen << " --> "<< dR << std::endl;
    printf("\t\t%g --> %g\n", wt_gen[0], partial*wt_gen[0]);

    for (unsigned order=7-NWT(N)+1; order<7; ++order){
        unsigned i = order - (7-NWT(N));
        partial *= wtsum;

        printf("\tTransfer fill %u: (i = %u)\n", order, i);
        printf("\t\t%g --> %g\n", wt_gen[i], partial*wt_gen[i]);

        transfer.fill(
            order,
            dR, dR_gen,
            partial*wt_gen[i], wt_gen[i]
        );
    }
}

template <class ResultType, class JetType, bool doUnmatched, class TransferResultType, bool doTransfer>
inline void proj_mainloop(
        ResultType& result,

        [[maybe_unused]] ResultType* unmatched_gen,

        [[maybe_unused]] TransferResultType* transfer,

        [[maybe_unused]] const std::shared_ptr<const JetType> thisjet_reco,
        const std::shared_ptr<const JetType> thisjet_gen,
        [[maybe_unused]] const std::vector<bool>* const matched,

        [[maybe_unused]] const std::shared_ptr<const EEC::Adjacency> adj){

    for (unsigned i1=0; i1<thisjet_gen->N; ++i1){
        const double& E1 = thisjet_gen->singles.getE(i1);
        [[maybe_unused]] bool matched1;
        if constexpr(doUnmatched){
            matched1 = matched->at(i1);
        }
        [[maybe_unused]] EEC::neighborhood const * n1 = nullptr;
        if constexpr(doTransfer){
            n1 = &(adj->get_neighborhood(i1));

        }

        /*
         * The one-particle part
         * dR = 0
         * symmetry factor = 1 for all orders
         */

        double wt2 = E1*E1;
        double wt3 = E1*wt2;
        double wt4 = E1*wt3;
        double wt5 = E1*wt4;
        double wt6 = E1*wt5;

        result.template fill<2>(0, wt2);
        result.template fill<3>(0, wt3);
        result.template fill<4>(0, wt4);
        result.template fill<5>(0, wt5);
        result.template fill<6>(0, wt6);

        if constexpr (doUnmatched){
            if (!matched1){
                unmatched_gen->template fill<2>(0, wt2);
                unmatched_gen->template fill<3>(0, wt3);
                unmatched_gen->template fill<4>(0, wt4);
                unmatched_gen->template fill<5>(0, wt5);
                unmatched_gen->template fill<6>(0, wt6);
            }
        }

        if constexpr(doTransfer){
            proj_transferloop<1>(
                *transfer,
                thisjet_reco,
                {n1},
                0,
                {wt2, wt3, wt4, wt5, wt6}
            );
        }

        for(unsigned i2=i1+1; i2<thisjet_gen->N; ++i2){
            const double& E2 = thisjet_gen->singles.getE(i2);
            [[maybe_unused]] bool matched2;
            if constexpr(doUnmatched){
                matched2 = matched1 && matched->at(i2);
            }
            [[maybe_unused]] EEC::neighborhood const * n2 = nullptr;
            if constexpr(doTransfer){
                n2 = &(adj->get_neighborhood(i2));
            }

            const typename ResultType::T& dR12 = thisjet_gen->pairs.getDR(i1, i2);

            const double E12 = E1*E2;

            /*
             * two-particle part
             */
            wt2 = 2*E12;
            wt3 = 3*E12*(E1+E2);
            wt4 = 2*E12*(2*(E1*E1+E2*E2) + 3*E12);
            wt5 = 5*E12*(E1+E2)*(E1*E1 + E1*E2 + E2*E2);
            wt6 = E12*(6*(E1*E1*E1*E1+E2*E2*E2*E2) + 15*(E1*E1*E1*E2+E1*E2*E2*E2) + 20*E1*E1*E2*E2);

            result.template fill<2>(dR12, wt2);
            result.template fill<3>(dR12, wt3);
            result.template fill<4>(dR12, wt4);
            result.template fill<5>(dR12, wt5);
            result.template fill<6>(dR12, wt6);

            if constexpr (doUnmatched){
                if (!matched2){
                    unmatched_gen->template fill<2>(dR12, wt2);
                    unmatched_gen->template fill<3>(dR12, wt3);
                    unmatched_gen->template fill<4>(dR12, wt4);
                    unmatched_gen->template fill<5>(dR12, wt5);
                    unmatched_gen->template fill<6>(dR12, wt6);
                }
            }

            if constexpr(doTransfer){
                proj_transferloop<2>(
                    *transfer,
                    thisjet_reco,
                    {n1, n2},
                    dR12,
                    {wt2, wt3, wt4, wt5, wt6}
                );
            }

            for(unsigned i3=i2+1; i3<thisjet_gen->N; ++i3){
                const double& E3 = thisjet_gen->singles.getE(i3);
                [[maybe_unused]] bool matched3;
                if constexpr(doUnmatched){
                    matched3 = matched2 && matched->at(i3);
                }
                [[maybe_unused]] EEC::neighborhood const * n3 = nullptr;
                if constexpr(doTransfer){
                    n3 = &(adj->get_neighborhood(i3));
                }

                const typename ResultType::T& dR13 = thisjet_gen->pairs.getDR(i1, i3);
                const typename ResultType::T& dR23 = thisjet_gen->pairs.getDR(i2, i3);

                const double E123 = E12*E3;
                const typename ResultType::T& dR123 = std::max({
                    dR12,
                    dR13,
                    dR23
                });

                /*
                 * three-particle part
                 */
                wt3 = 6*E123;
                wt4 = 12*E123*(E1+E2+E3);
                wt5 = 10*E123*(2*(E1*E1+E2*E2+E3*E3) + 3*(E1*E2+E1*E3+E2*E3));
                wt6 = 30*E123*(E1+E2+E3)*(E1*E1+E2*E2+E3*E3+E1*E2+E1*E3+E2*E3);

                result.template fill<3>(dR123, wt3);
                result.template fill<4>(dR123, wt4);
                result.template fill<5>(dR123, wt5);
                result.template fill<6>(dR123, wt6);

                if constexpr (doUnmatched){
                    if (!matched3){
                        unmatched_gen->template fill<3>(dR123, wt3);
                        unmatched_gen->template fill<4>(dR123, wt4);
                        unmatched_gen->template fill<5>(dR123, wt5);
                        unmatched_gen->template fill<6>(dR123, wt6);
                    }
                }

                if constexpr(doTransfer){
                    proj_transferloop<3>(
                        *transfer,
                        thisjet_reco,
                        {n1, n2, n3},
                        dR123,
                        {wt3, wt4, wt5, wt6}
                    );
                }

                for(unsigned i4=i3+1; i4<thisjet_gen->N; ++i4){
                    const double& E4 = thisjet_gen->singles.getE(i4);
                    [[maybe_unused]] bool matched4;
                    if constexpr(doUnmatched){
                        matched4 = matched3 && matched->at(i4);
                    }
                    [[maybe_unused]] EEC::neighborhood const * n4 = nullptr;
                    if constexpr(doTransfer){
                        n4 = &(adj->get_neighborhood(i4));
                    }

                    const typename ResultType::T& dR14 = thisjet_gen->pairs.getDR(i1, i4);  
                    const typename ResultType::T& dR24 = thisjet_gen->pairs.getDR(i2, i4);
                    const typename ResultType::T& dR34 = thisjet_gen->pairs.getDR(i3, i4);

                    const double E1234 = E123*E4;
                    const typename ResultType::T& dR1234 = std::max({
                        dR123,
                        dR14,
                        dR24,
                        dR34
                    });

                    /*
                     * four-particle part
                     */
                    wt4 = 24*E1234;
                    wt5 = 60*E1234*(E1+E2+E3+E4);
                    wt6 = 60*E1234*(2*(E1*E1+E2*E2+E3*E3+E4*E4)+3*(E1*E2+E1*E3+E1*E4+E2*E3+E2*E4+E3*E4));

                    result.template fill<4>(dR1234, wt4);
                    result.template fill<5>(dR1234, wt5);
                    result.template fill<6>(dR1234, wt6);

                    if constexpr (doUnmatched){
                        if (!matched4){
                            unmatched_gen->template fill<4>(dR1234, wt4);
                            unmatched_gen->template fill<5>(dR1234, wt5);
                            unmatched_gen->template fill<6>(dR1234, wt6);
                        }
                    }

                    if constexpr(doTransfer){
                        proj_transferloop<4>(
                            *transfer,
                            thisjet_reco,
                            {n1, n2, n3, n4},
                            dR1234,
                            {wt4, wt5, wt6}
                        );
                    }

                    for(unsigned i5=i4+1; i5<thisjet_gen->N; ++i5){
                        const double& E5 = thisjet_gen->singles.getE(i5);
                        [[maybe_unused]] bool matched5;
                        if constexpr(doUnmatched){
                            matched5 = matched4 && matched->at(i5);
                        }
                        [[maybe_unused]] EEC::neighborhood const * n5 = nullptr;
                        if constexpr(doTransfer){
                            n5 = &(adj->get_neighborhood(i5));
                        }

                        const typename ResultType::T& dR15 = thisjet_gen->pairs.getDR(i1, i5);
                        const typename ResultType::T& dR25 = thisjet_gen->pairs.getDR(i2, i5);
                        const typename ResultType::T& dR35 = thisjet_gen->pairs.getDR(i3, i5);
                        const typename ResultType::T& dR45 = thisjet_gen->pairs.getDR(i4, i5);

                        const double E12345 = E1234*E5;
                        const typename ResultType::T& dR12345 = std::max({
                            dR1234,
                            dR15,
                            dR25,
                            dR35,
                            dR45
                        });

                        /*
                         * five-particle part
                         */

                        wt5 = 120*E12345;
                        wt6 = 360*E12345*(E1+E2+E3+E4+E5);

                        result.template fill<5>(dR12345, wt5);
                        result.template fill<6>(dR12345, wt6);

                        if constexpr (doUnmatched){
                            if (!matched5){
                                unmatched_gen->template fill<5>(dR12345, wt5);
                                unmatched_gen->template fill<6>(dR12345, wt6);
                            }
                        }

                        if constexpr(doTransfer){
                            proj_transferloop<5>(
                                *transfer,
                                thisjet_reco,
                                {n1, n2, n3, n4, n5},
                                dR12345,
                                {wt5, wt6}
                            );
                        }

                        for(unsigned i6=i5+1; i6<thisjet_gen->N; ++i6){
                            const double& E6 = thisjet_gen->singles.getE(i6);
                            [[maybe_unused]] bool matched6;
                            if constexpr(doUnmatched){
                                matched6 = matched5 && matched->at(i6);
                            }
                            [[maybe_unused]] EEC::neighborhood const * n6 = nullptr;
                            if constexpr(doTransfer){
                                n6 = &(adj->get_neighborhood(i6));
                            }

                            const typename ResultType::T& dR16 = thisjet_gen->pairs.getDR(i1, i6);
                            const typename ResultType::T& dR26 = thisjet_gen->pairs.getDR(i2, i6);
                            const typename ResultType::T& dR36 = thisjet_gen->pairs.getDR(i3, i6);
                            const typename ResultType::T& dR46 = thisjet_gen->pairs.getDR(i4, i6);
                            const typename ResultType::T& dR56 = thisjet_gen->pairs.getDR(i5, i6);

                            const double E123456 = E12345*E6;
                            const typename ResultType::T& dR123456 = std::max({
                                dR12345,
                                dR16,
                                dR26,
                                dR36,
                                dR46,
                                dR56
                            });

                            /*
                             * six-particle part
                             */
                            wt6 = 720*E123456;
                            result.template fill<6>(dR123456, wt6);
                            if constexpr (doUnmatched){
                                if (!matched6){
                                    unmatched_gen->template fill<6>(dR123456, wt6);
                                }
                            }

                            if constexpr(doTransfer){
                                proj_transferloop<6>(
                                    *transfer,
                                    thisjet_reco,
                                    {n1, n2, n3, n4, n5, n6},
                                    dR123456,
                                    {wt6}
                                );
                            }
                        }
                    }
                }
            }
        }
    }
}

#endif
