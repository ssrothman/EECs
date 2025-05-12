#ifndef SROTHMN_EEC_PROJ_BACKEND_H
#define SROTHMN_EEC_PROJ_BACKEND_H

#include "Adjacency.h"

#include "SRothman/SimonTools/src/histutil.h"
#include "SRothman/SimonTools/src/util.h"

#include <memory>

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

        [[maybe_unused]] unsigned j1=0;
        [[maybe_unused]] double Ej1=0;
        if constexpr(doTransfer){
            const auto& n1 = adj->get_neighborhood(i1);
            if constexpr(!doUnmatched){
                matched1 = !(n1.empty());
            }
            if (matched1){
                j1 = n1[0].idx;
                Ej1 = thisjet_reco->singles.getE(j1);

                double jwt2 = Ej1*Ej1;
                double jwt3 = Ej1*jwt2;
                double jwt4 = Ej1*jwt3;
                double jwt5 = Ej1*jwt4;
                double jwt6 = Ej1*jwt5;

                transfer->template fill<2>(0, 0, jwt2, wt2);
                transfer->template fill<3>(0, 0, jwt3, wt3);
                transfer->template fill<4>(0, 0, jwt4, wt4);
                transfer->template fill<5>(0, 0, jwt5, wt5);
                transfer->template fill<6>(0, 0, jwt6, wt6);
            }
        }

        for(unsigned i2=i1+1; i2<thisjet_gen->N; ++i2){
            const double& E2 = thisjet_gen->singles.getE(i2);
            [[maybe_unused]] bool matched2;
            if constexpr(doUnmatched){
                matched2 = matched1 && matched->at(i2);
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

            [[maybe_unused]] unsigned j2=0;
            [[maybe_unused]] double Ej2=0, Ej12=0;
            [[maybe_unused]] typename ResultType::T dRj12=0;
            if constexpr(doTransfer){
                const auto& n2 = adj->get_neighborhood(i2);
                if constexpr(!doUnmatched){
                    matched2 = matched1 && !(n2.empty());
                }
                if (matched2){
                    j2 = n2[0].idx;
                    Ej2 = thisjet_reco->singles.getE(j2);

                    Ej12 = Ej1 * Ej2;

                    dRj12 = thisjet_reco->pairs.getDR(j1, j2);

                    double jwt2 = 2*Ej12;
                    double jwt3 = 3*Ej12*(Ej1+Ej2);
                    double jwt4 = 2*Ej12*(2*(Ej1*Ej1+Ej2*Ej2) + 3*Ej12);
                    double jwt5 = 5*Ej12*(Ej1+Ej2)*(Ej1*Ej1 + Ej1*Ej2 + Ej2*Ej2);
                    double jwt6 = Ej12*(6*(Ej1*Ej1*Ej1*Ej1+Ej2*Ej2*Ej2*Ej2) + 15*(Ej1*Ej1*Ej1*Ej2+Ej1*Ej2*Ej2*Ej2) + 20*Ej1*Ej1*Ej2*Ej2);   

                    transfer->template fill<2>(dRj12, dR12, jwt2, wt2);
                    transfer->template fill<3>(dRj12, dR12, jwt3, wt3);
                    transfer->template fill<4>(dRj12, dR12, jwt4, wt4);
                    transfer->template fill<5>(dRj12, dR12, jwt5, wt5);
                    transfer->template fill<6>(dRj12, dR12, jwt6, wt6);
                }
            }

            for(unsigned i3=i2+1; i3<thisjet_gen->N; ++i3){
                const double& E3 = thisjet_gen->singles.getE(i3);
                [[maybe_unused]] bool matched3;
                if constexpr(doUnmatched){
                    matched3 = matched2 && matched->at(i3);
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
                [[maybe_unused]] unsigned j3=0;
                [[maybe_unused]] double Ej3=0, Ej123=0;
                [[maybe_unused]] typename ResultType::T dRj123=0;
                if constexpr(doTransfer){
                    const auto& n3 = adj->get_neighborhood(i3);
                    if constexpr(!doUnmatched){
                        matched3 = matched2 && !(n3.empty());
                    }
                    if (matched3){
                        j3 = n3[0].idx;
                        Ej3 = thisjet_reco->singles.getE(j3);

                        Ej123 = Ej12 * Ej3;

                        const typename ResultType::T& dRj13 = thisjet_reco->pairs.getDR(j1, j3);
                        const typename ResultType::T& dRj23 = thisjet_reco->pairs.getDR(j2, j3);
                        dRj123 = std::max({
                            dRj12,
                            dRj13,
                            dRj23
                        });

                        double jwt3 = 6*Ej123;
                        double jwt4 = 12*Ej123*(Ej1+Ej2+Ej3);
                        double jwt5 = 10*Ej123*(2*(Ej1*Ej1+Ej2*Ej2+Ej3*Ej3) + 3*(Ej1*Ej2+Ej1*Ej3+Ej2*Ej3));
                        double jwt6 = 30*Ej123*(Ej1+Ej2+Ej3)*(Ej1*Ej1+Ej2*Ej2+Ej3*Ej3+Ej1*Ej2+Ej1*Ej3+Ej2*Ej3);

                        transfer->template fill<3>(dRj123, dR123, jwt3, wt3);
                        transfer->template fill<4>(dRj123, dR123, jwt4, wt4);
                        transfer->template fill<5>(dRj123, dR123, jwt5, wt5);
                        transfer->template fill<6>(dRj123, dR123, jwt6, wt6);
                    }
                }

                for(unsigned i4=i3+1; i4<thisjet_gen->N; ++i4){
                    const double& E4 = thisjet_gen->singles.getE(i4);
                    [[maybe_unused]] bool matched4;
                    if constexpr(doUnmatched){
                        matched4 = matched3 && matched->at(i4);
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

                    [[maybe_unused]] unsigned j4=0;
                    [[maybe_unused]] double Ej4=0, Ej1234=0;
                    [[maybe_unused]] typename ResultType::T dRj1234=0;
                    if constexpr(doTransfer){
                        const auto& n4 = adj->get_neighborhood(i4);
                        if constexpr(!doUnmatched){
                            matched4 = matched3 && !(n4.empty());
                        }
                        if (matched4){
                            j4 = n4[0].idx;
                            Ej4 = thisjet_reco->singles.getE(j4);

                            Ej1234 = Ej123 * Ej4;
                            const typename ResultType::T& dRj14 = thisjet_reco->pairs.getDR(j1, j4);
                            const typename ResultType::T& dRj24 = thisjet_reco->pairs.getDR(j2, j4);
                            const typename ResultType::T& dRj34 = thisjet_reco->pairs.getDR(j3, j4);

                            dRj1234 = std::max({
                                dRj123,
                                dRj14,
                                dRj24,
                                dRj34
                            });

                            double jwt4 = 24*Ej1234;
                            double jwt5 = 60*Ej1234*(Ej1+Ej2+Ej3+Ej4);
                            double jwt6 = 60*Ej1234*(2*(Ej1*Ej1+Ej2*Ej2+Ej3*Ej3+Ej4*Ej4)+3*(Ej1*Ej2+Ej1*Ej3+Ej1*Ej4+Ej2*Ej3+Ej2*Ej4+Ej3*Ej4));

                            transfer->template fill<4>(dRj1234, dR1234, jwt4, wt4);
                            transfer->template fill<5>(dRj1234, dR1234, jwt5, wt5);
                            transfer->template fill<6>(dRj1234, dR1234, jwt6, wt6);
                        }
                    }

                    for(unsigned i5=i4+1; i5<thisjet_gen->N; ++i5){
                        const double& E5 = thisjet_gen->singles.getE(i5);
                        [[maybe_unused]] bool matched5;
                        if constexpr(doUnmatched){
                            matched5 = matched4 && matched->at(i5);
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

                        [[maybe_unused]] unsigned j5=0;
                        [[maybe_unused]] double Ej5=0, Ej12345=0;
                        [[maybe_unused]] typename ResultType::T dRj12345=0;
                        if constexpr(doTransfer){
                            auto n5 = adj->get_neighborhood(i5);
                            if constexpr(!doUnmatched){
                                matched5 = matched4 && !(n5.empty());
                            }
                            if (matched5){
                                j5 = n5[0].idx;
                                Ej5 = thisjet_reco->singles.getE(j5);
                                Ej12345 = Ej1234 * Ej5;

                                const typename ResultType::T& dRj15 = thisjet_reco->pairs.getDR(j1, j5);
                                const typename ResultType::T& dRj25 = thisjet_reco->pairs.getDR(j2, j5);
                                const typename ResultType::T& dRj35 = thisjet_reco->pairs.getDR(j3, j5);
                                const typename ResultType::T& dRj45 = thisjet_reco->pairs.getDR(j4, j5);

                                dRj12345 = std::max({
                                    dRj1234,
                                    dRj15,
                                    dRj25,
                                    dRj35,
                                    dRj45
                                });

                                double jwt5 = 120*Ej12345;
                                double jwt6 = 360*Ej12345*(Ej1+Ej2+Ej3+Ej4+Ej5);

                                transfer->template fill<5>(dRj12345, dR12345, jwt5, wt5);
                                transfer->template fill<6>(dRj12345, dR12345, jwt6, wt6);
                            }
                        }

                        for(unsigned i6=i5+1; i6<thisjet_gen->N; ++i6){
                            const double& E6 = thisjet_gen->singles.getE(i6);
                            [[maybe_unused]] bool matched6;
                            if constexpr(doUnmatched){
                                matched6 = matched5 && matched->at(i6);
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

                            [[maybe_unused]] unsigned j6=0;
                            [[maybe_unused]] double Ej6=0, Ej123456=0;
                            [[maybe_unused]] typename ResultType::T dRj123456=0;
                            if constexpr(doTransfer){
                                const auto& n6 = adj->get_neighborhood(i6);
                                if constexpr(!doUnmatched){
                                    matched6 = matched5 && !(n6.empty());
                                }
                                if (matched6){
                                    j6 = n6[0].idx;
                                    Ej6 = thisjet_reco->singles.getE(j6);
                                    Ej123456 = Ej12345 * Ej6;

                                    const typename ResultType::T& dRj16 = thisjet_reco->pairs.getDR(j1, j6);
                                    const typename ResultType::T& dRj26 = thisjet_reco->pairs.getDR(j2, j6);
                                    const typename ResultType::T& dRj36 = thisjet_reco->pairs.getDR(j3, j6);
                                    const typename ResultType::T& dRj46 = thisjet_reco->pairs.getDR(j4, j6);
                                    const typename ResultType::T& dRj56 = thisjet_reco->pairs.getDR(j5, j6);

                                    dRj123456 = std::max({
                                        dRj12345,
                                        dRj16,
                                        dRj26,
                                        dRj36,
                                        dRj46,
                                        dRj56
                                    });

                                    double jwt6 = 720*Ej123456;

                                    transfer->template fill<6>(dRj123456, dR123456, jwt6, wt6);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

#endif
