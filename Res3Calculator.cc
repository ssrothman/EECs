#include "Res3Calculator.h"
#include "Adjacency.h"
#include "Res3Result.h"
#include "Res3TransferResult.h"

#include "SRothman/SimonTools/src/deltaR.h"
#include "SRothman/SimonTools/src/histutil.h"
#include "SRothman/SimonTools/src/util.h"

#include <memory>

#include <array>

template <bool distances_squared>
static void get_indices(const double dR12, const double dR13, const double dR23,
                        unsigned& R_idx, unsigned& r_idx, unsigned& c_idx,
                        const EEC::Res3Axes& axes) noexcept {

    std::array<double, 3> dRs = {dR12, dR13, dR23};
    if constexpr(distances_squared){
        for(auto& dR : dRs){
            dR = std::sqrt(dR);
        }
    }
    std::sort(dRs.begin(), dRs.end());
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

template <class PairsType, class TransferContainer>
static void res3_transferloop(
        EEC::Res3TransferResult<TransferContainer>& transfer,

        const std::shared_ptr<const EEC::EECjet<PairsType>> thisjet_reco,

        const EEC::Res3Axes& axes_reco,

        const unsigned R_idx_gen,
        const unsigned r_idx_gen,
        const unsigned c_idx_gen,

        const EEC::neighborhood& n1,
        const EEC::neighborhood& n2,
        const EEC::neighborhood& n3,

        const double wt_gen) noexcept {

    for(const EEC::neighbor& j1 : n1){
        const double twt1 = wt_gen * j1.wt;

        for(const EEC::neighbor& j2 : n2){
            const double twt2 = twt1 * j2.wt;

            const auto&[deta12, dphi12, dR12] = thisjet_reco->pairs.get(j1.idx, j2.idx);

            for(const EEC::neighbor& j3 : n3){
                const double twt3 = twt2 * j3.wt;

                const auto&[deta13, dphi13, dR13] = thisjet_reco->pairs.get(j1.idx, j3.idx);
                const auto&[deta23, dphi23, dR23] = thisjet_reco->pairs.get(j2.idx, j3.idx);

                unsigned R_idx_reco, r_idx_reco, c_idx_reco;
                get_indices<PairsType::distances_squared>(dR12, dR13, dR23, R_idx_reco, r_idx_reco, c_idx_reco, axes_reco);

                transfer.fill(R_idx_reco, r_idx_reco, c_idx_reco,
                              R_idx_gen, r_idx_gen, c_idx_gen, 
                              twt3, wt_gen);
            }
        }
    }
}

template <class BasicContainer, class PairsType, bool doUnmatched, class TransferContainer, bool doTransfer>
static void res3_mainloop(
        EEC::Res3Result<BasicContainer>& result,
        [[maybe_unused]] EEC::Res3Result<BasicContainer>* unmatched_gen,
        [[maybe_unused]] EEC::Res3TransferResult<TransferContainer>* transfer,

        [[maybe_unused]] const std::shared_ptr<const EEC::EECjet<PairsType>> thisjet_reco,
        const std::shared_ptr<const EEC::EECjet<PairsType>> thisjet_gen,
        [[maybe_unused]] const std::vector<bool>* matched,

        [[maybe_unused]] const std::shared_ptr<const EEC::Adjacency> adj,

        [[maybe_unused]] const EEC::Res3Axes* axes_reco,
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
        result.fill(0, 0, 0, wt);
        if constexpr(doUnmatched){
            if(!matched1){
                unmatched_gen->fill(0, 0, 0, wt);
            }
        }
        if constexpr(doTransfer){
            res3_transferloop(
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
            if constexpr(PairsType::distances_squared){
                R_idx = simon::getIndex(std::sqrt(dR12), axes_gen.R);
            } else {
                R_idx = simon::getIndex(dR12, axes_gen.R);
            }
            wt = 3 * E12 * (E1 + E2);
            result.fill(R_idx, 0, 0, wt);
            if constexpr(doUnmatched){
                if(!matched2){
                    unmatched_gen->fill(R_idx, 0, 0, wt);
                }
            }
            if constexpr(doTransfer){
                double wt1 = 3*E12*E1;
                res3_transferloop(
                        *transfer,
                        thisjet_reco,
                        *axes_reco,
                        R_idx, 0, 0,
                        *n1,
                        *n1,
                        *n2,
                        wt1);
                double wt2 = 3*E12*E2;
                res3_transferloop(
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
                get_indices<PairsType::distances_squared>(dR12, dR13, dR23, R_idx, r_idx, c_idx, axes_gen);
                
                wt = E123*6;
                result.fill(R_idx, r_idx, c_idx, wt);
                if constexpr(doUnmatched){
                    if(!matched3){
                        unmatched_gen->fill(R_idx, r_idx, c_idx, wt);
                    }
                }
                if constexpr(doTransfer){
                    res3_transferloop(
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

template <class BasicContainer, class PairsType>
static void call_res3(
        EEC::Res3Result<BasicContainer>& result,

        const simon::jet& J,
        const EEC::normType& nt,

        const EEC::Res3Axes& axes) noexcept {

    auto thisjet = std::make_shared<EEC::EECjet<PairsType>>(
            J, nt);
    result.set_pt_denom(thisjet->singles.get_pt_denom());

    res3_mainloop<BasicContainer, PairsType, false, BasicContainer, false>(
            result,
            nullptr,
            nullptr, 
            nullptr,
            thisjet, 
            nullptr,
            nullptr,
            nullptr,
            axes);
}

void EEC::Res3Calculator::compute_JIT(
        const simon::jet& J,
        Res3Result_Vector& result) const noexcept {
    
    call_res3<Res3VectorContainer, JITPairs>(
            result, 
            J, nt, 
            axes);
}

void EEC::Res3Calculator::compute_JIT(
        const simon::jet& J,
        Res3Result_MultiArray& result) const noexcept {

    call_res3<Res3MultiArrayContainer, JITPairs>(
            result, 
            J, nt, 
            axes);
}

void EEC::Res3Calculator::compute_precomputed(
        const simon::jet& J,
        Res3Result_Vector& result) const noexcept {
    
    call_res3<Res3VectorContainer, PrecomputedPairs>(
            result, 
            J, nt, 
            axes);
}

void EEC::Res3Calculator::compute_precomputed(
        const simon::jet& J,
        Res3Result_MultiArray& result) const noexcept {

    call_res3<Res3MultiArrayContainer, PrecomputedPairs>(
            result, 
            J, nt, 
            axes);
}

template <class BasicContainer, class PairsType>
static void call_res3_matched(
        EEC::Res3Result<BasicContainer>& result,
        EEC::Res3Result<BasicContainer>& unmatched,

        const simon::jet& J,
        const std::vector<bool>& matched,
        const EEC::normType& nt,

        const EEC::Res3Axes& axes) noexcept {

    auto thisjet = std::make_shared<const EEC::EECjet<PairsType>>(J, nt);
    result.set_pt_denom(thisjet->singles.get_pt_denom());
    unmatched.set_pt_denom(thisjet->singles.get_pt_denom());

    res3_mainloop<BasicContainer, PairsType, true, BasicContainer, false>(
            result, 
            &unmatched,
            nullptr, 
            nullptr,
            thisjet,
            &matched,
            nullptr,
            nullptr,
            axes);
}

void EEC::Res3Calculator::compute_JIT_matched(
        const simon::jet& J,
        const std::vector<bool>& matched,
        Res3Result_Vector& result,
        Res3Result_Vector& unmatched) const noexcept {

    call_res3_matched<Res3VectorContainer, JITPairs>(
            result, 
            unmatched,
            J, matched, nt, 
            axes);
}

void EEC::Res3Calculator::compute_JIT_matched(
        const simon::jet& J,
        const std::vector<bool>& matched,
        Res3Result_MultiArray& result,
        Res3Result_MultiArray& unmatched) const noexcept {

    call_res3_matched<Res3MultiArrayContainer, JITPairs>(
            result, 
            unmatched,
            J, matched, nt, 
            axes); 
}

void EEC::Res3Calculator::compute_precomputed_matched(
        const simon::jet& J,
        const std::vector<bool>& matched,
        Res3Result_Vector& result,
        Res3Result_Vector& unmatched) const noexcept {

    call_res3_matched<Res3VectorContainer, PrecomputedPairs>(
            result, 
            unmatched,
            J, matched, nt, 
            axes);
}

void EEC::Res3Calculator::compute_precomputed_matched(
        const simon::jet& J,
        const std::vector<bool>& matched,
        Res3Result_MultiArray& result,
        Res3Result_MultiArray& unmatched) const noexcept {

    call_res3_matched<Res3MultiArrayContainer, PrecomputedPairs>(
            result, 
            unmatched,
            J, matched, nt, 
            axes);
}

template <class PairsType, class TransferContainer, class BasicContainer>
static void call_res3_transfer(
        EEC::Res3Result<BasicContainer>& result_gen,
        EEC::Res3Result<BasicContainer>& unmatched_gen,
        EEC::Res3TransferResult<TransferContainer>& transfer,

        const simon::jet& J_reco,
        const simon::jet& J_gen,
        const EEC::normType& nt,

        const Eigen::MatrixXd& adjmat,

        const EEC::Res3Axes& axes_reco,
        const EEC::Res3Axes& axes_gen) noexcept {

    auto thisjet_reco = std::make_shared<const EEC::EECjet<PairsType>>(J_reco, nt);
    auto thisjet_gen = std::make_shared<const EEC::EECjet<PairsType>>(J_gen, nt);

    result_gen.set_pt_denom(thisjet_gen->singles.get_pt_denom());
    transfer.set_pt_denom(thisjet_reco->singles.get_pt_denom(),
                          thisjet_gen->singles.get_pt_denom());

    double denom_reco = thisjet_reco->singles.get_pt_denom();
    double denom_gen = thisjet_gen->singles.get_pt_denom();
    Eigen::VectorXd ptgen = J_gen.ptvec()/denom_gen;
    Eigen::VectorXd ptreco = J_reco.ptvec()/denom_reco;
    Eigen::VectorXd ptfwd = adjmat*ptgen;

    Eigen::MatrixXd adjmat_copy(adjmat);
    for(unsigned iRecoPart=0; iRecoPart<J_reco.nPart; ++iRecoPart){
        if(ptfwd(iRecoPart) == 0){
            continue;
        }
        double factor = ptreco(iRecoPart) / ptfwd(iRecoPart);
        for(unsigned iGenPart=0; iGenPart<J_gen.nPart; ++iGenPart){
            if(adjmat(iRecoPart, iGenPart) > 0){
                adjmat_copy(iRecoPart, iGenPart) *= factor;
            }
        }
    }

    auto adj = std::make_shared<const EEC::Adjacency>(adjmat_copy);

    std::vector<bool> matched(J_gen.nPart, false);
    for(unsigned iGenPart=0; iGenPart<J_gen.nPart; ++iGenPart){
        for(unsigned iRecoPart=0; iRecoPart<J_reco.nPart; ++iRecoPart){
            if(adjmat_copy(iRecoPart, iGenPart) > 0){
                matched[iGenPart] = true;
                break;
            }
        }
    }

    res3_mainloop<BasicContainer, PairsType, true, TransferContainer, true>(
            result_gen,
            &unmatched_gen,
            &transfer,
            thisjet_reco,
            thisjet_gen,
            &matched,
            adj,
            &axes_reco,
            axes_gen);
}

void EEC::Res3TransferCalculator::compute_JIT(
        const simon::jet& J_reco,
        const simon::jet& J_gen, 
        const Eigen::MatrixXd& tmat,
        Res3Result_Vector& result,
        EEC::Res3Result_Vector& unmatched_gen,
        Res3TransferResult_Vector& tresult) const noexcept {

    call_res3_transfer<JITPairs>(
        result,
        unmatched_gen,
        tresult,

        J_reco,
        J_gen,
        nt,

        tmat,

        axes_reco,
        axes_gen);
}

void EEC::Res3TransferCalculator::compute_JIT(
        const simon::jet& J_reco,
        const simon::jet& J_gen, 
        const Eigen::MatrixXd& tmat,
        Res3Result_MultiArray& result,
        EEC::Res3Result_MultiArray& unmatched_gen,
        Res3TransferResult_Vector& tresult) const noexcept {

    call_res3_transfer<JITPairs>(
        result,
        unmatched_gen,
        tresult,

        J_reco,
        J_gen,
        nt,

        tmat,

        axes_reco,
        axes_gen);
}

void EEC::Res3TransferCalculator::compute_JIT(
        const simon::jet& J_reco,
        const simon::jet& J_gen, 
        const Eigen::MatrixXd& tmat,
        Res3Result_Vector& result,
        EEC::Res3Result_Vector& unmatched_gen,
        Res3TransferResult_MultiArray& tresult) const noexcept {

    call_res3_transfer<JITPairs>(
        result,
        unmatched_gen,
        tresult,

        J_reco,
        J_gen,
        nt,

        tmat,

        axes_reco,
        axes_gen);
}

void EEC::Res3TransferCalculator::compute_JIT(
        const simon::jet& J_reco,
        const simon::jet& J_gen, 
        const Eigen::MatrixXd& tmat,
        Res3Result_MultiArray& result,
        EEC::Res3Result_MultiArray& unmatched_gen,
        Res3TransferResult_MultiArray& tresult) const noexcept {

    call_res3_transfer<JITPairs>(
        result,
        unmatched_gen,
        tresult,

        J_reco,
        J_gen,
        nt,

        tmat,

        axes_reco,
        axes_gen);
}

void EEC::Res3TransferCalculator::compute_precomputed(
        const simon::jet& J_reco,
        const simon::jet& J_gen, 
        const Eigen::MatrixXd& tmat,
        Res3Result_Vector& result,
        EEC::Res3Result_Vector& unmatched_gen,
        Res3TransferResult_Vector& tresult) const noexcept {

    call_res3_transfer<PrecomputedPairs>(
        result,
        unmatched_gen,
        tresult,

        J_reco,
        J_gen,
        nt,

        tmat,

        axes_reco,
        axes_gen);
}

void EEC::Res3TransferCalculator::compute_precomputed(
        const simon::jet& J_reco,
        const simon::jet& J_gen, 
        const Eigen::MatrixXd& tmat,
        Res3Result_MultiArray& result,
        EEC::Res3Result_MultiArray& unmatched_gen,
        Res3TransferResult_Vector& tresult) const noexcept {

    call_res3_transfer<PrecomputedPairs>(
        result,
        unmatched_gen,
        tresult,

        J_reco,
        J_gen,
        nt,

        tmat,

        axes_reco,
        axes_gen);
}

void EEC::Res3TransferCalculator::compute_precomputed(
        const simon::jet& J_reco,
        const simon::jet& J_gen, 
        const Eigen::MatrixXd& tmat,
        Res3Result_Vector& result,
        EEC::Res3Result_Vector& unmatched_gen,
        Res3TransferResult_MultiArray& tresult) const noexcept {

    call_res3_transfer<PrecomputedPairs>(
        result,
        unmatched_gen,
        tresult,

        J_reco,
        J_gen,
        nt,

        tmat,

        axes_reco,
        axes_gen);
}

void EEC::Res3TransferCalculator::compute_precomputed(
        const simon::jet& J_reco,
        const simon::jet& J_gen, 
        const Eigen::MatrixXd& tmat,
        Res3Result_MultiArray& result,
        EEC::Res3Result_MultiArray& unmatched_gen,
        Res3TransferResult_MultiArray& tresult) const noexcept {

    call_res3_transfer<PrecomputedPairs>(
        result,
        unmatched_gen,
        tresult,

        J_reco,
        J_gen,
        nt,

        tmat,

        axes_reco,
        axes_gen);
}


