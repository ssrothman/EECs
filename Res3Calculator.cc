#include "Res3Calculator.h"
#include "Adjacency.h"
#include "Res3Result.h"
#include "Res3TransferResult.h"

#include "SRothman/SimonTools/src/deltaR.h"
#include "SRothman/SimonTools/src/histutil.h"
#include "SRothman/SimonTools/src/util.h"

#include <memory>

#include <array>

template <class BasicContainer, class PairsType, bool doUnmatched>
static void res3_mainloop(
        EEC::Res3Result<BasicContainer>& result,
        [[maybe_unused]] EEC::Res3Result<BasicContainer>* unmatched,

        const std::shared_ptr<const EEC::EECjet<PairsType>> thisjet,
        [[maybe_unused]] const std::vector<bool>* matched,

        [[maybe_unused]] const EEC::Res3Axes& axes) noexcept {//TODO: REMOVE [[maybe_unused]] here

    for(unsigned i1=0; i1<thisjet->N; ++i1){
        const auto&[E1, eta1, phi1] = thisjet->singles.get(i1);
        [[maybe_unused]] bool matched1;
        if constexpr(doUnmatched){
            matched1 = matched->at(i1);
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
                unmatched->fill(0, 0, 0, wt);
            }
        }

        for(unsigned i2=i1+1; i2<thisjet->N; ++i2){
            const auto&[E2, eta2, phi2] = thisjet->singles.get(i2);
            [[maybe_unused]] bool matched2;
            if constexpr(doUnmatched){
                matched2 = matched1 && matched->at(i2);
            }

            const double E12 = E1*E2;

            const auto&[deta12, dphi12, dR12] = thisjet->pairs.get(i1, i2);

            /*
             * Two-particle part
             * This has symmetry factor 3!/2! = 3
             * And RL = R12, r=c=0
             */
            int R_idx;
            if constexpr(PairsType::distances_squared){
                R_idx = simon::getIndex(std::sqrt(dR12), axes.R);
            } else {
                R_idx = simon::getIndex(dR12, axes.R);
            }
            wt = 3 * E12 * (E1 + E2);
            result.fill(R_idx, 0, 0, wt);
            if constexpr(doUnmatched){
                if(!matched2){
                    unmatched->fill(R_idx, 0, 0, wt);
                }
            }

            for(unsigned i3=i2+1; i3<thisjet->N; ++i3){
                const auto&[E3, eta3, phi3] = thisjet->singles.get(i3);
                [[maybe_unused]] bool matched3;
                if constexpr(doUnmatched){
                    matched3 = matched2 && matched->at(i3);
                }

                const double E123 = E12*E3;

                const auto&[deta13, dphi13, dR13] = thisjet->pairs.get(i1, i3);
                const auto&[deta23, dphi23, dR23] = thisjet->pairs.get(i2, i3);

                /*
                 * Three-particle part
                 * This has symmetry factor 3! = 6
                 * and full coordinates
                 */
                std::array<double, 3> dRs = {dR12, dR13, dR23};
                if constexpr(PairsType::distances_squared){
                    for(auto& dR : dRs){
                        dR = std::sqrt(dR);
                    }
                }
                std::sort(dRs.begin(), dRs.end());
                double R = dRs[0];
                double r = dRs[2]/dRs[1];
                double c = std::asin(std::sqrt(1 - simon::square((dRs[0] - dRs[1])/dRs[2])));
                R_idx = simon::getIndex(R, axes.R);
                int r_idx = simon::getIndex(r, axes.r);
                int c_idx = simon::getIndex(c, axes.c);
                wt = E123*6;
                result.fill(R_idx, r_idx, c_idx, wt);
                if constexpr(doUnmatched){
                    if(!matched3){
                        unmatched->fill(R_idx, r_idx, c_idx, wt);
                    }
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

    res3_mainloop<BasicContainer, PairsType, false>(
            result,
            nullptr, 
            thisjet, 
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

    res3_mainloop<BasicContainer, PairsType, true>(
            result, 
            &unmatched,
            thisjet,
            &matched,
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
