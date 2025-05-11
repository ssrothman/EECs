#include "ProjCalculator.h"
#include "Adjacency.h"
#include "ProjResult.h"
#include "ProjTransferResult.h"
#include "ProjBackend.h"

template <class JetType, class ResultType>
static void call_proj(
        ResultType& result,

        const simon::jet& J,
        const EEC::normType& nt,

        const EEC::ProjAxes& axes) noexcept {

    auto thisjet = std::make_shared<JetType>(
        J, nt, axes
    );
    result.set_pt_denom(thisjet->singles.get_pt_denom());

    proj_mainloop<ResultType, JetType, false, ResultType, false>(
        result,
        nullptr, 
        nullptr, 
        nullptr, thisjet, nullptr,
        nullptr
    ); 
}

void EEC::ProjCalculator::compute(
        const simon::jet& J,
        ProjResult_Array& result) const noexcept {

    call_proj<EEC::EECjet_ProjBinned>(
        result,
        J, nt, axes
    );
}

void EEC::ProjCalculator::compute(
        const simon::jet& J,
        ProjResult_Vector& result) const noexcept {

    call_proj<EEC::EECjet_ProjBinned>(
        result, 
        J, nt, axes
    );
}

void EEC::ProjCalculator::compute(
        const simon::jet& J,
        ProjResult_Unbinned& result) const noexcept {

    call_proj<EEC::EECjet_ProjUnbinned>(
        result,
        J, nt, axes
    );
}

template <class JetType, class ResultType>
static void call_proj_matched(
        ResultType& result,
        ResultType& unmatched,

        const simon::jet& J,
        const std::vector<bool>& matched,
        const EEC::normType& nt,

        const EEC::ProjAxes& axes) noexcept {

    auto thisjet = std::make_shared<JetType>(
        J, nt, axes
    );
    result.set_pt_denom(thisjet->singles.get_pt_denom());
    unmatched.set_pt_denom(thisjet->singles.get_pt_denom());

    proj_mainloop<ResultType, JetType, true, ResultType, false>(
        result,
        &unmatched,
        nullptr, 
        nullptr, thisjet, &matched,
        nullptr
    );
}

void EEC::ProjCalculator::compute_matched(
        const simon::jet& J,
        const std::vector<bool>& matched,
        ProjResult_Array& result,
        ProjResult_Array& unmatched) const noexcept {

    call_proj_matched<EEC::EECjet_ProjBinned>(
        result, unmatched,
        J, matched,
        nt, axes
    );
}

void EEC::ProjCalculator::compute_matched(
        const simon::jet& J,
        const std::vector<bool>& matched,
        ProjResult_Vector& result,
        ProjResult_Vector& unmatched) const noexcept {

    call_proj_matched<EEC::EECjet_ProjBinned>(
        result, unmatched,
        J, matched,
        nt, axes
    );
}

void EEC::ProjCalculator::compute_matched(
        const simon::jet& J,
        const std::vector<bool>& matched,
        ProjResult_Unbinned& result,
        ProjResult_Unbinned& unmatched) const noexcept {

    call_proj_matched<EEC::EECjet_ProjUnbinned>(
        result, unmatched,
        J, matched,
        nt, axes
    );
}
