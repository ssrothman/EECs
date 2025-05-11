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

template <class JetType, class ResultType, class TransferResultType>
static void call_transfer(
        ResultType& result,
        ResultType& unmatched,
        TransferResultType& tresult,

        const simon::jet& J_reco,
        const simon::jet& J_gen,
        const Eigen::MatrixXd& tmat,
        const EEC::normType& nt,

        const EEC::ProjAxes& axes_reco,
        const EEC::ProjAxes& axes_gen) noexcept {

    auto thisjet_reco = std::make_shared<JetType>(
        J_reco, nt, axes_reco
    );
    auto thisjet_gen = std::make_shared<JetType>(
        J_gen, nt, axes_gen
    );

    result.set_pt_denom(thisjet_gen->singles.get_pt_denom());
    unmatched.set_pt_denom(thisjet_gen->singles.get_pt_denom());
    tresult.set_pt_denom(thisjet_reco->singles.get_pt_denom(),
                         thisjet_gen->singles.get_pt_denom());

    double denom_reco = thisjet_reco->singles.get_pt_denom();
    double denom_gen = thisjet_gen->singles.get_pt_denom();
    Eigen::VectorXd ptgen = J_gen.ptvec()/denom_gen;
    Eigen::VectorXd ptreco = J_reco.ptvec()/denom_reco;
    Eigen::VectorXd ptfwd = tmat * ptgen;

    Eigen::MatrixXd tmat_copy(tmat);
    for(unsigned iRecoPart=0; iRecoPart<J_reco.nPart; ++iRecoPart){
        if(ptfwd(iRecoPart) == 0){
            continue;
        }
        double factor = ptreco(iRecoPart) / ptfwd(iRecoPart);
        for(unsigned iGenPart=0; iGenPart<J_gen.nPart; ++iGenPart){
            if(tmat(iRecoPart, iGenPart) > 0){
                tmat_copy(iRecoPart, iGenPart) *= factor;
            }
        }
    }

    auto adj = std::make_shared<EEC::Adjacency>(tmat_copy);

    std::vector<bool> matched(J_gen.nPart, false);
    for(unsigned iGenPart=0; iGenPart<J_gen.nPart; ++iGenPart){
        for(unsigned iRecoPart=0; iRecoPart<J_reco.nPart; ++iRecoPart){
            if(tmat_copy(iRecoPart, iGenPart) > 0){
                matched[iGenPart] = true;
                break;
            }
        }
    }
    
    proj_mainloop<ResultType, JetType, true, TransferResultType, true>(
        result,
        &unmatched,
        &tresult,
        thisjet_reco,
        thisjet_gen,
        &matched, 
        adj);
}

void EEC::ProjTransferCalculator::compute(
        const simon::jet& J_reco,
        const simon::jet& J_gen,
        const Eigen::MatrixXd& tmat,
        ProjResult_Vector& result,
        ProjResult_Vector& unmatched,
        ProjTransferResult_Vector& tresult) const noexcept {

    call_transfer<EEC::EECjet_ProjBinned>(
        result, unmatched, tresult,
        J_reco, J_gen, tmat,
        nt, axes_reco, axes_gen
    );
}

void EEC::ProjTransferCalculator::compute(
        const simon::jet& J_reco,
        const simon::jet& J_gen,
        const Eigen::MatrixXd& tmat,
        ProjResult_Unbinned& result,
        ProjResult_Unbinned& unmatched,
        ProjTransferResult_Unbinned& tresult) const noexcept {

    call_transfer<EEC::EECjet_ProjUnbinned>(
        result, unmatched, tresult,
        J_reco, J_gen, tmat,
        nt, axes_reco, axes_gen
    );
}

void EEC::ProjTransferCalculator::compute(
        const simon::jet& J_reco,
        const simon::jet& J_gen,
        const Eigen::MatrixXd& tmat,
        ProjResult_Array& result,
        ProjResult_Array& unmatched,
        ProjTransferResult_Array& tresult) const noexcept {

    call_transfer<EEC::EECjet_ProjBinned>(
        result, unmatched, tresult,
        J_reco, J_gen, tmat,
        nt, axes_reco, axes_gen
    );
}

void EEC::ProjTransferCalculator::compute(
        const simon::jet& J_reco,
        const simon::jet& J_gen,
        const Eigen::MatrixXd& tmat,
        ProjResult_Array& result,
        ProjResult_Array& unmatched,
        ProjTransferResult_Vector& tresult) const noexcept {

    call_transfer<EEC::EECjet_ProjBinned>(
        result, unmatched, tresult,
        J_reco, J_gen, tmat,
        nt, axes_reco, axes_gen
    );
}

void EEC::ProjTransferCalculator::compute(
        const simon::jet& J_reco,
        const simon::jet& J_gen,
        const Eigen::MatrixXd& tmat,
        ProjResult_Vector& result,
        ProjResult_Vector& unmatched,
        ProjTransferResult_Array& tresult) const noexcept {

    call_transfer<EEC::EECjet_ProjBinned>(
        result, unmatched, tresult,
        J_reco, J_gen, tmat,
        nt, axes_reco, axes_gen
    );
}

