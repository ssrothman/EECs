#include "CARes3Calculator.h"
#include "Adjacency.h"
#include "Res3Result.h"
#include "Res3TransferResult.h"
#include "Res3Backend.h"
#include "Res3Axes.h"

#include "SRothman/SimonTools/src/deltaR.h"
#include "SRothman/SimonTools/src/histutil.h"
#include "SRothman/SimonTools/src/util.h"

#include <memory>

#include <array>

template <class BasicContainer, class JetType>
static void call_res3(
        EEC::Res3Result<BasicContainer>& result,

        const simon::jet& J,
        const EEC::normType& nt,

        const EEC::Res3Axes& axes) noexcept {

    auto thisjet = std::make_shared<JetType>(
            J, nt);
    result.set_pt_denom(thisjet->singles.get_pt_denom());

    res3_mainloop<true, BasicContainer, JetType, false, BasicContainer, false>(
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

void EEC::CARes3Calculator::compute_JIT(
        const simon::jet& J,
        Res3Result_Vector& result) const noexcept {
    
    call_res3<ResVectorContainer, EECjet_JIT>(
            result, 
            J, nt, 
            axes);
}

void EEC::CARes3Calculator::compute_JIT(
        const simon::jet& J,
        Res3Result_MultiArray& result) const noexcept {

    call_res3<ResMultiArrayContainer, EECjet_JIT>(
            result, 
            J, nt, 
            axes);
}

void EEC::CARes3Calculator::compute_precomputed(
        const simon::jet& J,
        Res3Result_Vector& result) const noexcept {
    
    call_res3<ResVectorContainer, EECjet_Precomputed>(
            result, 
            J, nt, 
            axes);
}

void EEC::CARes3Calculator::compute_precomputed(
        const simon::jet& J,
        Res3Result_MultiArray& result) const noexcept {

    call_res3<ResMultiArrayContainer, EECjet_Precomputed>(
            result, 
            J, nt, 
            axes);
}

template <class BasicContainer, class JetType>
static void call_res3_matched(
        EEC::Res3Result<BasicContainer>& result,
        EEC::Res3Result<BasicContainer>& unmatched,

        const simon::jet& J,
        const std::vector<bool>& matched,
        const EEC::normType& nt,

        const EEC::Res3Axes& axes) noexcept {

    auto thisjet = std::make_shared<const JetType>(J, nt);
    result.set_pt_denom(thisjet->singles.get_pt_denom());
    unmatched.set_pt_denom(thisjet->singles.get_pt_denom());

    res3_mainloop<true, BasicContainer, JetType, true, BasicContainer, false>(
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

void EEC::CARes3Calculator::compute_JIT_matched(
        const simon::jet& J,
        const std::vector<bool>& matched,
        Res3Result_Vector& result,
        Res3Result_Vector& unmatched) const noexcept {

    call_res3_matched<ResVectorContainer, EECjet_JIT>(
            result, 
            unmatched,
            J, matched, nt, 
            axes);
}

void EEC::CARes3Calculator::compute_JIT_matched(
        const simon::jet& J,
        const std::vector<bool>& matched,
        Res3Result_MultiArray& result,
        Res3Result_MultiArray& unmatched) const noexcept {

    call_res3_matched<ResMultiArrayContainer, EECjet_JIT>(
            result, 
            unmatched,
            J, matched, nt, 
            axes); 
}

void EEC::CARes3Calculator::compute_precomputed_matched(
        const simon::jet& J,
        const std::vector<bool>& matched,
        Res3Result_Vector& result,
        Res3Result_Vector& unmatched) const noexcept {

    call_res3_matched<ResVectorContainer, EECjet_Precomputed>(
            result, 
            unmatched,
            J, matched, nt, 
            axes);
}

void EEC::CARes3Calculator::compute_precomputed_matched(
        const simon::jet& J,
        const std::vector<bool>& matched,
        Res3Result_MultiArray& result,
        Res3Result_MultiArray& unmatched) const noexcept {

    call_res3_matched<ResMultiArrayContainer, EECjet_Precomputed>(
            result, 
            unmatched,
            J, matched, nt, 
            axes);
}

template <class JetType, class TransferContainer, class BasicContainer>
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

    auto thisjet_reco = std::make_shared<const JetType>(J_reco, nt);
    auto thisjet_gen = std::make_shared<const JetType>(J_gen, nt);

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

    res3_mainloop<true, BasicContainer, JetType, true, TransferContainer, true>(
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

void EEC::CARes3TransferCalculator::compute_JIT(
        const simon::jet& J_reco,
        const simon::jet& J_gen, 
        const Eigen::MatrixXd& tmat,
        Res3Result_Vector& result,
        EEC::Res3Result_Vector& unmatched_gen,
        Res3TransferResult_Vector& tresult) const noexcept {

    call_res3_transfer<EECjet_JIT>(
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

void EEC::CARes3TransferCalculator::compute_JIT(
        const simon::jet& J_reco,
        const simon::jet& J_gen, 
        const Eigen::MatrixXd& tmat,
        Res3Result_MultiArray& result,
        EEC::Res3Result_MultiArray& unmatched_gen,
        Res3TransferResult_Vector& tresult) const noexcept {

    call_res3_transfer<EECjet_JIT>(
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

void EEC::CARes3TransferCalculator::compute_JIT(
        const simon::jet& J_reco,
        const simon::jet& J_gen, 
        const Eigen::MatrixXd& tmat,
        Res3Result_Vector& result,
        EEC::Res3Result_Vector& unmatched_gen,
        Res3TransferResult_MultiArray& tresult) const noexcept {

    call_res3_transfer<EECjet_JIT>(
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

void EEC::CARes3TransferCalculator::compute_JIT(
        const simon::jet& J_reco,
        const simon::jet& J_gen, 
        const Eigen::MatrixXd& tmat,
        Res3Result_MultiArray& result,
        EEC::Res3Result_MultiArray& unmatched_gen,
        Res3TransferResult_MultiArray& tresult) const noexcept {

    call_res3_transfer<EECjet_JIT>(
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

void EEC::CARes3TransferCalculator::compute_precomputed(
        const simon::jet& J_reco,
        const simon::jet& J_gen, 
        const Eigen::MatrixXd& tmat,
        Res3Result_Vector& result,
        EEC::Res3Result_Vector& unmatched_gen,
        Res3TransferResult_Vector& tresult) const noexcept {

    call_res3_transfer<EECjet_Precomputed>(
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

void EEC::CARes3TransferCalculator::compute_precomputed(
        const simon::jet& J_reco,
        const simon::jet& J_gen, 
        const Eigen::MatrixXd& tmat,
        Res3Result_MultiArray& result,
        EEC::Res3Result_MultiArray& unmatched_gen,
        Res3TransferResult_Vector& tresult) const noexcept {

    call_res3_transfer<EECjet_Precomputed>(
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

void EEC::CARes3TransferCalculator::compute_precomputed(
        const simon::jet& J_reco,
        const simon::jet& J_gen, 
        const Eigen::MatrixXd& tmat,
        Res3Result_Vector& result,
        EEC::Res3Result_Vector& unmatched_gen,
        Res3TransferResult_MultiArray& tresult) const noexcept {

    call_res3_transfer<EECjet_Precomputed>(
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

void EEC::CARes3TransferCalculator::compute_precomputed(
        const simon::jet& J_reco,
        const simon::jet& J_gen, 
        const Eigen::MatrixXd& tmat,
        Res3Result_MultiArray& result,
        EEC::Res3Result_MultiArray& unmatched_gen,
        Res3TransferResult_MultiArray& tresult) const noexcept {

    call_res3_transfer<EECjet_Precomputed>(
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

#ifdef CMSSW_GIT_HASH

EEC::CARes3Calculator::CARes3Calculator(const edm::ParameterSet& iConfig) :
    axes(iConfig.getParameter<edm::ParameterSet>("bins")),
    nt(normType_from_string(iConfig.getParameter<std::string>("normType"))) {}

void EEC::CARes3Calculator::fillPSetDescription(edm::ParameterSetDescription& desc) {
    edm::ParameterSetDescription binsDesc; 
    Res3Axes::fillPSetDescription(binsDesc);
    desc.add<edm::ParameterSetDescription>("bins", binsDesc);

    desc.add<std::string>("normType");
}

EEC::CARes3TransferCalculator::CARes3TransferCalculator(const edm::ParameterSet& iConfig):
    axes_reco(iConfig.getParameter<edm::ParameterSet>("bins_reco")),
    axes_gen(iConfig.getParameter<edm::ParameterSet>("bins_gen")),
    nt(normType_from_string(iConfig.getParameter<std::string>("normType"))) {}

void EEC::CARes3TransferCalculator::fillPSetDescription(edm::ParameterSetDescription& desc) {

    edm::ParameterSetDescription binsRecoDesc;
    Res3Axes::fillPSetDescription(binsRecoDesc);
    desc.add<edm::ParameterSetDescription>("bins_reco", binsRecoDesc);

    edm::ParameterSetDescription binsGenDesc;
    Res3Axes::fillPSetDescription(binsGenDesc);
    desc.add<edm::ParameterSetDescription>("bins_gen", binsGenDesc);

    desc.add<std::string>("normType");
}

#endif
