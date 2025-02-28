#include "Res4Calculator.h"
#include "Adjacency.h"
#include "Res4Result.h"
#include "Res4TransferResult.h"
#include "Res4Backend.h"

#include "SRothman/SimonTools/src/deltaR.h"
#include "SRothman/SimonTools/src/histutil.h"
#include "SRothman/SimonTools/src/util.h"

#include <memory>

#include <array>

template <class ResultType, class JetType>
static void call_res4(
        ResultType& ans,

        const simon::jet& J,
        const EEC::normType& nt,

        const EEC::Res4Axes& axes,
        const double tolerance2,
        const double tri_tolerance) noexcept {

    auto thisjet = std::make_shared<const JetType>(J, nt);
    ans.set_pt_denom(thisjet->singles.get_pt_denom());
    
    res4_mainloop<ResultType, JetType, false, ResultType, false>(
            ans,                                //ans
            nullptr,                            //unmatched_gen
            nullptr,                            //transfer_ans
            nullptr,                            //untransfered_reco
            nullptr,                            //untransfered_gen
            nullptr,                            //thisjet_reco
            thisjet,                            //thisjet_gen
            nullptr,                            //matched
            nullptr,                            //adj
            nullptr,                            //axes_reco
            axes,                               //axes_gen
            tolerance2,                         //tolerance2
            tri_tolerance);                     //tri_tolerance
}

void EEC::Res4Calculator::compute_JIT(
        const simon::jet& J,
        Res4Result_Vector& result) const noexcept {
    
    call_res4<Res4Result_Vector, EECjet_JIT>(
            result, 
            J, nt, 
            axes, 
            tolerance2, 
            tri_tolerance);
}

void EEC::Res4Calculator::compute_JIT(
        const simon::jet& J,
        Res4Result_MultiArray& result) const noexcept {

    call_res4<Res4Result_MultiArray, EECjet_JIT>(
            result, 
            J, nt, 
            axes, 
            tolerance2,
            tri_tolerance);
}

void EEC::Res4Calculator::compute_precomputed(
        const simon::jet& J,
        Res4Result_Vector& result) const noexcept {
    
    call_res4<Res4Result_Vector, EECjet_Precomputed>(
            result, 
            J, nt, 
            axes, 
            tolerance2, 
            tri_tolerance);
}

void EEC::Res4Calculator::compute_precomputed(
        const simon::jet& J,
        Res4Result_MultiArray& result) const noexcept {

    call_res4<Res4Result_MultiArray, EECjet_Precomputed>(
            result, 
            J, nt, 
            axes, 
            tolerance2, 
            tri_tolerance);
}

template <class ResultType, class JetType>
static void call_res4_matched(
        ResultType& ans,
        ResultType& unmatched,

        const simon::jet& J,
        const std::vector<bool>& matched,
        const EEC::normType& nt,

        const EEC::Res4Axes& axes,
        const double tolerance2,
        const double tri_tolerance) noexcept {

    auto thisjet = std::make_shared<const JetType>(J, nt);
    ans.set_pt_denom(thisjet->singles.get_pt_denom());
    unmatched.set_pt_denom(thisjet->singles.get_pt_denom());

    res4_mainloop<ResultType, JetType, true, ResultType, false>(
            ans,                        //ans
            &unmatched,                 //unmatched_gen
            nullptr,                    //transfer_ans
            nullptr,                    //untransfered_reco
            nullptr,                    //untransfered_gen
            nullptr,                    //thisjet_reco
            thisjet,                    //thisjet_gen
            &matched,                   //matched
            nullptr,                    //adj
            nullptr,                    //axes_reco
            axes,                       //axes_gen
            tolerance2,                 //tolerance2
            tri_tolerance);             //tri_tolerance
}

void EEC::Res4Calculator::compute_JIT_matched(
        const simon::jet& J,
        const std::vector<bool>& matched,
        Res4Result_Vector& result,
        Res4Result_Vector& unmatched) const noexcept {

    call_res4_matched<Res4Result_Vector, EECjet_JIT>(
            result, 
            unmatched,
            J, matched, nt, 
            axes, 
            tolerance2, 
            tri_tolerance);
}

void EEC::Res4Calculator::compute_JIT_matched(
        const simon::jet& J,
        const std::vector<bool>& matched,
        Res4Result_MultiArray& result,
        Res4Result_MultiArray& unmatched) const noexcept {

    call_res4_matched<Res4Result_MultiArray, EECjet_JIT>(
            result, 
            unmatched,
            J, matched, nt, 
            axes, 
            tolerance2, 
            tri_tolerance);
}

void EEC::Res4Calculator::compute_precomputed_matched(
        const simon::jet& J,
        const std::vector<bool>& matched,
        Res4Result_Vector& result,
        Res4Result_Vector& unmatched) const noexcept {

    call_res4_matched<Res4Result_Vector, EECjet_Precomputed>(
            result, 
            unmatched,
            J, matched, nt, 
            axes, 
            tolerance2, 
            tri_tolerance);
}

void EEC::Res4Calculator::compute_precomputed_matched(
        const simon::jet& J,
        const std::vector<bool>& matched,
        Res4Result_MultiArray& result,
        Res4Result_MultiArray& unmatched) const noexcept {

    call_res4_matched<Res4Result_MultiArray, EECjet_Precomputed>(
            result, 
            unmatched,
            J, matched, nt, 
            axes, 
            tolerance2, 
            tri_tolerance);
}
 
template <class TransferResultType, class ResultType, class JetType>
static void call_res4_transfer(
        ResultType& ans,
        ResultType& unmatched_gen,
        TransferResultType& transfer_ans,
        ResultType& untransfered_reco,
        ResultType& untransfered_gen,

        const simon::jet& J_reco,
        const simon::jet& J_gen,
        const EEC::normType& nt,

        const Eigen::MatrixXd& adjmat,

        const EEC::Res4Axes& axes_reco,
        const EEC::Res4Axes& axes_gen,

        const double tolerance2,
        const double tri_tolerance) noexcept {

    auto thisjet_reco = std::make_shared<const JetType>(J_reco, nt);
    auto thisjet_gen = std::make_shared<const JetType>(J_gen, nt);

    ans.set_pt_denom(thisjet_gen->singles.get_pt_denom());
    transfer_ans.set_pt_denom(thisjet_reco->singles.get_pt_denom(),
                              thisjet_gen->singles.get_pt_denom());
    untransfered_reco.set_pt_denom(thisjet_reco->singles.get_pt_denom());
    untransfered_gen.set_pt_denom(thisjet_gen->singles.get_pt_denom());
    unmatched_gen.set_pt_denom(thisjet_gen->singles.get_pt_denom());

    //rescale adjmat
    double denom_gen = thisjet_gen->singles.get_pt_denom();
    double denom_reco = thisjet_reco->singles.get_pt_denom();
    Eigen::VectorXd ptgen = J_gen.ptvec()/denom_gen;
    Eigen::VectorXd ptreco = J_reco.ptvec()/denom_reco;
    Eigen::VectorXd ptfwd = adjmat * ptgen;

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

    res4_mainloop<ResultType, JetType, true, TransferResultType, true>(
            ans,                            //ans
            &unmatched_gen,                 //unmatched_gen
            &transfer_ans,                  //transfer_ans
            &untransfered_reco,             //untransfered_reco
            &untransfered_gen,              //untransfered_gen
            thisjet_reco,                   //thisjet_reco
            thisjet_gen,                    //thisjet_gen
            &matched,                       //matched
            adj,                            //adj
            &axes_reco,                     //axes_reco
            axes_gen,                       //axes_gen
            tolerance2,                     //tolerance2
            tri_tolerance);                 //tri_tolerance
}

void EEC::Res4TransferCalculator::compute_JIT(
        const simon::jet& J_reco,
        const simon::jet& J_gen, 
        const Eigen::MatrixXd& tmat,
        Res4Result_Vector& result,
        EEC::Res4Result_Vector& unmatched_gen,
        Res4TransferResult_Vector& tresult,
        Res4Result_Vector& untransfered_reco,
        Res4Result_Vector& untransfered_gen) const noexcept{

    call_res4_transfer<Res4TransferResult_Vector, Res4Result_Vector, EECjet_JIT>(
            result, unmatched_gen,
            tresult,
            untransfered_reco, untransfered_gen,
            J_reco, J_gen, nt,
            tmat,
            axes_reco, axes_gen,
            tolerance2, tri_tolerance);
}

void EEC::Res4TransferCalculator::compute_JIT(
        const simon::jet& J_reco,
        const simon::jet& J_gen, 
        const Eigen::MatrixXd& tmat,
        Res4Result_MultiArray& result,
        EEC::Res4Result_MultiArray& unmatched_gen,
        Res4TransferResult_Vector& tresult,
        Res4Result_MultiArray& untransfered_reco,
        Res4Result_MultiArray& untransfered_gen) const noexcept{

    call_res4_transfer<Res4TransferResult_Vector, Res4Result_MultiArray, EECjet_JIT>(
            result, unmatched_gen,
            tresult,
            untransfered_reco, untransfered_gen,
            J_reco, J_gen, nt,
            tmat,
            axes_reco, axes_gen,
            tolerance2, tri_tolerance);
}

void EEC::Res4TransferCalculator::compute_JIT(
        const simon::jet& J_reco,
        const simon::jet& J_gen, 
        const Eigen::MatrixXd& tmat,
        Res4Result_Vector& result,
        EEC::Res4Result_Vector& unmatched_gen,
        Res4TransferResult_MultiArray& tresult,
        Res4Result_Vector& untransfered_reco,
        Res4Result_Vector& untransfered_gen) const noexcept{

    call_res4_transfer<Res4TransferResult_MultiArray, Res4Result_Vector, EECjet_JIT>(
            result, unmatched_gen,
            tresult,
            untransfered_reco, untransfered_gen,
            J_reco, J_gen, nt,
            tmat,
            axes_reco, axes_gen,
            tolerance2, tri_tolerance);
}

void EEC::Res4TransferCalculator::compute_JIT(
        const simon::jet& J_reco,
        const simon::jet& J_gen, 
        const Eigen::MatrixXd& tmat,
        Res4Result_MultiArray& result,
        EEC::Res4Result_MultiArray& unmatched_gen,
        Res4TransferResult_MultiArray& tresult,
        Res4Result_MultiArray& untransfered_reco,
        Res4Result_MultiArray& untransfered_gen) const noexcept{

    call_res4_transfer<Res4TransferResult_MultiArray, Res4Result_MultiArray, EECjet_JIT>(
            result, unmatched_gen,
            tresult,
            untransfered_reco, untransfered_gen,
            J_reco, J_gen, nt,
            tmat,
            axes_reco, axes_gen,
            tolerance2, tri_tolerance);
}

void EEC::Res4TransferCalculator::compute_precomputed(
        const simon::jet& J_reco,
        const simon::jet& J_gen, 
        const Eigen::MatrixXd& tmat,
        Res4Result_Vector& result,
        EEC::Res4Result_Vector& unmatched_gen,
        Res4TransferResult_Vector& tresult,
        Res4Result_Vector& untransfered_reco,
        Res4Result_Vector& untransfered_gen) const noexcept{

    call_res4_transfer<Res4TransferResult_Vector, Res4Result_Vector, EECjet_Precomputed>(
            result, unmatched_gen,
            tresult,
            untransfered_reco, untransfered_gen,
            J_reco, J_gen, nt,
            tmat,
            axes_reco, axes_gen,
            tolerance2, tri_tolerance);
}

void EEC::Res4TransferCalculator::compute_precomputed(
        const simon::jet& J_reco,
        const simon::jet& J_gen, 
        const Eigen::MatrixXd& tmat,
        Res4Result_MultiArray& result,
        EEC::Res4Result_MultiArray& unmatched_gen,
        Res4TransferResult_Vector& tresult,
        Res4Result_MultiArray& untransfered_reco,
        Res4Result_MultiArray& untransfered_gen) const noexcept{

    call_res4_transfer<Res4TransferResult_Vector, Res4Result_MultiArray, EECjet_Precomputed>(
            result, unmatched_gen,
            tresult,
            untransfered_reco, untransfered_gen,
            J_reco, J_gen, nt,
            tmat,
            axes_reco, axes_gen,
            tolerance2, tri_tolerance);
}

void EEC::Res4TransferCalculator::compute_precomputed(
        const simon::jet& J_reco,
        const simon::jet& J_gen, 
        const Eigen::MatrixXd& tmat,
        Res4Result_Vector& result,
        EEC::Res4Result_Vector& unmatched_gen,
        Res4TransferResult_MultiArray& tresult,
        Res4Result_Vector& untransfered_reco,
        Res4Result_Vector& untransfered_gen) const noexcept{

    call_res4_transfer<Res4TransferResult_MultiArray, Res4Result_Vector, EECjet_Precomputed>(
            result, unmatched_gen,
            tresult,
            untransfered_reco, untransfered_gen,
            J_reco, J_gen, nt,
            tmat,
            axes_reco, axes_gen,
            tolerance2, tri_tolerance);
}

void EEC::Res4TransferCalculator::compute_precomputed(
        const simon::jet& J_reco,
        const simon::jet& J_gen, 
        const Eigen::MatrixXd& tmat,
        Res4Result_MultiArray& result,
        EEC::Res4Result_MultiArray& unmatched_gen,
        Res4TransferResult_MultiArray& tresult,
        Res4Result_MultiArray& untransfered_reco,
        Res4Result_MultiArray& untransfered_gen) const noexcept{

    call_res4_transfer<Res4TransferResult_MultiArray, Res4Result_MultiArray, EECjet_Precomputed>(
            result, unmatched_gen,
            tresult,
            untransfered_reco, untransfered_gen,
            J_reco, J_gen, nt,
            tmat,
            axes_reco, axes_gen,
            tolerance2, tri_tolerance);
}

#ifdef CMSSW_GIT_HASH

EEC::Res4Calculator::Res4Calculator(const edm::ParameterSet& iConfig) :
    axes(iConfig.getParameter<edm::ParameterSet>("bins")),
    tolerance2(simon::square(iConfig.getParameter<double>("tolerance"))),
    tri_tolerance(iConfig.getParameter<double>("tri_tolerance")),
    nt(normType_from_string(iConfig.getParameter<std::string>("normType"))) {}

void EEC::Res4Calculator::fillPSetDescription(edm::ParameterSetDescription& desc) {

    edm::ParameterSetDescription binsDesc; 
    Res4Axes::fillPSetDescription(binsDesc);
    desc.add<edm::ParameterSetDescription>("bins", binsDesc);

    desc.add<double>("tolerance");
    desc.add<double>("tri_tolerance");
    desc.add<std::string>("normType");
}

EEC::Res4TransferCalculator::Res4TransferCalculator(const edm::ParameterSet& iConfig):
    axes_reco(iConfig.getParameter<edm::ParameterSet>("bins_reco")),
    axes_gen(iConfig.getParameter<edm::ParameterSet>("bins_gen")),
    tolerance2(simon::square(iConfig.getParameter<double>("tolerance"))),
    tri_tolerance(iConfig.getParameter<double>("tri_tolerance")),
    nt(normType_from_string(iConfig.getParameter<std::string>("normType"))) {}

void EEC::Res4TransferCalculator::fillPSetDescription(edm::ParameterSetDescription& desc) {

    edm::ParameterSetDescription binsRecoDesc;
    Res4Axes::fillPSetDescription(binsRecoDesc);
    desc.add<edm::ParameterSetDescription>("bins_reco", binsRecoDesc);

    edm::ParameterSetDescription binsGenDesc;
    Res4Axes::fillPSetDescription(binsGenDesc);
    desc.add<edm::ParameterSetDescription>("bins_gen", binsGenDesc);

    desc.add<double>("tolerance");
    desc.add<double>("tri_tolerance");
    desc.add<std::string>("normType");
}

#endif
