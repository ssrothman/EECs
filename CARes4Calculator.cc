#include "CARes4Calculator.h"
#include "Adjacency.h"
#include "CARes4Result.h"
//#include "CARes4TransferResult.h"
#include "CARes4Backend.h"

#include "SRothman/SimonTools/src/deltaR.h"
#include "SRothman/SimonTools/src/histutil.h"
#include "SRothman/SimonTools/src/util.h"

#include <memory>

#include <array>

template <class ResultType, class JetType>
static void call_CAres4(
        ResultType& ans,

        const simon::jet& J,
        const EEC::normType& nt,

        const EEC::CARes4Axes& axes) noexcept {

    const JetType thisjet(J, nt);
    ans.set_pt_denom(thisjet.singles.get_pt_denom());
    
    CAres4_mainloop<ResultType, JetType, false, ResultType, false>(
            ans,                                //ans
            nullptr,                            //unmatched_gen
            nullptr,                            //transfer_ans
            nullptr,                            //thisjet_reco
            thisjet,                            //thisjet_gen
            nullptr,                            //matched
            nullptr,                            //adj
            nullptr,                            //axes_reco
            axes);                              //axes_gen
}

void EEC::CARes4Calculator::compute_JIT(
        const simon::jet& J,
        CARes4Result_Vector& result) const noexcept {
    
    call_CAres4<CARes4Result_Vector, EECjet_JIT>(
            result, 
            J, nt, 
            axes); 
}

void EEC::CARes4Calculator::compute_JIT(
        const simon::jet& J,
        CARes4Result_MultiArray& result) const noexcept {

    call_CAres4<CARes4Result_MultiArray, EECjet_JIT>(
            result, 
            J, nt, 
            axes); 
}

void EEC::CARes4Calculator::compute_precomputed(
        const simon::jet& J,
        CARes4Result_Vector& result) const noexcept {
    
    call_CAres4<CARes4Result_Vector, EECjet_Precomputed>(
            result, 
            J, nt, 
            axes); 
}

void EEC::CARes4Calculator::compute_precomputed(
        const simon::jet& J,
        CARes4Result_MultiArray& result) const noexcept {

    call_CAres4<CARes4Result_MultiArray, EECjet_Precomputed>(
            result, 
            J, nt, 
            axes); 
}

template <class ResultType, class JetType>
static void call_CAres4_matched(
        ResultType& ans,
        ResultType& unmatched,

        const simon::jet& J,
        const std::vector<bool>& matched,
        const EEC::normType& nt,

        const EEC::CARes4Axes& axes) noexcept {

    const JetType thisjet(J, nt);
    ans.set_pt_denom(thisjet.singles.get_pt_denom());
    unmatched.set_pt_denom(thisjet.singles.get_pt_denom());

    CAres4_mainloop<ResultType, JetType, true, ResultType, false>(
            ans,                        //ans
            &unmatched,                 //unmatched_gen
            nullptr,                    //transfer_ans
            nullptr,                    //thisjet_reco
            thisjet,                    //thisjet_gen
            &matched,                   //matched
            nullptr,                    //adj
            nullptr,                    //axes_reco
            axes);                       //axes_gen
}

void EEC::CARes4Calculator::compute_JIT_matched(
        const simon::jet& J,
        const std::vector<bool>& matched,
        CARes4Result_Vector& result,
        CARes4Result_Vector& unmatched) const noexcept {

    call_CAres4_matched<CARes4Result_Vector, EECjet_JIT>(
            result, 
            unmatched,
            J, matched, nt, 
            axes); 
}

void EEC::CARes4Calculator::compute_JIT_matched(
        const simon::jet& J,
        const std::vector<bool>& matched,
        CARes4Result_MultiArray& result,
        CARes4Result_MultiArray& unmatched) const noexcept {

    call_CAres4_matched<CARes4Result_MultiArray, EECjet_JIT>(
            result, 
            unmatched,
            J, matched, nt, 
            axes); 
}

void EEC::CARes4Calculator::compute_precomputed_matched(
        const simon::jet& J,
        const std::vector<bool>& matched,
        CARes4Result_Vector& result,
        CARes4Result_Vector& unmatched) const noexcept {

    call_CAres4_matched<CARes4Result_Vector, EECjet_Precomputed>(
            result, 
            unmatched,
            J, matched, nt, 
            axes); 
}

void EEC::CARes4Calculator::compute_precomputed_matched(
        const simon::jet& J,
        const std::vector<bool>& matched,
        CARes4Result_MultiArray& result,
        CARes4Result_MultiArray& unmatched) const noexcept {

    call_CAres4_matched<CARes4Result_MultiArray, EECjet_Precomputed>(
            result, 
            unmatched,
            J, matched, nt, 
            axes); 
}
 
template <class TransferResultType, class ResultType, class JetType>
static void call_CAres4_transfer(
        ResultType& ans,
        ResultType& unmatched_gen,
        TransferResultType& transfer_ans,

        const simon::jet& J_reco,
        const simon::jet& J_gen,
        const EEC::normType& nt,

        const Eigen::MatrixXd& adjmat,

        const EEC::CARes4Axes& axes_reco,
        const EEC::CARes4Axes& axes_gen) noexcept {

    const JetType thisjet_reco(J_reco, nt);
    const JetType thisjet_gen(J_gen, nt);

    ans.set_pt_denom(thisjet_gen.singles.get_pt_denom());
    transfer_ans.set_pt_denom(thisjet_reco.singles.get_pt_denom(),
                              thisjet_gen.singles.get_pt_denom());
    unmatched_gen.set_pt_denom(thisjet_gen.singles.get_pt_denom());

    //rescale adjmat
    double denom_gen = thisjet_gen.singles.get_pt_denom();
    double denom_reco = thisjet_reco.singles.get_pt_denom();
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

    const EEC::Adjacency adj(adjmat_copy);

    std::vector<bool> matched(J_gen.nPart, false);
    for(unsigned iGenPart=0; iGenPart<J_gen.nPart; ++iGenPart){
        for(unsigned iRecoPart=0; iRecoPart<J_reco.nPart; ++iRecoPart){
            if(adjmat_copy(iRecoPart, iGenPart) > 0){
                matched[iGenPart] = true;
                break;
            }
        }
    }

    CAres4_mainloop<ResultType, JetType, true, TransferResultType, true>(
            ans,                            //ans
            &unmatched_gen,                 //unmatched_gen
            &transfer_ans,                  //transfer_ans
            &thisjet_reco,                   //thisjet_reco
            thisjet_gen,                    //thisjet_gen
            &matched,                       //matched
            &adj,                            //adj
            &axes_reco,                     //axes_reco
            axes_gen);                       //axes_gen
}

void EEC::CARes4TransferCalculator::compute_JIT(
        const simon::jet& J_reco,
        const simon::jet& J_gen, 
        const Eigen::MatrixXd& tmat,
        CARes4Result_Vector& result,
        EEC::CARes4Result_Vector& unmatched_gen,
        CARes4TransferResult_Vector& tresult) const noexcept {

    call_CAres4_transfer<CARes4TransferResult_Vector, CARes4Result_Vector, EECjet_JIT>(
            result, unmatched_gen,
            tresult,
            J_reco, J_gen, nt,
            tmat,
            axes_reco, axes_gen);
}

void EEC::CARes4TransferCalculator::compute_JIT(
        const simon::jet& J_reco,
        const simon::jet& J_gen, 
        const Eigen::MatrixXd& tmat,
        CARes4Result_MultiArray& result,
        EEC::CARes4Result_MultiArray& unmatched_gen,
        CARes4TransferResult_Vector& tresult) const noexcept {

    call_CAres4_transfer<CARes4TransferResult_Vector, CARes4Result_MultiArray, EECjet_JIT>(
            result, unmatched_gen,
            tresult,
            J_reco, J_gen, nt,
            tmat,
            axes_reco, axes_gen);
}

void EEC::CARes4TransferCalculator::compute_JIT(
        const simon::jet& J_reco,
        const simon::jet& J_gen, 
        const Eigen::MatrixXd& tmat,
        CARes4Result_Vector& result,
        EEC::CARes4Result_Vector& unmatched_gen,
        CARes4TransferResult_MultiArray& tresult) const noexcept {

    call_CAres4_transfer<CARes4TransferResult_MultiArray, CARes4Result_Vector, EECjet_JIT>(
            result, unmatched_gen,
            tresult,
            J_reco, J_gen, nt,
            tmat,
            axes_reco, axes_gen);
}

void EEC::CARes4TransferCalculator::compute_JIT(
        const simon::jet& J_reco,
        const simon::jet& J_gen, 
        const Eigen::MatrixXd& tmat,
        CARes4Result_MultiArray& result,
        EEC::CARes4Result_MultiArray& unmatched_gen,
        CARes4TransferResult_MultiArray& tresult) const noexcept {

    call_CAres4_transfer<CARes4TransferResult_MultiArray, CARes4Result_MultiArray, EECjet_JIT>(
            result, unmatched_gen,
            tresult,
            J_reco, J_gen, nt,
            tmat,
            axes_reco, axes_gen);
}

void EEC::CARes4TransferCalculator::compute_precomputed(
        const simon::jet& J_reco,
        const simon::jet& J_gen, 
        const Eigen::MatrixXd& tmat,
        CARes4Result_Vector& result,
        EEC::CARes4Result_Vector& unmatched_gen,
        CARes4TransferResult_Vector& tresult) const noexcept {

    call_CAres4_transfer<CARes4TransferResult_Vector, CARes4Result_Vector, EECjet_Precomputed>(
            result, unmatched_gen,
            tresult,
            J_reco, J_gen, nt,
            tmat,
            axes_reco, axes_gen);
}

void EEC::CARes4TransferCalculator::compute_precomputed(
        const simon::jet& J_reco,
        const simon::jet& J_gen, 
        const Eigen::MatrixXd& tmat,
        CARes4Result_MultiArray& result,
        EEC::CARes4Result_MultiArray& unmatched_gen,
        CARes4TransferResult_Vector& tresult) const noexcept {

    call_CAres4_transfer<CARes4TransferResult_Vector, CARes4Result_MultiArray, EECjet_Precomputed>(
            result, unmatched_gen,
            tresult,
            J_reco, J_gen, nt,
            tmat,
            axes_reco, axes_gen);
}

void EEC::CARes4TransferCalculator::compute_precomputed(
        const simon::jet& J_reco,
        const simon::jet& J_gen, 
        const Eigen::MatrixXd& tmat,
        CARes4Result_Vector& result,
        EEC::CARes4Result_Vector& unmatched_gen,
        CARes4TransferResult_MultiArray& tresult) const noexcept {

    call_CAres4_transfer<CARes4TransferResult_MultiArray, CARes4Result_Vector, EECjet_Precomputed>(
            result, unmatched_gen,
            tresult,
            J_reco, J_gen, nt,
            tmat,
            axes_reco, axes_gen);
}

void EEC::CARes4TransferCalculator::compute_precomputed(
        const simon::jet& J_reco,
        const simon::jet& J_gen, 
        const Eigen::MatrixXd& tmat,
        CARes4Result_MultiArray& result,
        EEC::CARes4Result_MultiArray& unmatched_gen,
        CARes4TransferResult_MultiArray& tresult) const noexcept {

    call_CAres4_transfer<CARes4TransferResult_MultiArray, CARes4Result_MultiArray, EECjet_Precomputed>(
            result, unmatched_gen,
            tresult,
            J_reco, J_gen, nt,
            tmat,
            axes_reco, axes_gen);
}

#ifdef CMSSW_GIT_HASH

EEC::CARes4Calculator::CARes4Calculator(const edm::ParameterSet& iConfig) :
    axes(iConfig.getParameter<edm::ParameterSet>("bins")),
    nt(normType_from_string(iConfig.getParameter<std::string>("normType"))) {}

void EEC::CARes4Calculator::fillPSetDescription(edm::ParameterSetDescription& desc) {

    edm::ParameterSetDescription binsDesc; 
    CARes4Axes::fillPSetDescription(binsDesc);
    desc.add<edm::ParameterSetDescription>("bins", binsDesc);

    desc.add<std::string>("normType");
}

EEC::CARes4TransferCalculator::CARes4TransferCalculator(const edm::ParameterSet& iConfig):
    axes_reco(iConfig.getParameter<edm::ParameterSet>("bins_reco")),
    axes_gen(iConfig.getParameter<edm::ParameterSet>("bins_gen")),
    nt(normType_from_string(iConfig.getParameter<std::string>("normType"))) {}

void EEC::CARes4TransferCalculator::fillPSetDescription(edm::ParameterSetDescription& desc) {

    edm::ParameterSetDescription binsRecoDesc;
    CARes4Axes::fillPSetDescription(binsRecoDesc);
    desc.add<edm::ParameterSetDescription>("bins_reco", binsRecoDesc);

    edm::ParameterSetDescription binsGenDesc;
    CARes4Axes::fillPSetDescription(binsGenDesc);
    desc.add<edm::ParameterSetDescription>("bins_gen", binsGenDesc);

    desc.add<std::string>("normType");
}

#endif
