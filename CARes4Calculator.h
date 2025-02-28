#ifndef SROTHMAN_EECS_V2_RES4CALCULATOR_H
#define SROTHMAN_EECS_V2_RES4CALCULATOR_H

#include "CARes4Axes.h"
#include "EECjet.h"

#include "CARes4Result.h"
#include "CARes4TransferResult.h"

#ifdef CMSSW_GIT_HASH
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#endif

namespace EEC{
    class CARes4Calculator {
    public:
        CARes4Calculator(
                const std::vector<double>& R,
                const std::vector<double>& r_chain,
                const std::vector<double>& c_chain,
                const std::vector<double>& r_symmetric,
                const std::vector<double>& c_symmetric,
                const normType nt) :
            axes(R, 
                 r_chain, c_chain,
                 r_symmetric, c_symmetric),
            nt(nt) {}

#ifdef CMSSW_GIT_HASH
        CARes4Calculator(const edm::ParameterSet& iConfig);
        static void fillPSetDescription(edm::ParameterSetDescription& desc);
#endif

        void compute_JIT(
                const simon::jet& J, 
                CARes4Result_Vector& result) const noexcept;

        void compute_JIT(
                const simon::jet& J, 
                CARes4Result_MultiArray& result) const noexcept;

        void compute_precomputed(
                const simon::jet& J, 
                CARes4Result_Vector& result) const noexcept;

        void compute_precomputed(
                const simon::jet& J, 
                CARes4Result_MultiArray& result) const noexcept;

        void compute_JIT_matched(
                const simon::jet& J,
                const std::vector<bool>& matched,
                CARes4Result_Vector& result,
                CARes4Result_Vector& unmatched) const noexcept;

        void compute_JIT_matched(
                const simon::jet& J,
                const std::vector<bool>& matched,
                CARes4Result_MultiArray& result,
                CARes4Result_MultiArray& unmatched) const noexcept;

        void compute_precomputed_matched(
                const simon::jet& J,
                const std::vector<bool>& matched,
                CARes4Result_Vector& result,
                CARes4Result_Vector& unmatched) const noexcept;

        void compute_precomputed_matched(
                const simon::jet& J,
                const std::vector<bool>& matched,
                CARes4Result_MultiArray& result,
                CARes4Result_MultiArray& unmatched) const noexcept;

        const CARes4Axes& get_axes() const noexcept{
            return axes;
        }
    private:
        const CARes4Axes axes;

        const EEC::normType nt;
    };

    class CARes4TransferCalculator {
    public:
        CARes4TransferCalculator(
                const std::vector<double>& R_reco,
                const std::vector<double>& r_chain_reco,
                const std::vector<double>& c_chain_reco,
                const std::vector<double>& r_symmetric_reco,
                const std::vector<double>& c_symmetric_reco,
                const std::vector<double>& R_gen,
                const std::vector<double>& r_chain_gen,
                const std::vector<double>& c_chain_gen,
                const std::vector<double>& r_symmetric_gen,
                const std::vector<double>& c_symmetric_gen,
                const normType nt) :
            axes_reco(R_reco, 
                      r_chain_reco, c_chain_reco, 
                      r_symmetric_reco, c_symmetric_reco),
            axes_gen(R_gen, 
                     r_chain_gen, c_chain_gen, 
                     r_symmetric_gen, c_symmetric_gen),
            nt(nt) {}

#ifdef CMSSW_GIT_HASH
        CARes4TransferCalculator(const edm::ParameterSet& iConfig);
        static void fillPSetDescription(edm::ParameterSetDescription& desc);
#endif

        void compute_JIT(
                const simon::jet& J_reco,
                const simon::jet& J_gen, 
                const Eigen::MatrixXd& tmat,
                CARes4Result_Vector& result,
                EEC::CARes4Result_Vector& unmatched_gen,
                CARes4TransferResult_Vector& tresult) const noexcept;

        void compute_JIT(
                const simon::jet& J_reco,
                const simon::jet& J_gen, 
                const Eigen::MatrixXd& tmat,
                CARes4Result_MultiArray& result,
                EEC::CARes4Result_MultiArray& unmatched_gen,
                CARes4TransferResult_Vector& tresult) const noexcept;

        void compute_JIT(
                const simon::jet& J_reco,
                const simon::jet& J_gen, 
                const Eigen::MatrixXd& tmat,
                CARes4Result_Vector& result,
                EEC::CARes4Result_Vector& unmatched_gen,
                CARes4TransferResult_MultiArray& tresult) const noexcept;

        void compute_JIT(
                const simon::jet& J_reco,
                const simon::jet& J_gen, 
                const Eigen::MatrixXd& tmat,
                CARes4Result_MultiArray& result,
                EEC::CARes4Result_MultiArray& unmatched_gen,
                CARes4TransferResult_MultiArray& tresult) const noexcept;

        void compute_precomputed(
                const simon::jet& J_reco,
                const simon::jet& J_gen, 
                const Eigen::MatrixXd& tmat,
                CARes4Result_Vector& result,
                EEC::CARes4Result_Vector& unmatched_gen,
                CARes4TransferResult_Vector& tresult) const noexcept;

        void compute_precomputed(
                const simon::jet& J_reco,
                const simon::jet& J_gen, 
                const Eigen::MatrixXd& tmat,
                CARes4Result_MultiArray& result,
                EEC::CARes4Result_MultiArray& unmatched_gen,
                CARes4TransferResult_Vector& tresult) const noexcept;

        void compute_precomputed(
                const simon::jet& J_reco,
                const simon::jet& J_gen, 
                const Eigen::MatrixXd& tmat,
                CARes4Result_Vector& result,
                EEC::CARes4Result_Vector& unmatched_gen,
                CARes4TransferResult_MultiArray& tresult) const noexcept;

        void compute_precomputed(
                const simon::jet& J_reco,
                const simon::jet& J_gen, 
                const Eigen::MatrixXd& tmat,
                CARes4Result_MultiArray& result,
                EEC::CARes4Result_MultiArray& unmatched_gen,
                CARes4TransferResult_MultiArray& tresult) const noexcept;

        const CARes4Axes& get_axes_reco() const noexcept{
            return axes_reco;
        }

        const CARes4Axes& get_axes_gen() const noexcept{
            return axes_gen;
        }

    private:
        const CARes4Axes axes_reco, axes_gen;

        const EEC::normType nt;
    };
};

#endif
