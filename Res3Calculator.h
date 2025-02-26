#ifndef SROTHMAN_ECS_V2_RES3CALCULATOR_H
#define SROTHMAN_ECS_V2_RES3CALCULATOR_H

#include "Res3Axes.h"
#include "EECjet.h"

#include "Res3Result_forward.h"
#include "Res3TransferResult_forward.h"

#ifdef CMSSW_GIT_HASH
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#endif

namespace EEC{
    class Res3Calculator {
    public:
        Res3Calculator(
                const std::vector<double>& R,
                const std::vector<double>& r,
                const std::vector<double>& c,
                const normType nt) :
            axes(R, r, c),
            nt(nt) {}

#ifdef CMSSW_GIT_HASH
        Res3Calculator(const edm::ParameterSet& iConfig);
        static void fillPSetDescription(edm::ParameterSetDescription& desc);
#endif

        void compute_JIT(
                const simon::jet& J, 
                Res3Result_Vector& result) const noexcept;

        void compute_JIT(
                const simon::jet& J, 
                Res3Result_MultiArray& result) const noexcept;

        void compute_precomputed(
                const simon::jet& J, 
                Res3Result_Vector& result) const noexcept;

        void compute_precomputed(
                const simon::jet& J, 
                Res3Result_MultiArray& result) const noexcept;

        void compute_JIT_matched(
                const simon::jet& J,
                const std::vector<bool>& matched,
                Res3Result_Vector& result,
                Res3Result_Vector& unmatched) const noexcept;

        void compute_JIT_matched(
                const simon::jet& J,
                const std::vector<bool>& matched,
                Res3Result_MultiArray& result,
                Res3Result_MultiArray& unmatched) const noexcept;

        void compute_precomputed_matched(
                const simon::jet& J,
                const std::vector<bool>& matched,
                Res3Result_Vector& result,
                Res3Result_Vector& unmatched) const noexcept;

        void compute_precomputed_matched(
                const simon::jet& J,
                const std::vector<bool>& matched,
                Res3Result_MultiArray& result,
                Res3Result_MultiArray& unmatched) const noexcept;

        const Res3Axes& get_axes() const noexcept{
            return axes;
        }
    private:
        const Res3Axes axes;
        const EEC::normType nt;
    };

    class Res3TransferCalculator {
    public:
        Res3TransferCalculator(
                const std::vector<double>& R_reco,
                const std::vector<double>& r_reco,
                const std::vector<double>& c_reco,
                const std::vector<double>& R_gen,
                const std::vector<double>& r_gen,
                const std::vector<double>& c_gen,
                const normType nt) :
            axes_reco(R_reco, r_reco, c_reco),
            axes_gen(R_gen, r_gen, c_gen),
            nt(nt) {}

#ifdef CMSSW_GIT_HASH
        Res3TransferCalculator(const edm::ParameterSet& iConfig);
        static void fillPSetDescription(edm::ParameterSetDescription& desc);
#endif

        void compute_JIT(
                const simon::jet& J_reco,
                const simon::jet& J_gen, 
                const Eigen::MatrixXd& tmat,
                Res3Result_Vector& result,
                EEC::Res3Result_Vector& unmatched_gen,
                Res3TransferResult_Vector& tresult) const noexcept;

        void compute_JIT(
                const simon::jet& J_reco,
                const simon::jet& J_gen, 
                const Eigen::MatrixXd& tmat,
                Res3Result_MultiArray& result,
                EEC::Res3Result_MultiArray& unmatched_gen,
                Res3TransferResult_Vector& tresult) const noexcept;

        void compute_JIT(
                const simon::jet& J_reco,
                const simon::jet& J_gen, 
                const Eigen::MatrixXd& tmat,
                Res3Result_Vector& result,
                EEC::Res3Result_Vector& unmatched_gen,
                Res3TransferResult_MultiArray& tresult) const noexcept;

        void compute_JIT(
                const simon::jet& J_reco,
                const simon::jet& J_gen, 
                const Eigen::MatrixXd& tmat,
                Res3Result_MultiArray& result,
                EEC::Res3Result_MultiArray& unmatched_gen,
                Res3TransferResult_MultiArray& tresult) const noexcept;

        void compute_precomputed(
                const simon::jet& J_reco,
                const simon::jet& J_gen, 
                const Eigen::MatrixXd& tmat,
                Res3Result_Vector& result,
                EEC::Res3Result_Vector& unmatched_gen,
                Res3TransferResult_Vector& tresult) const noexcept;

        void compute_precomputed(
                const simon::jet& J_reco,
                const simon::jet& J_gen, 
                const Eigen::MatrixXd& tmat,
                Res3Result_MultiArray& result,
                EEC::Res3Result_MultiArray& unmatched_gen,
                Res3TransferResult_Vector& tresult) const noexcept;

        void compute_precomputed(
                const simon::jet& J_reco,
                const simon::jet& J_gen, 
                const Eigen::MatrixXd& tmat,
                Res3Result_Vector& result,
                EEC::Res3Result_Vector& unmatched_gen,
                Res3TransferResult_MultiArray& tresult) const noexcept;

        void compute_precomputed(
                const simon::jet& J_reco,
                const simon::jet& J_gen, 
                const Eigen::MatrixXd& tmat,
                Res3Result_MultiArray& result,
                EEC::Res3Result_MultiArray& unmatched_gen,
                Res3TransferResult_MultiArray& tresult) const noexcept;

        const Res3Axes& get_axes_reco() const noexcept{
            return axes_reco;
        }

        const Res3Axes& get_axes_gen() const noexcept{
            return axes_gen;
        }

    private:
        const Res3Axes axes_reco, axes_gen;
        const EEC::normType nt;
    };
};

#endif

