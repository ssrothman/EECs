#ifndef SROTHMAN_EECS_V2_RES4_CALCULATOR_H
#define SROTHMAN_EECS_V2_RES4_CALCULATOR_H

#include "Res4Axes.h"
#include "EECjet.h"

#include "Res4Result.h"
#include "Res4TransferResult.h"

#ifdef CMSSW_GIT_HASH
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#endif

namespace EEC{
    //this is in principle configurable
    //but it is important that the target ratios not be 1
    static constexpr double TRIANGLE_RATIO_LM = 5.0/4.0;
    static constexpr double TRIANGLE_RATIO_MS = 4.0/3.0;

    class Res4Calculator {
    public:
        Res4Calculator(
                const std::vector<double>& R,
                const std::vector<double>& r_dipole,
                const std::vector<double>& c_dipole,
                const std::vector<double>& r_tee,
                const std::vector<double>& c_tee,
                const std::vector<double>& r_triangle,
                const std::vector<double>& c_triangle,
                const double tolerance,
                const double tri_tolerance,
                const normType nt) :
            axes(R, 
                 r_dipole, c_dipole, 
                 r_tee, c_tee, 
                 r_triangle, c_triangle),
            tolerance2(tolerance*tolerance),
            tri_tolerance(tri_tolerance),
            nt(nt) {}

#ifdef CMSSW_GIT_HASH
        Res4Calculator(const edm::ParameterSet& iConfig);
        static void fillPSetDescription(edm::ParameterSetDescription& desc);
#endif

        void compute_JIT(
                const simon::jet& J, 
                Res4Result_Vector& result) const noexcept;

        void compute_JIT(
                const simon::jet& J, 
                Res4Result_Unbinned& result) const noexcept;

        void compute_JIT(
                const simon::jet& J, 
                Res4Result_MultiArray& result) const noexcept;

        void compute_precomputed(
                const simon::jet& J, 
                Res4Result_Vector& result) const noexcept;

        void compute_precomputed(
                const simon::jet& J, 
                Res4Result_Unbinned& result) const noexcept;

        void compute_precomputed(
                const simon::jet& J, 
                Res4Result_MultiArray& result) const noexcept;

        void compute_JIT_matched(
                const simon::jet& J,
                const std::vector<bool>& matched,
                Res4Result_Vector& result,
                Res4Result_Vector& unmatched) const noexcept;

        void compute_JIT_matched(
                const simon::jet& J,
                const std::vector<bool>& matched,
                Res4Result_Unbinned& result,
                Res4Result_Unbinned& unmatched) const noexcept;

        void compute_JIT_matched(
                const simon::jet& J,
                const std::vector<bool>& matched,
                Res4Result_MultiArray& result,
                Res4Result_MultiArray& unmatched) const noexcept;

        void compute_precomputed_matched(
                const simon::jet& J,
                const std::vector<bool>& matched,
                Res4Result_Vector& result,
                Res4Result_Vector& unmatched) const noexcept;

        void compute_precomputed_matched(
                const simon::jet& J,
                const std::vector<bool>& matched,
                Res4Result_Unbinned& result,
                Res4Result_Unbinned& unmatched) const noexcept;

        void compute_precomputed_matched(
                const simon::jet& J,
                const std::vector<bool>& matched,
                Res4Result_MultiArray& result,
                Res4Result_MultiArray& unmatched) const noexcept;

        const Res4Axes& get_axes() const noexcept{
            return axes;
        }
    private:
        const Res4Axes axes;
        const double tolerance2;
        const double tri_tolerance;

        const EEC::normType nt;
    };

    class Res4TransferCalculator {
    public:
        constexpr static bool HAS_UNTRANSFERED = true;

        Res4TransferCalculator(
                const std::vector<double>& R_reco,
                const std::vector<double>& r_dipole_reco,
                const std::vector<double>& c_dipole_reco,
                const std::vector<double>& r_tee_reco,
                const std::vector<double>& c_tee_reco,
                const std::vector<double>& r_triangle_reco,
                const std::vector<double>& c_triangle_reco,
                const std::vector<double>& R_gen,
                const std::vector<double>& r_dipole_gen,
                const std::vector<double>& c_dipole_gen,
                const std::vector<double>& r_tee_gen,
                const std::vector<double>& c_tee_gen,
                const std::vector<double>& r_triangle_gen,
                const std::vector<double>& c_triangle_gen,
                const double tolerance,
                const double tri_tolerance,
                const normType nt) :
            axes_reco(R_reco, 
                      r_dipole_reco, c_dipole_reco, 
                      r_tee_reco, c_tee_reco, 
                      r_triangle_reco, c_triangle_reco),
            axes_gen(R_gen, 
                     r_dipole_gen, c_dipole_gen, 
                     r_tee_gen, c_tee_gen, 
                     r_triangle_gen, c_triangle_gen),
            tolerance2(tolerance*tolerance),
            tri_tolerance(tri_tolerance),
            nt(nt) {}

#ifdef CMSSW_GIT_HASH
        Res4TransferCalculator(const edm::ParameterSet& iConfig);
        static void fillPSetDescription(edm::ParameterSetDescription& desc);
#endif

        void compute_JIT(
                const simon::jet& J_reco,
                const simon::jet& J_gen, 
                const Eigen::MatrixXd& tmat,
                Res4Result_Vector& result,
                EEC::Res4Result_Vector& unmatched_gen,
                Res4TransferResult_Vector& tresult,
                Res4Result_Vector& untransfered_reco,
                Res4Result_Vector& untransfered_gen) const noexcept;

        void compute_JIT(
                const simon::jet& J_reco,
                const simon::jet& J_gen, 
                const Eigen::MatrixXd& tmat,
                Res4Result_Unbinned& result,
                EEC::Res4Result_Unbinned& unmatched_gen,
                Res4TransferResult_Unbinned& tresult,
                Res4Result_Unbinned& untransfered_reco,
                Res4Result_Unbinned& untransfered_gen) const noexcept;

        void compute_JIT(
                const simon::jet& J_reco,
                const simon::jet& J_gen, 
                const Eigen::MatrixXd& tmat,
                Res4Result_MultiArray& result,
                EEC::Res4Result_MultiArray& unmatched_gen,
                Res4TransferResult_Vector& tresult,
                Res4Result_MultiArray& untransfered_reco,
                Res4Result_MultiArray& untransfered_gen) const noexcept;

        void compute_JIT(
                const simon::jet& J_reco,
                const simon::jet& J_gen, 
                const Eigen::MatrixXd& tmat,
                Res4Result_Vector& result,
                EEC::Res4Result_Vector& unmatched_gen,
                Res4TransferResult_MultiArray& tresult,
                Res4Result_Vector& untransfered_reco,
                Res4Result_Vector& untransfered_gen) const noexcept;

        void compute_JIT(
                const simon::jet& J_reco,
                const simon::jet& J_gen, 
                const Eigen::MatrixXd& tmat,
                Res4Result_MultiArray& result,
                EEC::Res4Result_MultiArray& unmatched_gen,
                Res4TransferResult_MultiArray& tresult,
                Res4Result_MultiArray& untransfered_reco,
                Res4Result_MultiArray& untransfered_gen) const noexcept;

        void compute_precomputed(
                const simon::jet& J_reco,
                const simon::jet& J_gen, 
                const Eigen::MatrixXd& tmat,
                Res4Result_Vector& result,
                EEC::Res4Result_Vector& unmatched_gen,
                Res4TransferResult_Vector& tresult,
                Res4Result_Vector& untransfered_reco,
                Res4Result_Vector& untransfered_gen) const noexcept;

        void compute_precomputed(
                const simon::jet& J_reco,
                const simon::jet& J_gen, 
                const Eigen::MatrixXd& tmat,
                Res4Result_Unbinned& result,
                EEC::Res4Result_Unbinned& unmatched_gen,
                Res4TransferResult_Unbinned& tresult,
                Res4Result_Unbinned& untransfered_reco,
                Res4Result_Unbinned& untransfered_gen) const noexcept;

        void compute_precomputed(
                const simon::jet& J_reco,
                const simon::jet& J_gen, 
                const Eigen::MatrixXd& tmat,
                Res4Result_MultiArray& result,
                EEC::Res4Result_MultiArray& unmatched_gen,
                Res4TransferResult_Vector& tresult,
                Res4Result_MultiArray& untransfered_reco,
                Res4Result_MultiArray& untransfered_gen) const noexcept;

        void compute_precomputed(
                const simon::jet& J_reco,
                const simon::jet& J_gen, 
                const Eigen::MatrixXd& tmat,
                Res4Result_Vector& result,
                EEC::Res4Result_Vector& unmatched_gen,
                Res4TransferResult_MultiArray& tresult,
                Res4Result_Vector& untransfered_reco,
                Res4Result_Vector& untransfered_gen) const noexcept;

        void compute_precomputed(
                const simon::jet& J_reco,
                const simon::jet& J_gen, 
                const Eigen::MatrixXd& tmat,
                Res4Result_MultiArray& result,
                EEC::Res4Result_MultiArray& unmatched_gen,
                Res4TransferResult_MultiArray& tresult,
                Res4Result_MultiArray& untransfered_reco,
                Res4Result_MultiArray& untransfered_gen) const noexcept;

        const Res4Axes& get_axes_reco() const noexcept{
            return axes_reco;
        }

        const Res4Axes& get_axes_gen() const noexcept{
            return axes_gen;
        }

    private:
        const Res4Axes axes_reco, axes_gen;
        const double tolerance2;
        const double tri_tolerance;

        const EEC::normType nt;
    };
};

#endif
