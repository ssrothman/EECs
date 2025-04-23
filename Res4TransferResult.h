#ifndef SROTHMAN_EECS_RES4_TRANSFER_RESULT_H
#define SROTHMAN_EECS_RES4_TRANSFER_RESULT_H

#include "usings.h"
#include "Res4Axes.h"

#include "SRothman/SimonTools/src/histutil.h"
#include "ResTransferResultContainers.h"


namespace EEC{
    template <class TransferContainer>
    class Res4TransferResult{
    public:
        using T = typename TransferContainer::TYPE;
        constexpr static bool SHOULD_BIN = TransferContainer::SHOULD_BIN;

        Res4TransferResult() noexcept  :
            dipole(), tee(), triangle(),
            pt_denom_set(false),
            pt_denom_reco(-1), pt_denom_gen(-1) {}

        Res4TransferResult(
                const size_t nR_reco,
                const size_t nr_dipole_reco,
                const size_t nc_dipole_reco,
                const size_t nr_tee_reco,
                const size_t nc_tee_reco,
                const size_t nr_triangle_reco,
                const size_t nc_triangle_reco,
                const size_t nR_gen,
                const size_t nr_dipole_gen,
                const size_t nc_dipole_gen,
                const size_t nr_tee_gen,
                const size_t nc_tee_gen,
                const size_t nr_triangle_gen,
                const size_t nc_triangle_gen) noexcept  :
            dipole(nR_reco, nr_dipole_reco, nc_dipole_reco, 
                    nR_gen, nr_dipole_gen, nc_dipole_gen),
            tee(nR_reco, nr_tee_reco, nc_tee_reco,
                nR_gen, nr_tee_gen, nc_tee_gen),
            triangle(nR_reco, nr_triangle_reco, nc_triangle_reco,
                     nR_gen, nr_triangle_gen, nc_triangle_gen),
            pt_denom_set(false),
            pt_denom_reco(-1), pt_denom_gen(-1) {}

        Res4TransferResult(
                const Res4Axes& axes_reco,
                const Res4Axes& axes_gen) noexcept  :
            Res4TransferResult(
                    simon::AXextent(axes_reco.R),
                    simon::AXextent(axes_reco.r_dipole),
                    simon::AXextent(axes_reco.c_dipole),
                    simon::AXextent(axes_reco.r_tee),
                    simon::AXextent(axes_reco.c_tee),
                    simon::AXextent(axes_reco.r_triangle),
                    simon::AXextent(axes_reco.c_triangle),
                    simon::AXextent(axes_gen.R),
                    simon::AXextent(axes_gen.r_dipole),
                    simon::AXextent(axes_gen.c_dipole),
                    simon::AXextent(axes_gen.r_tee),
                    simon::AXextent(axes_gen.c_tee),
                    simon::AXextent(axes_gen.r_triangle),
                    simon::AXextent(axes_gen.c_triangle)) {}

        template <typename CALCULATOR>
        Res4TransferResult(const CALCULATOR& calculator) noexcept:
            Res4TransferResult(calculator.get_axes_reco(), calculator.get_axes_gen()) {}

        Res4TransferResult(const TransferContainer& dipole_,
                           const TransferContainer& tee_,
                           const TransferContainer& triangle_) noexcept :
            dipole(dipole_), tee(tee_), triangle(triangle_),
            pt_denom_set(false),
            pt_denom_reco(-1), pt_denom_gen(-1) {}

        Res4TransferResult(TransferContainer&& dipole_,
                           TransferContainer&& tee_,
                           TransferContainer&& triangle_) noexcept :
            dipole(std::move(dipole_)), tee(std::move(tee_)), triangle(std::move(triangle_)),
            pt_denom_set(false),
            pt_denom_reco(-1), pt_denom_gen(-1) {}
        
        void fill_dipole(T iR_reco, T ir_reco, T ic_reco,
                         T iR_gen, T ir_gen, T ic_gen,
                         double wt_reco, double wt_gen) noexcept {
            dipole.fill(iR_reco, ir_reco, ic_reco, iR_gen, ir_gen, ic_gen, wt_reco, wt_gen);
        }

        void fill_tee(T iR_reco, T ir_reco, T ic_reco,
                      T iR_gen, T ir_gen, T ic_gen,
                      double wt_reco, double wt_gen) noexcept {
            tee.fill(iR_reco, ir_reco, ic_reco, iR_gen, ir_gen, ic_gen, wt_reco, wt_gen);
        }

        void fill_triangle(T iR_reco, T ir_reco, T ic_reco,
                           T iR_gen, T ir_gen, T ic_gen,
                           double wt_reco, double wt_gen) noexcept {
            triangle.fill(iR_reco, ir_reco, ic_reco, iR_gen, ir_gen, ic_gen, wt_reco, wt_gen);
        }

        void set_pt_denom(double pt_denom_reco_, double pt_denom_gen_) {
            if (pt_denom_set){
                throw std::runtime_error("pt_denom already set");
            } else {
                pt_denom_reco = pt_denom_reco_;
                pt_denom_gen = pt_denom_gen_;
                pt_denom_set = true;
            }
        }

        double get_pt_denom_reco() const {
            if (!pt_denom_set){
                throw std::runtime_error("pt_denom not set");
            } else {
                return pt_denom_reco;
            }
        }

        double get_pt_denom_gen() const {
            if (!pt_denom_set){
                throw std::runtime_error("pt_denom not set");
            } else {
                return pt_denom_gen;
            }
        }

        const TransferContainer& get_dipole() const  noexcept {
            return dipole;
        }

        const TransferContainer& get_tee() const  noexcept {
            return tee;
        }

        const TransferContainer& get_triangle() const  noexcept {
            return triangle;
        }

        Res4TransferResult<TransferContainer>& operator+=(const Res4TransferResult<TransferContainer>& other) noexcept {
            dipole += other.dipole;
            tee += other.tee;
            triangle += other.triangle;
            return *this;
        }

        Res4Result_MultiArray get_sum_over_gen() const noexcept {
            Res4Result_MultiArray sum(
                    std::move(dipole.get_sum_over_gen()),
                    std::move(tee.get_sum_over_gen()),
                    std::move(triangle.get_sum_over_gen()));
            sum.set_pt_denom(pt_denom_reco);
            return sum;
        }

        Res4Result_MultiArray get_sum_over_reco() const noexcept {
            Res4Result_MultiArray sum(
                    std::move(dipole.get_sum_over_reco()),
                    std::move(tee.get_sum_over_reco()),
                    std::move(triangle.get_sum_over_reco()));
            sum.set_pt_denom(pt_denom_gen);
            return sum;
        }

        template <class OtherContainer>
        bool operator==(const Res4TransferResult<OtherContainer>& other) const noexcept{
            return dipole == other.get_dipole() &&
                   tee == other.get_tee() &&
                   triangle == other.get_triangle();
        }

        double total_dipole_weight_gen() const noexcept {
            return dipole.total_weight_gen();
        }

        double total_tee_weight_gen() const noexcept {
            return tee.total_weight_gen();
        }

        double total_triangle_weight_gen() const noexcept {
            return triangle.total_weight_gen();
        }

        double total_dipole_weight_reco() const noexcept {
            return dipole.total_weight_reco();
        }

        double total_tee_weight_reco() const noexcept {
            return tee.total_weight_reco();
        }

        double total_triangle_weight_reco() const noexcept {
            return triangle.total_weight_reco();
        }
    private:
        TransferContainer dipole;
        TransferContainer tee;
        TransferContainer triangle;

        bool pt_denom_set;
        double pt_denom_reco, pt_denom_gen;
    };

    using Res4TransferResult_MultiArray = Res4TransferResult<ResTransferMultiArrayContainer>;
    using Res4TransferResult_Vector = Res4TransferResult<ResTransferVectorContainer<unsigned>>;
    using Res4TransferResult_Unbinned = Res4TransferResult<ResTransferVectorContainer<double>>;
};

#endif
