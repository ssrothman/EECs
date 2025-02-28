#ifndef SROTHMAN_EECS_Res3TransferResult_H
#define SROTHMAN_EECS_Res3TransferResult_H

#include "usings.h"
#include "Res3Axes.h"

#include "SRothman/SimonTools/src/histutil.h"
#include "ResTransferResultContainers.h"


namespace EEC{
    template <class TransferContainer>
    class Res3TransferResult{
    public:
        Res3TransferResult() noexcept  :
            data(),
            pt_denom_set(false),
            pt_denom_reco(-1), pt_denom_gen(-1) {}

        Res3TransferResult(
                const size_t nR_reco,
                const size_t nr_reco,
                const size_t nc_reco,
                const size_t nR_gen,
                const size_t nr_gen,
                const size_t nc_gen) noexcept  :
            data(nR_reco, nr_reco, nc_reco,
                     nR_gen, nr_gen, nc_gen),
            pt_denom_set(false),
            pt_denom_reco(-1), pt_denom_gen(-1) {}

        Res3TransferResult(
                const Res3Axes& axes_reco,
                const Res3Axes& axes_gen) noexcept  :
            Res3TransferResult(
                    simon::AXextent(axes_reco.R),
                    simon::AXextent(axes_reco.r),
                    simon::AXextent(axes_reco.c),
                    simon::AXextent(axes_gen.R),
                    simon::AXextent(axes_gen.r),
                    simon::AXextent(axes_gen.c)) {}

        template <typename CALCULATOR>
        Res3TransferResult(const CALCULATOR& calculator) noexcept:
            Res3TransferResult(calculator.get_axes_reco(), calculator.get_axes_gen()) {}
        
        void fill(unsigned iR_reco, unsigned ir_reco, unsigned ic_reco,
                           unsigned iR_gen, unsigned ir_gen, unsigned ic_gen,
                           double wt_reco, double wt_gen) noexcept {
            data.fill(iR_reco, ir_reco, ic_reco, iR_gen, ir_gen, ic_gen, wt_reco, wt_gen);
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

        const TransferContainer& get_data() const  noexcept {
            return data;
        }

        Res3TransferResult<TransferContainer>& operator+=(const Res3TransferResult<TransferContainer>& other) noexcept {
            data += other.data;
            return *this;
        }

        double total_weight_gen() const noexcept {
            return data.total_weight_gen();
        }

        double total_weight_reco() const noexcept {
            return data.total_weight_reco();
        }
    private:
        TransferContainer data;

        bool pt_denom_set;
        double pt_denom_reco, pt_denom_gen;
    };

    using Res3TransferResult_MultiArray = Res3TransferResult<ResTransferMultiArrayContainer>;
    using Res3TransferResult_Vector = Res3TransferResult<ResTransferVectorContainer>;
};

#endif
