#ifndef SROTHMAN_EEC_PROJTRANSFERRESULT_H
#define SROTHMAN_EEC_PROJTRANSFERRESULT_H

#include "usings.h"
#include "ProjAxes.h"

#include "SRothman/SimonTools/src/histutil.h"
#include "ProjTransferResultContainers.h"
#include "ProjResult.h"

namespace EEC{
    template <class TransferContainer>
    class ProjTransferResult{
    public:
        ProjTransferResult() noexcept :
            data(),
            pt_denom_set(false),
            pt_denom_reco(-1), pt_denom_gen(-1) {}

        using T = typename TransferContainer::TYPE;
        constexpr static bool SHOULD_BIN = TransferContainer::SHOULD_BIN;

        ProjTransferResult(
                const size_t nR_reco,
                const size_t nR_gen) noexcept :
            data({{TransferContainer(nR_reco, nR_gen),
                   TransferContainer(nR_reco, nR_gen),
                   TransferContainer(nR_reco, nR_gen),
                   TransferContainer(nR_reco, nR_gen),
                   TransferContainer(nR_reco, nR_gen)}}),
            pt_denom_set(false),
            pt_denom_reco(-1), pt_denom_gen(-1) {}

        ProjTransferResult(
                const ProjAxes& axes_reco,
                const ProjAxes& axes_gen) noexcept :
            ProjTransferResult(
                    simon::AXextent(axes_reco.R),
                    simon::AXextent(axes_gen.R)) {}

        template <typename CALCULATOR>
        ProjTransferResult(const CALCULATOR& calculator) noexcept :
            ProjTransferResult(calculator.get_axes_reco(), 
                               calculator.get_axes_gen()) {}

        template <unsigned ORDER>
        void fill(T iR_reco, T iR_gen,
                  double wt_reco, double wt_gen) noexcept {
            data[ORDER-2].fill(iR_reco, iR_gen, wt_reco, wt_gen);
        }

        void set_pt_denom(double pt_denom_reco_, double pt_denom_gen_) {
            if (pt_denom_set){
                throw std::runtime_error("pt_denom already set");
            } else {
                pt_denom_set = true;
                pt_denom_reco = pt_denom_reco_;
                pt_denom_gen = pt_denom_gen_;
            }
        }

        double get_pt_denom_reco() const {
            if (!pt_denom_set) {
                throw std::runtime_error("pt_denom_reco not set");
            }
            return pt_denom_reco;
        }

        double get_pt_denom_gen() const {
            if (!pt_denom_set) {
                throw std::runtime_error("pt_denom_gen not set");
            }
            return pt_denom_gen;
        }

        const std::array<TransferContainer,5>& get_data() const noexcept {
            return data;
        }

        double total_weight_reco(unsigned order) const noexcept {
            return data[order-2].total_weight_reco();
        }

        double total_weight_gen(unsigned order) const noexcept {
            return data[order-2].total_weight_gen();
        }

        ProjResult_Array get_sum_over_gen() const noexcept {
            ProjResult_Array sum({{std::move(data[0].get_sum_over_gen())}},
                                 {{std::move(data[1].get_sum_over_gen())}},
                                 {{std::move(data[2].get_sum_over_gen())}},
                                 {{std::move(data[3].get_sum_over_gen())}},
                                 {{std::move(data[4].get_sum_over_gen())}});
            sum.set_pt_denom(get_pt_denom_reco());
            return sum;
        }

        ProjResult_Array get_sum_over_reco() const noexcept {
            ProjResult_Array sum({{std::move(data[0].get_sum_over_reco())}},
                                 {{std::move(data[1].get_sum_over_reco())}},
                                 {{std::move(data[2].get_sum_over_reco())}},
                                 {{std::move(data[3].get_sum_over_reco())}},
                                 {{std::move(data[4].get_sum_over_reco())}});
            sum.set_pt_denom(get_pt_denom_gen());
            return sum;
        }
        
        template <typename OtherContainer>
        ProjTransferResult<TransferContainer>& operator+=(
                const ProjTransferResult<OtherContainer>& other) noexcept {
            data[0] += other.get_data()[0];
            data[1] += other.get_data()[1];
            data[2] += other.get_data()[2];
            data[3] += other.get_data()[3];
            data[4] += other.get_data()[4];
            return *this;
        }

        template <typename OtherContainer>
        ProjTransferResult<ProjArrayContainer> operator+(
                const ProjTransferResult<OtherContainer>& other) const noexcept {
            ProjTransferResult<ProjArrayContainer> result(*this);
            result += other;
            return result;
        }

        template <typename OtherContainer>
        ProjTransferResult<TransferContainer>& operator-=(
                const ProjTransferResult<OtherContainer>& other) noexcept {
            data[0] -= other.get_data()[0];  
            data[1] -= other.get_data()[1];
            data[2] -= other.get_data()[2];
            data[3] -= other.get_data()[3];
            data[4] -= other.get_data()[4];
            return *this;
        }

        template <typename OtherContainer>
        ProjTransferResult<ProjArrayContainer> operator-(
                const ProjTransferResult<OtherContainer>& other) const noexcept {
            ProjTransferResult<ProjArrayContainer> result(*this);
            result -= other;
            return result;
        }

        template <typename OtherContainer>
        bool operator==(const ProjTransferResult<OtherContainer>& other) const noexcept {
            return data[0] == other.get_data()[0] &&
                   data[1] == other.get_data()[1] &&
                   data[2] == other.get_data()[2] &&
                   data[3] == other.get_data()[3] &&
                   data[4] == other.get_data()[4];
        }

    private:
        std::array<TransferContainer, 5> data;

        bool pt_denom_set;
        double pt_denom_reco;
        double pt_denom_gen;
    };

    using ProjTransferResult_Array = ProjTransferResult<ProjTransferArrayContainer>;
    using ProjTransferResult_Vector = ProjTransferResult<ProjTransferVectorContainer<unsigned>>;
    using ProjTransferResult_Unbinned = ProjTransferResult<ProjTransferVectorContainer<double>>;
};

#endif
