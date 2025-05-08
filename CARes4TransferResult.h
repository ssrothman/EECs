#ifndef SROTHMAN_EEC_CA_RES4_TRANSFER_RESULT_H
#define SROTHMAN_EEC_CA_RES4_TRANSFER_RESULT_H

#include "usings.h"
#include "CARes4Axes.h"
#include "ResTransferResultContainers.h"

#include "SRothman/SimonTools/src/histutil.h"


namespace EEC{
    template <class TransferContainer>
    class CARes4TransferResult{
    public:
        CARes4TransferResult() noexcept  :
            chain(), 
            symmetric_wrtR(),
            symmetric_wrtr(),
            chain_to_symmetric_wrtR(),
            chain_to_symmetric_wrtr(),
            symmetric_wrtR_to_chain(),
            symmetric_wrtr_to_chain(),
            pt_denom_set(false),
            pt_denom_reco(-1), pt_denom_gen(-1) {}

        constexpr static bool SHOULD_BIN = TransferContainer::SHOULD_BIN;
        using T = typename TransferContainer::TYPE;

        CARes4TransferResult(
                const size_t nR_reco,
                const size_t nr_chain_reco,
                const size_t nc_chain_reco,
                const size_t nr_symmetric_reco,
                const size_t nc_symmetric_reco,
                const size_t nR_gen,
                const size_t nr_chain_gen,
                const size_t nc_chain_gen,
                const size_t nr_symmetric_gen,
                const size_t nc_symmetric_gen) noexcept  :
            chain(
                nR_reco, nr_chain_reco, nc_chain_reco,
                nR_gen, nr_chain_gen, nc_chain_gen
            ),
            symmetric_wrtR(
                nR_reco, nr_symmetric_reco, nc_symmetric_reco,
                nR_gen, nr_symmetric_gen, nc_symmetric_gen
            ),
            symmetric_wrtr(
                nR_reco, nr_symmetric_reco, nc_symmetric_reco,
                nR_gen, nr_symmetric_gen, nc_symmetric_gen
            ),
            chain_to_symmetric_wrtR(
                nR_reco, nr_symmetric_reco, nc_symmetric_reco,
                nR_gen, nr_chain_gen, nc_chain_gen
            ),
            chain_to_symmetric_wrtr(
                nR_reco, nr_symmetric_reco, nc_symmetric_reco,
                nR_gen, nr_chain_gen, nc_chain_gen
            ),
            symmetric_wrtR_to_chain(
                nR_reco, nr_chain_reco, nc_chain_reco,
                nR_gen, nr_symmetric_gen, nc_symmetric_gen
            ),
            symmetric_wrtr_to_chain(
                nR_reco, nr_chain_reco, nc_chain_reco,
                nR_gen, nr_symmetric_gen, nc_symmetric_gen
            ),

            pt_denom_set(false),
            pt_denom_reco(-1), pt_denom_gen(-1) {}

        CARes4TransferResult(
                const CARes4Axes& axes_reco,
                const CARes4Axes& axes_gen) noexcept  :
            CARes4TransferResult(
                    simon::AXextent(axes_reco.R),
                    simon::AXextent(axes_reco.r_chain),
                    simon::AXextent(axes_reco.c_chain),
                    simon::AXextent(axes_reco.r_symmetric),
                    simon::AXextent(axes_reco.c_symmetric),
                    simon::AXextent(axes_gen.R),
                    simon::AXextent(axes_gen.r_chain),
                    simon::AXextent(axes_gen.c_chain),
                    simon::AXextent(axes_gen.r_symmetric),
                    simon::AXextent(axes_gen.c_symmetric)) {}

        template <typename CALCULATOR>
        CARes4TransferResult(const CALCULATOR& calculator) noexcept:
            CARes4TransferResult(
                    calculator.get_axes_reco(), 
                    calculator.get_axes_gen()) {}
        
        CARes4TransferResult(const TransferContainer& chain_,
                            const TransferContainer& symmetric_wrtR_,
                            const TransferContainer& symmetric_wrtr_,
                            const TransferContainer& chain_to_symmetric_wrtR_,
                            const TransferContainer& chain_to_symmetric_wrtr_,
                            const TransferContainer& symmetric_wrtR_to_chain_,
                            const TransferContainer& symmetric_wrtr_to_chain_) noexcept :
            chain(chain_),
            symmetric_wrtR(symmetric_wrtR_),
            symmetric_wrtr(symmetric_wrtr_),
            chain_to_symmetric_wrtR(chain_to_symmetric_wrtR_),
            chain_to_symmetric_wrtr(chain_to_symmetric_wrtr_),
            symmetric_wrtR_to_chain(symmetric_wrtR_to_chain_),
            symmetric_wrtr_to_chain(symmetric_wrtr_to_chain_),
            pt_denom_set(false),
            pt_denom_reco(-1), pt_denom_gen(-1) {}

        CARes4TransferResult(TransferContainer&& chain_,
                            TransferContainer&& symmetric_wrtR_,
                            TransferContainer&& symmetric_wrtr_,
                            TransferContainer&& chain_to_symmetric_wrtR_,
                            TransferContainer&& chain_to_symmetric_wrtr_,
                            TransferContainer&& symmetric_wrtR_to_chain_,
                            TransferContainer&& symmetric_wrtr_to_chain_) noexcept :
            chain(std::move(chain_)),
            symmetric_wrtR(std::move(symmetric_wrtR_)),
            symmetric_wrtr(std::move(symmetric_wrtr_)),
            chain_to_symmetric_wrtR(std::move(chain_to_symmetric_wrtR_)),
            chain_to_symmetric_wrtr(std::move(chain_to_symmetric_wrtr_)),
            symmetric_wrtR_to_chain(std::move(symmetric_wrtR_to_chain_)),
            symmetric_wrtr_to_chain(std::move(symmetric_wrtr_to_chain_)),
            pt_denom_set(false),
            pt_denom_reco(-1), pt_denom_gen(-1) {}

        template <typename T>
        void fill_chain(const CAres4_entry<T>& entry_reco,
                        const CAres4_entry<T>& entry_gen,
                        double wt_reco, 
                        double wt_gen) noexcept {
#ifdef CHECK_BY_HAND
            printf("filling chain -> chain\n");
            printf("\tweight %g -> %g\n", wt_gen, wt_reco);
            //printf("\t(%u, %u, %u) -> (%u, %u, %u)\n", 
            //        entry_gen.R_idx, entry_gen.chain_r_idx, entry_gen.chain_c_idx,
            //        entry_reco.R_idx, entry_reco.chain_r_idx, entry_reco.chain_c_idx);
#endif
            chain.fill(
                entry_reco.R_idx, 
                entry_reco.chain_r_idx, 
                entry_reco.chain_c_idx, 
                entry_gen.R_idx, 
                entry_gen.chain_r_idx, 
                entry_gen.chain_c_idx, 
                wt_reco, wt_gen
            );
        }

        template <typename T>
        void fill_symmetric(const CAres4_entry<T>& entry_reco,
                            const CAres4_entry<T>& entry_gen,
                            double wt_reco,
                            double wt_gen) noexcept {
#ifdef CHECK_BY_HAND
            printf("filling symmetric -> symmetric\n");
            printf("\tweight %g -> %g\n", wt_gen, wt_reco);
            //printf("\twrt_r: (%u, %u, %u) -> (%u, %u, %u)\n", 
            //        entry_gen.R_idx,
            //        entry_gen.symmetric_wrtr_r_idx,
            //        entry_gen.symmetric_wrtr_c_idx,
            //        entry_reco.R_idx,
            //        entry_reco.symmetric_wrtr_r_idx,
            //        entry_reco.symmetric_wrtr_c_idx);
            //printf("\twrtR1: (%u, %u, %u) -> (%u, %u, %u)\n",
            //        entry_gen.R_idx,
            //        entry_gen.symmetric_wrtR1_r_idx,
            //        entry_gen.symmetric_wrtR1_c_idx,
            //        entry_reco.R_idx,
            //        entry_reco.symmetric_wrtR1_r_idx,
            //        entry_reco.symmetric_wrtR1_c_idx);
            //printf("\twrtR2: (%u, %u, %u) -> (%u, %u, %u)\n",
            //        entry_gen.R_idx,
            //        entry_gen.symmetric_wrtR2_r_idx,
            //        entry_gen.symmetric_wrtR2_c_idx,
            //        entry_reco.R_idx,
            //        entry_reco.symmetric_wrtR2_r_idx,
            //        entry_reco.symmetric_wrtR2_c_idx);
#endif
            symmetric_wrtR.fill(
                entry_reco.R_idx, 
                entry_reco.symmetric_wrtR1_r_idx, 
                entry_reco.symmetric_wrtR1_c_idx,
                entry_gen.R_idx, 
                entry_gen.symmetric_wrtR1_r_idx, 
                entry_gen.symmetric_wrtR1_c_idx,
                wt_reco, wt_gen
            );

            symmetric_wrtR.fill(
                entry_reco.R_idx, 
                entry_reco.symmetric_wrtR2_r_idx, 
                entry_reco.symmetric_wrtR2_c_idx,
                entry_gen.R_idx, 
                entry_gen.symmetric_wrtR2_r_idx, 
                entry_gen.symmetric_wrtR2_c_idx,
                wt_reco, wt_gen
            );

            symmetric_wrtr.fill(
                entry_reco.R_idx, 
                entry_reco.symmetric_wrtr_r_idx, 
                entry_reco.symmetric_wrtr_c_idx,
                entry_gen.R_idx, 
                entry_gen.symmetric_wrtr_r_idx, 
                entry_gen.symmetric_wrtr_c_idx,
                wt_reco, wt_gen
            );
        }

        template <typename T>
        void fill_chain_to_symmetric(
                const CAres4_entry<T>& entry_reco,
                const CAres4_entry<T>& entry_gen,
                double wt_reco,
                double wt_gen) noexcept {
            
            chain_to_symmetric_wrtR.fill(
                entry_reco.R_idx,
                entry_reco.symmetric_wrtR1_r_idx,
                entry_reco.symmetric_wrtR1_c_idx,
                entry_gen.R_idx,
                entry_gen.chain_r_idx,
                entry_gen.chain_c_idx,
                wt_reco, wt_gen
            );

            chain_to_symmetric_wrtR.fill(
                entry_reco.R_idx,
                entry_reco.symmetric_wrtR2_r_idx,
                entry_reco.symmetric_wrtR2_c_idx,
                entry_gen.R_idx,
                entry_gen.chain_r_idx,
                entry_gen.chain_c_idx,
                wt_reco, wt_gen
            );

            chain_to_symmetric_wrtr.fill(
                entry_reco.R_idx,
                entry_reco.symmetric_wrtr_r_idx,
                entry_reco.symmetric_wrtr_c_idx,
                entry_gen.R_idx,
                entry_gen.chain_r_idx,
                entry_gen.chain_c_idx,
                wt_reco, wt_gen
            );
        }

        template <typename T>
        void fill_symmetric_to_chain(
                const CAres4_entry<T>& entry_reco,
                const CAres4_entry<T>& entry_gen,
                double wt_reco, 
                double wt_gen) noexcept {
            
            symmetric_wrtR_to_chain.fill(
                entry_reco.R_idx,
                entry_reco.chain_r_idx,
                entry_reco.chain_c_idx,
                entry_gen.R_idx,
                entry_gen.symmetric_wrtR1_r_idx,
                entry_gen.symmetric_wrtR1_c_idx,
                wt_reco/2, wt_gen
            );

            symmetric_wrtR_to_chain.fill(
                entry_reco.R_idx,
                entry_reco.chain_r_idx,
                entry_reco.chain_c_idx,
                entry_gen.R_idx,
                entry_gen.symmetric_wrtR2_r_idx,
                entry_gen.symmetric_wrtR2_c_idx,
                wt_reco/2, wt_gen
            );

            symmetric_wrtr_to_chain.fill(
                entry_reco.R_idx,
                entry_reco.chain_r_idx,
                entry_reco.chain_c_idx,
                entry_gen.R_idx,
                entry_gen.symmetric_wrtr_r_idx,
                entry_gen.symmetric_wrtr_c_idx,
                wt_reco, wt_gen
            );
        }

        template <typename T>
        void fill(
                const CAres4_entry<T>& entry_reco,
                const CAres4_entry<T>& entry_gen,
                double wt_reco,
                double wt_gen) noexcept {
            if (entry_gen.is_chain){
                if(entry_reco.is_chain){
                    fill_chain(
                        entry_reco, entry_gen, 
                        wt_reco, wt_gen
                    );
                } else {
                    fill_chain_to_symmetric(
                        entry_reco, entry_gen, 
                        wt_reco, wt_gen
                    );
                }
            } else {
                if(entry_reco.is_chain){
                    fill_symmetric_to_chain(
                        entry_reco, entry_gen, 
                        wt_reco, wt_gen
                    );
                } else {
                    fill_symmetric(
                        entry_reco, entry_gen, 
                        wt_reco, wt_gen
                    );
                }
            }
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

        const TransferContainer& get_chain() const  noexcept {
            return chain;
        }

        const TransferContainer& get_symmetric_wrtR() const  noexcept {
            return symmetric_wrtR;
        }

        const TransferContainer& get_symmetric_wrtr() const  noexcept {
            return symmetric_wrtr;
        }

        const TransferContainer& get_chain_to_symmetric_wrtR() const  noexcept {
            return chain_to_symmetric_wrtR;
        }

        const TransferContainer& get_chain_to_symmetric_wrtr() const  noexcept {
            return chain_to_symmetric_wrtr;
        }

        const TransferContainer& get_symmetric_wrtR_to_chain() const  noexcept {
            return symmetric_wrtR_to_chain;
        }

        const TransferContainer& get_symmetric_wrtr_to_chain() const  noexcept {
            return symmetric_wrtr_to_chain;
        }

        CARes4TransferResult<TransferContainer>& operator+=(const CARes4TransferResult<TransferContainer>& other) noexcept {
            chain += other.chain;
            symmetric_wrtR += other.symmetric_wrtR;
            symmetric_wrtr += other.symmetric_wrtr;
            chain_to_symmetric_wrtR += other.chain_to_symmetric_wrtR;
            chain_to_symmetric_wrtr += other.chain_to_symmetric_wrtr;
            symmetric_wrtR_to_chain += other.symmetric_wrtR_to_chain;
            symmetric_wrtr_to_chain += other.symmetric_wrtr_to_chain;
            return *this;
        }

        CARes4Result_MultiArray get_sum_over_gen() const noexcept {
            CARes4Result_MultiArray sum(
                    std::move(chain.get_sum_over_gen() + symmetric_wrtR_to_chain.get_sum_over_gen()),
                    std::move(symmetric_wrtR.get_sum_over_gen() + chain_to_symmetric_wrtR.get_sum_over_gen()),
                    std::move(symmetric_wrtr.get_sum_over_gen() + chain_to_symmetric_wrtr.get_sum_over_gen()));
            sum.set_pt_denom(pt_denom_reco);
            return sum;
        }

        CARes4Result_MultiArray get_sum_over_reco() const noexcept {
            CARes4Result_MultiArray sum(
                    std::move(chain.get_sum_over_reco() + chain_to_symmetric_wrtr.get_sum_over_reco()),
                    std::move(symmetric_wrtR.get_sum_over_reco() + symmetric_wrtR_to_chain.get_sum_over_reco()),
                    std::move(symmetric_wrtr.get_sum_over_reco() + symmetric_wrtr_to_chain.get_sum_over_reco()));
            sum.set_pt_denom(pt_denom_gen);
            return sum;
        }

        template <class OtherContainer>
        bool operator==(const CARes4TransferResult<OtherContainer>& other) const noexcept{
            return chain == other.get_chain() &&
                   symmetric_wrtR == other.get_symmetric_wrtR() &&
                   symmetric_wrtr == other.get_symmetric_wrtr() &&
                   chain_to_symmetric_wrtR == other.get_chain_to_symmetric_wrtR() &&
                   chain_to_symmetric_wrtr == other.get_chain_to_symmetric_wrtr() &&
                   symmetric_wrtR_to_chain == other.get_symmetric_wrtR_to_chain();
                   symmetric_wrtr_to_chain == other.get_symmetric_wrtr_to_chain();
        }

        double total_chain_weight_gen() const noexcept {
            return chain.total_weight_gen();
        }

        double total_symmetric_wrtR_weight_gen() const noexcept {
            return symmetric_wrtR.total_weight_gen();
        }

        double total_symmetric_wrtr_weight_gen() const noexcept {
            return symmetric_wrtr.total_weight_gen();
        }

        double total_chain_to_symmetric_wrtR_weight_gen() const noexcept {
            return chain_to_symmetric_wrtR.total_weight_gen();
        }

        double total_chain_to_symmetric_wrtr_weight_gen() const noexcept {
            return chain_to_symmetric_wrtr.total_weight_gen();
        }

        double total_symmetric_wrtR_to_chain_weight_gen() const noexcept {
            return symmetric_wrtR_to_chain.total_weight_gen();
        }

        double total_symmetric_wrtr_to_chain_weight_gen() const noexcept {
            return symmetric_wrtr_to_chain.total_weight_gen();
        }

        double total_chain_weight_reco() const noexcept {
            return chain.total_weight_reco();
        }

        double total_symmetric_wrtR_weight_reco() const noexcept {
            return symmetric_wrtR.total_weight_reco();
        }

        double total_symmetric_wrtr_weight_reco() const noexcept {
            return symmetric_wrtr.total_weight_reco();
        }

        double total_chain_to_symmetric_wrtR_weight_reco() const noexcept {
            return chain_to_symmetric_wrtR.total_weight_reco();
        }

        double total_chain_to_symmetric_wrtr_weight_reco() const noexcept {
            return chain_to_symmetric_wrtr.total_weight_reco();
        }

        double total_symmetric_wrtR_to_chain_weight_reco() const noexcept {
            return symmetric_wrtR_to_chain.total_weight_reco();
        }

        double total_symmetric_wrtr_to_chain_weight_reco() const noexcept {
            return symmetric_wrtr_to_chain.total_weight_reco();
        }
    private:
        TransferContainer chain;
        TransferContainer symmetric_wrtR;
        TransferContainer symmetric_wrtr;

        TransferContainer chain_to_symmetric_wrtR;
        TransferContainer chain_to_symmetric_wrtr;

        TransferContainer symmetric_wrtR_to_chain;
        TransferContainer symmetric_wrtr_to_chain;

        bool pt_denom_set;
        double pt_denom_reco, pt_denom_gen;
    };

    using CARes4TransferResult_MultiArray = CARes4TransferResult<ResTransferMultiArrayContainer>;
    using CARes4TransferResult_Vector = CARes4TransferResult<ResTransferVectorContainer<unsigned>>;
    using CARes4TransferResult_Unbinned = CARes4TransferResult<ResTransferVectorContainer<double>>;
};

#endif
