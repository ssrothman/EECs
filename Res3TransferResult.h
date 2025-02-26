#ifndef SROTHMAN_EECS_Res3TransferResult_H
#define SROTHMAN_EECS_Res3TransferResult_H

#include "usings.h"
#include "Res3Axes.h"
#include "Res3Calculator.h"

#include "SRothman/SimonTools/src/histutil.h"


namespace EEC{

    class Res3TransferMultiArrayContainer{
    public:
        using data_t = multi_array<double, 6>;

        Res3TransferMultiArrayContainer(
                const size_t nR_reco, 
                const size_t nr_reco, 
                const size_t nc_reco,
                const size_t nR_gen, 
                const size_t nr_gen, 
                const size_t nc_gen) noexcept :
            nR_reco(nR_reco), nr_reco(nr_reco), nc_reco(nc_reco),
            nR_gen(nR_gen), nr_gen(nr_gen), nc_gen(nc_gen) {

            data.resize(boost::extents[nR_reco][nr_reco][nc_reco][nR_gen][nr_gen][nc_gen]);
            std::fill(data.data(), data.data() + data.num_elements(), 0);
        }

        Res3TransferMultiArrayContainer() noexcept :
            Res3TransferMultiArrayContainer(0, 0, 0, 0, 0, 0) {}

        void fill(unsigned iR_reco, unsigned ir_reco, unsigned ic_reco,
                  unsigned iR_gen, unsigned ir_gen, unsigned ic_gen,
                  double wt_reco, [[maybe_unused]] double wt_gen) noexcept {
            data[iR_reco][ir_reco][ic_reco][iR_gen][ir_gen][ic_gen] += wt_reco;
        }

        const data_t& get_data() const noexcept {
            return data;
        }

        Res3TransferMultiArrayContainer& operator+=(const Res3TransferMultiArrayContainer& other) noexcept {
            assert(data.num_elements() == other.data.num_elements());
            std::transform(data.data(), data.data() + data.num_elements(),
                           other.data.data(), data.data(),
                           std::plus<double>());
            return *this;
        }

        double total_weight() const noexcept {
            return std::accumulate(data.data(), data.data() + data.num_elements(), 0.0);
        }
        const size_t nR_reco, nr_reco, nc_reco;
        const size_t nR_gen, nr_gen, nc_gen;
    private:
        data_t data;
    };

    class Res3TransferVectorContainer{
    public:
        struct entry{
            unsigned iR_reco;
            unsigned ir_reco;
            unsigned ic_reco;
            unsigned iR_gen;
            unsigned ir_gen;
            unsigned ic_gen;
            double wt_reco;
            double wt_gen;
            entry(unsigned iR_reco, unsigned ir_reco, unsigned ic_reco,
                  unsigned iR_gen, unsigned ir_gen, unsigned ic_gen,
                  double wt_reco, double wt_gen) noexcept  :
                iR_reco(iR_reco), ir_reco(ir_reco), ic_reco(ic_reco),
                iR_gen(iR_gen), ir_gen(ir_gen), ic_gen(ic_gen),
                wt_reco(wt_reco), wt_gen(wt_gen) {}
        };
        using data_t = std::vector<entry>;

        Res3TransferVectorContainer(
                const size_t nR_reco, 
                const size_t nr_reco, 
                const size_t nc_reco,
                const size_t nR_gen, 
                const size_t nr_gen, 
                const size_t nc_gen) noexcept :
            nR_reco(nR_reco), nr_reco(nr_reco), nc_reco(nc_reco),
            nR_gen(nR_gen), nr_gen(nr_gen), nc_gen(nc_gen) {

            data.reserve(nR_reco);
        }

        Res3TransferVectorContainer() noexcept :
            Res3TransferVectorContainer(0, 0, 0, 0, 0, 0) {}

        void fill(unsigned iR_reco, unsigned ir_reco, unsigned ic_reco,
                  unsigned iR_gen, unsigned ir_gen, unsigned ic_gen,
                  double wt_reco, double wt_gen) noexcept {
            data.emplace_back(iR_reco, ir_reco, ic_reco, iR_gen, ir_gen, ic_gen, wt_reco, wt_gen);
        }

        const data_t& get_data() const noexcept {
            return data;
        }

        Res3TransferVectorContainer& operator+=(const Res3TransferVectorContainer& other) noexcept {
            data.insert(data.end(), other.data.begin(), other.data.end());
            return *this;
        }

        double total_weight_reco() const noexcept {
            return std::accumulate(data.begin(), data.end(), 0.0,
                                   [](double sum, const entry& e) { return sum + e.wt_reco; });
        }

        double total_weight_gen() const noexcept {
            return std::accumulate(data.begin(), data.end(), 0.0,
                                   [](double sum, const entry& e) { return sum + e.wt_gen; });
        }
        const size_t nR_reco, nr_reco, nc_reco;
        const size_t nR_gen, nr_gen, nc_gen;
    private:
        data_t data;
    };

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

        Res3TransferResult(const Res3TransferCalculator& calculator) noexcept:
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

    using Res3TransferResult_MultiArray = Res3TransferResult<Res3TransferMultiArrayContainer>;
    using Res3TransferResult_Vector = Res3TransferResult<Res3TransferVectorContainer>;
};

#endif
