#ifndef SROTHMAN_EEC_RES_TRANSFER_CONTAINERS_H
#define SROTHMAN_EEC_RES_TRANSFER_CONTAINERS_H

#include "usings.h"

#include "SRothman/SimonTools/src/histutil.h"

namespace EEC{
    class ResTransferMultiArrayContainer{
    public:
        using data_t = multi_array<double, 6>;

        ResTransferMultiArrayContainer(
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

        ResTransferMultiArrayContainer() noexcept :
            ResTransferMultiArrayContainer(0, 0, 0, 0, 0, 0) {}

        void fill(unsigned iR_reco, unsigned ir_reco, unsigned ic_reco,
                  unsigned iR_gen, unsigned ir_gen, unsigned ic_gen,
                  double wt_reco, [[maybe_unused]] double wt_gen) noexcept {
            data[iR_reco][ir_reco][ic_reco][iR_gen][ir_gen][ic_gen] += wt_reco;
        }

        const data_t& get_data() const noexcept {
            return data;
        }

        ResTransferMultiArrayContainer& operator+=(const ResTransferMultiArrayContainer& other) noexcept {
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

    class ResTransferVectorContainer{
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

        ResTransferVectorContainer(
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

        ResTransferVectorContainer() noexcept :
            ResTransferVectorContainer(0, 0, 0, 0, 0, 0) {}

        void fill(unsigned iR_reco, unsigned ir_reco, unsigned ic_reco,
                  unsigned iR_gen, unsigned ir_gen, unsigned ic_gen,
                  double wt_reco, double wt_gen) noexcept {
            data.emplace_back(iR_reco, ir_reco, ic_reco, iR_gen, ir_gen, ic_gen, wt_reco, wt_gen);
        }

        const data_t& get_data() const noexcept {
            return data;
        }

        ResTransferVectorContainer& operator+=(const ResTransferVectorContainer& other) noexcept {
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

   
};

#endif
