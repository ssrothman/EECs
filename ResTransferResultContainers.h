#ifndef SROTHMAN_EEC_RES_TRANSFER_CONTAINERS_H
#define SROTHMAN_EEC_RES_TRANSFER_CONTAINERS_H

#include "usings.h"
#include "ResResultContainers.h"

#include "SRothman/SimonTools/src/histutil.h"

namespace EEC{
    template <typename T>
    class ResTransferVectorContainer{
    public:
        struct entry{
            T iR_reco;
            T ir_reco;
            T ic_reco;
            T iR_gen;
            T ir_gen;
            T ic_gen;
            double wt_reco;
            double wt_gen;
            entry(T iR_reco, T ir_reco, T ic_reco,
                  T iR_gen, T ir_gen, T ic_gen,
                  double wt_reco, double wt_gen) noexcept  :
                iR_reco(iR_reco), ir_reco(ir_reco), ic_reco(ic_reco),
                iR_gen(iR_gen), ir_gen(ir_gen), ic_gen(ic_gen),
                wt_reco(wt_reco), wt_gen(wt_gen) {}
        };
        using data_t = std::vector<entry>;

        constexpr static bool SHOULD_BIN = std::is_same<T, unsigned>::value;

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

        void fill(T iR_reco, T ir_reco, T ic_reco,
                  T iR_gen, T ir_gen, T ic_gen,
                  double wt_reco, double wt_gen) noexcept {
            data.emplace_back(iR_reco, ir_reco, ic_reco, iR_gen, ir_gen, ic_gen, wt_reco, wt_gen);
        }

        const data_t& get_data() const noexcept {
            return data;
        }

        ResMultiArrayContainer get_sum_over_gen() const noexcept {
            ResMultiArrayContainer sum(nR_reco, nr_reco, nc_reco);
            for (const auto& entry : data){
                sum.fill(entry.iR_reco, entry.ir_reco, entry.ic_reco, entry.wt_reco);
            }
            return sum;
        }

        ResMultiArrayContainer get_sum_over_reco() const noexcept {
            ResMultiArrayContainer sum(nR_gen, nr_gen, nc_gen);
            for (const auto& entry : data){
                sum.fill(entry.iR_gen, entry.ir_gen, entry.ic_gen, entry.wt_gen);
            }
            return sum;
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

    class ResTransferMultiArrayContainer{
    public:
        using data_t = multi_array<double, 6>;

        static constexpr bool SHOULD_BIN = true;

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

        ResTransferMultiArrayContainer(const ResTransferVectorContainer<unsigned>& other) noexcept :
            ResTransferMultiArrayContainer(other.nR_reco, other.nr_reco, other.nc_reco,
                                           other.nR_gen, other.nr_gen, other.nc_gen) {
            for (const auto& entry : other.get_data()){
                fill(entry.iR_reco, entry.ir_reco, entry.ic_reco,
                     entry.iR_gen, entry.ir_gen, entry.ic_gen,
                     entry.wt_reco, entry.wt_gen);
            }
        }

        const data_t& get_data() const noexcept {
            return data;
        }

        ResMultiArrayContainer get_sum_over_gen() const noexcept {
            ResMultiArrayContainer sum(nR_reco, nr_reco, nc_reco);
            for (unsigned iR_reco = 0; iR_reco < nR_reco; ++iR_reco){
                for (unsigned ir_reco = 0; ir_reco < nr_reco; ++ir_reco){
                    for (unsigned ic_reco = 0; ic_reco < nc_reco; ++ic_reco){
                        for (unsigned iR_gen = 0; iR_gen < nR_gen; ++iR_gen){
                            for (unsigned ir_gen = 0; ir_gen < nr_gen; ++ir_gen){
                                for (unsigned ic_gen = 0; ic_gen < nc_gen; ++ic_gen){
                                    sum.fill(iR_reco, ir_reco, ic_reco, data[iR_reco][ir_reco][ic_reco][iR_gen][ir_gen][ic_gen]);
                                }
                            }
                        }
                    }
                }
            }
            return sum;
        }

        ResMultiArrayContainer get_sum_over_reco() const noexcept {
            ResMultiArrayContainer sum(nR_gen, nr_gen, nc_gen);
            for (unsigned iR_reco = 0; iR_reco < nR_reco; ++iR_reco){
                for (unsigned ir_reco = 0; ir_reco < nr_reco; ++ir_reco){
                    for (unsigned ic_reco = 0; ic_reco < nc_reco; ++ic_reco){
                        for (unsigned iR_gen = 0; iR_gen < nR_gen; ++iR_gen){
                            for (unsigned ir_gen = 0; ir_gen < nr_gen; ++ir_gen){
                                for (unsigned ic_gen = 0; ic_gen < nc_gen; ++ic_gen){
                                    sum.fill(iR_gen, ir_gen, ic_gen, data[iR_reco][ir_reco][ic_reco][iR_gen][ir_gen][ic_gen]);
                                }
                            }
                        }
                    }
                }
            }
            return sum;
        }

        ResTransferMultiArrayContainer& operator+=(const ResTransferMultiArrayContainer& other) noexcept {
            assert(data.num_elements() == other.data.num_elements());
            std::transform(data.data(), data.data() + data.num_elements(),
                           other.data.data(), data.data(),
                           std::plus<double>());
            return *this;
        }

        bool operator==(const ResTransferMultiArrayContainer& other) const noexcept {
            return std::equal(
                    data.data(), 
                    data.data() + data.num_elements(),
                    other.data.data(),
                    [](double a, double b) { 
                        return std::abs(a - b) < 1e-6; 
                    });
        }

        double total_weight_reco() const noexcept {
            return std::accumulate(data.data(), data.data() + data.num_elements(), 0.0);
        }

        double total_weight_gen() const noexcept {
            return std::accumulate(data.data(), data.data() + data.num_elements(), 0.0);
        }

        const size_t nR_reco, nr_reco, nc_reco;
        const size_t nR_gen, nr_gen, nc_gen;
    private:
        data_t data;
    };

    inline bool operator==(const ResTransferVectorContainer<unsigned>& a, const ResTransferVectorContainer<unsigned>& b) noexcept {
        return ResTransferMultiArrayContainer(a) == ResTransferMultiArrayContainer(b);
    }

    inline bool operator==(const ResTransferMultiArrayContainer& a, const ResTransferVectorContainer<unsigned>& b) noexcept {
        return a == ResTransferMultiArrayContainer(b);
    }

    inline bool operator==(const ResTransferVectorContainer<unsigned>& a, const ResTransferMultiArrayContainer& b) noexcept {
        return ResTransferMultiArrayContainer(a) == b;
    }

    inline void print_nonzero(const ResTransferMultiArrayContainer& a){
        for (size_t iR_reco = 0; iR_reco < a.nR_reco; ++iR_reco){
            for (size_t ir_reco = 0; ir_reco < a.nr_reco; ++ir_reco){
                for (size_t ic_reco = 0; ic_reco < a.nc_reco; ++ic_reco){
                    for (size_t iR_gen = 0; iR_gen < a.nR_gen; ++iR_gen){
                        for (size_t ir_gen = 0; ir_gen < a.nr_gen; ++ir_gen){
                            for (size_t ic_gen = 0; ic_gen < a.nc_gen; ++ic_gen){
                                double val = a.get_data()[iR_reco][ir_reco][ic_reco][iR_gen][ir_gen][ic_gen];
                                if (val > 0){
                                    printf("(%lu, %lu, %lu)_Gen -> (%lu, %lu, %lu)_Reco: %g\n", 
                                            iR_gen, ir_gen, ic_gen,
                                            iR_reco, ir_reco, ic_reco,
                                            a.get_data()[iR_reco][ir_reco][ic_reco][iR_gen][ir_gen][ic_gen]);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    inline void print_nonzero(const ResTransferVectorContainer<unsigned>& a){
        print_nonzero(ResTransferMultiArrayContainer(a));
    }
};

#endif
