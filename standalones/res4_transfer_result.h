#ifndef SROTHMAN_EECS_RES4_TRANSFER_RESULT_H
#define SROTHMAN_EECS_RES4_TRANSFER_RESULT_H

#include "usings.h"
#include "boost/multi_array.hpp"
#include <vector>

namespace standaloneEEC{
    class res4_transfer_multi_array_container{
    public:
        using data_t = multi_array<double, 6>;

        res4_transfer_multi_array_container(
                const size_t nR_1, 
                const size_t nr_1, 
                const size_t nc_1,
                const size_t nR_2, 
                const size_t nr_2, 
                const size_t nc_2) noexcept {
            data.resize(boost::extents[nR_1][nr_1][nc_1][nR_2][nr_2][nc_2]);
            std::fill(data.data(), data.data() + data.num_elements(), 0);
        }

        void fill(unsigned iR1, unsigned ir1, unsigned ic1,
                  unsigned iR2, unsigned ir2, unsigned ic2,
                  double wt) noexcept {
            data[iR1][ir1][ic1][iR2][ir2][ic2] += wt;
        }

        const data_t& get_data() const noexcept {
            return data;
        }

        res4_transfer_multi_array_container& operator+=(const res4_transfer_multi_array_container& other) noexcept {
            assert(data.num_elements() == other.data.num_elements());
            std::transform(data.data(), data.data() + data.num_elements(),
                           other.data.data(), data.data(),
                           std::plus<double>());
            return *this;
        }

        double total_weight() const noexcept {
            return std::accumulate(data.data(), data.data() + data.num_elements(), 0.0);
        }
    private:
        data_t data;
    };

    class res4_transfer_vector_container{
    public:
        struct entry{
            unsigned iR1;
            unsigned ir1;
            unsigned ic1;
            unsigned iR2;
            unsigned ir2;
            unsigned ic2;
            double wt;
            entry(unsigned iR1, unsigned ir1, unsigned ic1,
                  unsigned iR2, unsigned ir2, unsigned ic2,
                  double wt) noexcept  :
                iR1(iR1), ir1(ir1), ic1(ic1),
                iR2(iR2), ir2(ir2), ic2(ic2), wt(wt) {}
        };
        using data_t = std::vector<entry>;

        res4_transfer_vector_container(
                [[maybe_unused]] const size_t nR_1, 
                [[maybe_unused]] const size_t nr_1, 
                [[maybe_unused]] const size_t nc_1,
                [[maybe_unused]] const size_t nR_2, 
                [[maybe_unused]] const size_t nr_2, 
                [[maybe_unused]] const size_t nc_2) noexcept  {
            data.reserve(100);
        }

        void fill(unsigned iR1, unsigned ir1, unsigned ic1,
                  unsigned iR2, unsigned ir2, unsigned ic2,
                  double wt) noexcept {
            data.emplace_back(iR1, ir1, ic1, iR2, ir2, ic2, wt);
        }

        const data_t& get_data() const noexcept {
            return data;
        }

        res4_transfer_vector_container& operator+=(const res4_transfer_vector_container& other) noexcept {
            data.insert(data.end(), other.data.begin(), other.data.end());
            return *this;
        }

        double total_weight() const noexcept {
            return std::accumulate(data.begin(), data.end(), 0.0,
                                   [](double sum, const entry& e) { return sum + e.wt; });
        }
    private:
        data_t data;
    };

    template <class TransferContainer, class BasicContainer>
    class res4_transfer_result{
    public:
        res4_result<BasicContainer> unmatched1;
        res4_result<BasicContainer> unmatched2;

        res4_transfer_result(
                const size_t nR_1,
                const size_t nr_dipole_1,
                const size_t nc_dipole_1,
                const size_t nr_tee_1,
                const size_t nc_tee_1,
                const size_t nr_triangle_1,
                const size_t nc_triangle_1,
                const size_t nR_2,
                const size_t nr_dipole_2,
                const size_t nc_dipole_2,
                const size_t nr_tee_2,
                const size_t nc_tee_2,
                const size_t nr_triangle_2,
                const size_t nc_triangle_2) noexcept  :
            unmatched1(nR_1, nr_dipole_1, nc_dipole_1,
                       nr_tee_1, nc_tee_1,
                       nr_triangle_1, nc_triangle_1),
            unmatched2(nR_2, nr_dipole_2, nc_dipole_2,
                       nr_tee_2, nc_tee_2,
                       nr_triangle_2, nc_triangle_2),
            dipole(nR_1, nr_dipole_1, nc_dipole_1, 
                    nR_2, nr_dipole_2, nc_dipole_2),
            tee(nR_1, nr_tee_1, nc_tee_1,
                nR_2, nr_tee_2, nc_tee_2),
            triangle(nR_1, nr_triangle_1, nc_triangle_1,
                     nR_2, nr_triangle_2, nc_triangle_2) {}
        
        void fill_dipole(unsigned iR1, unsigned ir1, unsigned ic1,
                         unsigned iR2, unsigned ir2, unsigned ic2,
                         double wt) noexcept {
            dipole.fill(iR1, ir1, ic1, iR2, ir2, ic2, wt);
        }

        void fill_tee(unsigned iR1, unsigned ir1, unsigned ic1,
                      unsigned iR2, unsigned ir2, unsigned ic2,
                      double wt) noexcept {
            tee.fill(iR1, ir1, ic1, iR2, ir2, ic2, wt);
        }

        void fill_triangle(unsigned iR1, unsigned ir1, unsigned ic1,
                           unsigned iR2, unsigned ir2, unsigned ic2,
                           double wt) noexcept {
            triangle.fill(iR1, ir1, ic1, iR2, ir2, ic2, wt);
        }

        const typename TransferContainer::data_t& get_dipole_data() const  noexcept {
            return dipole.get_data();
        }

        const typename TransferContainer::data_t& get_tee_data() const  noexcept {
            return tee.get_data();
        }

        const typename TransferContainer::data_t& get_triangle_data() const  noexcept {
            return triangle.get_data();
        }

        res4_transfer_result<TransferContainer, BasicContainer>& operator+=(const res4_transfer_result<TransferContainer, BasicContainer>& other) noexcept {
            unmatched1 += other.unmatched1;
            unmatched2 += other.unmatched2;
            dipole += other.dipole;
            tee += other.tee;
            triangle += other.triangle;
            return *this;
        }

        double total_dipole_weight() const noexcept {
            return dipole.total_weight();
        }

        double total_tee_weight() const noexcept {
            return tee.total_weight();
        }

        double total_triangle_weight() const noexcept {
            return triangle.total_weight();
        }
    private:
        TransferContainer dipole;
        TransferContainer tee;
        TransferContainer triangle;
    };

    using res4_transfer_result_multi_array = res4_transfer_result<res4_transfer_multi_array_container, res4_multi_array_container>;
    using res4_transfer_result_vector = res4_transfer_result<res4_transfer_vector_container, res4_multi_array_container>;
};

#endif
