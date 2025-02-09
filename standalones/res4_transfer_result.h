#ifndef SROTHMAN_EECS_RES4_TRANSFER_RESULT_H
#define SROTHMAN_EECS_RES4_TRANSFER_RESULT_H

#include "usings.h"
#include "boost/multi_array.hpp"
#include <vector>

namespace standaloneEEC{
    class res4_transfer_multi_array_container{
    public:
        res4_transfer_multi_array_container(
                const size_t nR_1, 
                const size_t nr_1, 
                const size_t nc_1,
                const size_t nR_2, 
                const size_t nr_2, 
                const size_t nc_2){
            data.resize(boost::extents[nR_1][nr_1][nc_1][nR_2][nr_2][nc_2]);
            std::fill(data.data(), data.data() + data.num_elements(), 0);
        }

        void fill(unsigned iR1, unsigned ir1, unsigned ic1,
                  unsigned iR2, unsigned ir2, unsigned ic2,
                  double wt){
            data[iR1][ir1][ic1][iR2][ir2][ic2] += wt;
        }
    private:
        multi_array<double, 6> data;
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
                  double wt) :
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
                [[maybe_unused]] const size_t nc_2) {
            data.reserve(100);
        }

        void fill(unsigned iR1, unsigned ir1, unsigned ic1,
                  unsigned iR2, unsigned ir2, unsigned ic2,
                  double wt){
            data.emplace_back(iR1, ir1, ic1, iR2, ir2, ic2, wt);
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
                const size_t nc_triangle_2) :
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
                         double wt){
            dipole.fill(iR1, ir1, ic1, iR2, ir2, ic2, wt);
        }

        void fill_tee(unsigned iR1, unsigned ir1, unsigned ic1,
                      unsigned iR2, unsigned ir2, unsigned ic2,
                      double wt){
            tee.fill(iR1, ir1, ic1, iR2, ir2, ic2, wt);
        }

        void fill_triangle(unsigned iR1, unsigned ir1, unsigned ic1,
                           unsigned iR2, unsigned ir2, unsigned ic2,
                           double wt){
            triangle.fill(iR1, ir1, ic1, iR2, ir2, ic2, wt);
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
