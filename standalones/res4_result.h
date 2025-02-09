#ifndef SROTHMAN_EECS_RES4_RESULT_H
#define SROTHMAN_EECS_RES4_RESULT_H

#include "usings.h"
#include "boost/multi_array.hpp"

namespace standaloneEEC{
    class res4_multi_array_container{
    public:
        res4_multi_array_container(
                const size_t nR, 
                const size_t nr, 
                const size_t nc){
            data.resize(boost::extents[nR][nr][nc]);
            std::fill(data.data(), data.data() + data.num_elements(), 0);
        }

        void fill(unsigned iR, unsigned ir, 
                  unsigned ic, double wt){
            data[iR][ir][ic] += wt;
        }
    private:
        multi_array<double, 3> data;
    };

    class res4_vector_container{
    public:
        struct entry{
            unsigned iR;
            unsigned ir;
            unsigned ic;
            double wt;
            entry(unsigned iR, unsigned ir, unsigned ic, double wt) :
                iR(iR), ir(ir), ic(ic), wt(wt) {}
        };
        using data_t = std::vector<entry>;

        res4_vector_container(
                [[maybe_unused]] const size_t nR, 
                [[maybe_unused]] const size_t nr, 
                [[maybe_unused]] const size_t nc) {
            data.reserve(100);
        }

        void fill(unsigned iR, unsigned ir, 
                  unsigned ic, double wt){
            data.emplace_back(iR, ir, ic, wt);
        }
    private:
        data_t data;
    };

    template <class Container>
    class res4_result{
    public:
        res4_result(
                const size_t nR, 
                const size_t nr_dipole,
                const size_t nc_dipole,
                const size_t nr_tee, 
                const size_t nc_tee,
                const size_t nr_triangle, 
                const size_t nc_triangle) :
            dipole(nR, nr_dipole, nc_dipole),
            tee(nR, nr_tee, nc_tee),
            triangle(nR, nr_triangle, nc_triangle) {}

        void fill_dipole(
                unsigned iR, unsigned ir, 
                unsigned ic, double wt){
            dipole.fill(iR, ir, ic, wt);
        };
        void fill_tee(
                unsigned iR, unsigned ir, 
                unsigned ic, double wt){
            tee.fill(iR, ir, ic, wt);
        }
        void fill_triangle(
                unsigned iR, unsigned ir,
                unsigned ic, double wt){
            triangle.fill(iR, ir, ic, wt);
        }
    private:
        Container dipole;
        Container tee;
        Container triangle;
    };

    using res4_result_multi_array = res4_result<res4_multi_array_container>;
    using res4_result_vector = res4_result<res4_vector_container>;
};

#endif
