#ifndef SROTHMAN_EECS_RES4_RESULT_H
#define SROTHMAN_EECS_RES4_RESULT_H

#include "usings.h"
#include "boost/multi_array.hpp"

namespace standaloneEEC{
    class res4_vector_container{
    public:
        struct entry{
            unsigned iR;
            unsigned ir;
            unsigned ic;
            double wt;
            entry(unsigned iR, unsigned ir, unsigned ic, double wt) noexcept :
                iR(iR), ir(ir), ic(ic), wt(wt) {}
        };
        using data_t = std::vector<entry>;

        res4_vector_container(
                const size_t nR, 
                const size_t nr, 
                const size_t nc)  noexcept :
            nR(nR), nr(nr), nc(nc){

            data.reserve(100);
        }

        void fill(unsigned iR, unsigned ir, 
                  unsigned ic, double wt) noexcept{
            data.emplace_back(iR, ir, ic, wt);
        }

        const data_t& get_data() const noexcept{
            return data;
        }

        res4_vector_container& operator+=(const res4_vector_container& other) noexcept{
            data.insert(data.end(), other.data.begin(), other.data.end());
            return *this;
        }

        double total_weight() const noexcept{
            return std::accumulate(data.begin(), data.end(), 0.0,
                                   [](double acc, const entry& e){
                                       return acc + e.wt;
                                   });
        }

        const size_t nR, nr, nc;
    private:
        data_t data;
    };


    class res4_multi_array_container{
    public:
        using data_t = multi_array<double, 3>;

        res4_multi_array_container(
                const size_t nR, 
                const size_t nr, 
                const size_t nc) noexcept : 
            nR(nR), nr(nr), nc(nc){

            data.resize(boost::extents[nR][nr][nc]);
            std::fill(data.data(), data.data() + data.num_elements(), 0);
        }

        res4_multi_array_container(
                const res4_vector_container& other) noexcept:
            res4_multi_array_container(other.nR, other.nr, other.nc){

            for (const auto& e : other.get_data()){
                fill(e.iR, e.ir, e.ic, e.wt);
            }
        }

        void fill(unsigned iR, unsigned ir, 
                  unsigned ic, double wt) noexcept {
            data[iR][ir][ic] += wt;
        }

        const data_t& get_data() const noexcept{
            return data;
        }

        res4_multi_array_container& operator+=(const res4_multi_array_container& other) noexcept{
            assert(data.num_elements() == other.data.num_elements());
            std::transform(data.data(), data.data() + data.num_elements(),
                           other.data.data(), data.data(),
                           std::plus<double>());
            return *this;
        }

        double total_weight() const noexcept{
            return std::accumulate(data.data(), data.data() + data.num_elements(), 0.0);
        }

        bool operator==(const res4_multi_array_container& other) const noexcept{
            for(size_t iR=0; iR<nR; ++iR){
                for(size_t ir=0; ir<nr; ++ir){
                    for(size_t ic=0; ic<nc; ++ic){
                        if (data[iR][ir][ic] > 0 || other.data[iR][ir][ic] > 0){
                            //printf("(%lu, %lu, %lu): this = %g, other = %g\n", iR, ir, ic, data[iR][ir][ic], other.data[iR][ir][ic]);
                        }
                    }
                }
            }
            return std::equal(data.data(), data.data() + data.num_elements(),
                              other.data.data(),
                              [](double a, double b){
                                  return std::abs(a - b) < 1e-6;
                              });
        }

        const size_t nR, nr, nc;
    private:
        data_t data;
    };

    bool operator==(const res4_vector_container& a, const res4_vector_container& b) noexcept{
        return res4_multi_array_container(a) == res4_multi_array_container(b);
    }

    bool operator==(const res4_vector_container& a, const res4_multi_array_container& b) noexcept{
        return res4_multi_array_container(a) == b;
    }
    
    bool operator==(const res4_multi_array_container& a, const res4_vector_container& b) noexcept{
        return a == res4_multi_array_container(b);
    }

    void print_nonzero(const res4_multi_array_container& a){
        for(size_t iR=0; iR<a.nR; ++iR){
            for(size_t ir=0; ir<a.nr; ++ir){
                for(size_t ic=0; ic<a.nc; ++ic){
                    if (a.get_data()[iR][ir][ic] > 0){
                        printf("(%lu, %lu, %lu): %g\n", iR, ir, ic, a.get_data()[iR][ir][ic]);
                    }
                }
            }
        }
    }

    void print_nonzero(const res4_vector_container& a){
        print_nonzero(res4_multi_array_container(a));
    }

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
                const size_t nc_triangle) noexcept :
            dipole(nR, nr_dipole, nc_dipole),
            tee(nR, nr_tee, nc_tee),
            triangle(nR, nr_triangle, nc_triangle)  {}

        void fill_dipole(
                unsigned iR, unsigned ir, 
                unsigned ic, double wt) noexcept{
            dipole.fill(iR, ir, ic, wt);
        };
        void fill_tee(
                unsigned iR, unsigned ir, 
                unsigned ic, double wt) noexcept{
            tee.fill(iR, ir, ic, wt);
        }
        void fill_triangle(
                unsigned iR, unsigned ir,
                unsigned ic, double wt) noexcept{
            triangle.fill(iR, ir, ic, wt);
        }

        const Container& get_dipole() const noexcept{
            return dipole;
        }

        const Container& get_tee() const noexcept{
            return tee;
        }

        const Container& get_triangle() const noexcept{
            return triangle;
        }

        double total_dipole_weight() const noexcept{
            return dipole.total_weight();
        }

        double total_tee_weight() const noexcept{
            return tee.total_weight();
        }

        double total_triangle_weight() const noexcept{
            return triangle.total_weight();
        }

        template <class OtherContainer>
        bool operator==(const res4_result<OtherContainer>& other) const noexcept{
            return dipole == other.get_dipole() &&
                   tee == other.get_tee() &&
                   triangle == other.get_triangle();
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
