#ifndef SROTHMAN_EEC_V2_RESRESULT_CONTAINERS_H
#define SROTHMAN_EEC_V2_RESRESULT_CONTAINERS_H

#include "usings.h"

#include <vector>
#include <numeric>

#include <boost/multi_array.hpp>

namespace EEC{
    class ResVectorContainer{
    public:
        struct entry{
            unsigned iR;
            unsigned ir;
            unsigned ic;
            double wt;
            entry(unsigned iR, unsigned ir, 
                  unsigned ic, double wt) noexcept :
                iR(iR), ir(ir), ic(ic), wt(wt) {}
        };
        using data_t = std::vector<entry>;

        ResVectorContainer(
                const size_t nR, 
                const size_t nr, 
                const size_t nc)  noexcept :
            nR(nR), nr(nr), nc(nc){

            data.reserve(100);
        }

        ResVectorContainer() noexcept :
            ResVectorContainer(0, 0, 0) {}

        void fill(unsigned iR, unsigned ir, 
                  unsigned ic, double wt) noexcept{
            data.emplace_back(iR, ir, ic, wt);
        }

        const data_t& get_data() const noexcept{
            return data;
        }

        ResVectorContainer& operator+=(const ResVectorContainer& other) noexcept{
            data.insert(data.end(), 
                        other.data.begin(), 
                        other.data.end());
            return *this;
        }

        double total_weight() const noexcept{
            return std::accumulate(
                    data.begin(), data.end(), 0.0,
                    [](double acc, const entry& e){
                        return acc + e.wt;
                    });
        }

        const size_t nR, nr, nc;
    private:
        data_t data;
    };


    class ResMultiArrayContainer{
    public:
        using data_t = multi_array<double, 3>;

        ResMultiArrayContainer() noexcept :
            ResMultiArrayContainer(0, 0, 0) {}

        ResMultiArrayContainer(
                const size_t nR, 
                const size_t nr, 
                const size_t nc) noexcept : 
            nR(nR), nr(nr), nc(nc){

            data.resize(boost::extents[nR][nr][nc]);
            std::fill(data.data(), 
                      data.data() + data.num_elements(), 0);
        }

        ResMultiArrayContainer(
                const ResVectorContainer& other) noexcept:
            ResMultiArrayContainer(other.nR, 
                                    other.nr, 
                                    other.nc){

            for (const auto& e : other.get_data()){
                fill(e.iR, e.ir, e.ic, e.wt);
            }
        }

        void fill(unsigned iR, unsigned ir, 
                  unsigned ic, double wt) noexcept {
            //printf("ResMultiArrayContainer::fill(%u, %u, %u, %g)\n", iR, ir, ic, wt);
            //fflush(stdout);
            data[iR][ir][ic] += wt;
        }

        const data_t& get_data() const noexcept{
            return data;
        }

        ResMultiArrayContainer& operator+=(const ResMultiArrayContainer& other) noexcept{
            assert(data.num_elements() == other.data.num_elements());
            std::transform(data.data(), 
                           data.data() + data.num_elements(),
                           other.data.data(), data.data(),
                           std::plus<double>());
            return *this;
        }

        double total_weight() const noexcept{
            return std::accumulate(
                    data.data(), 
                    data.data() + data.num_elements(), 0.0);
        }

        bool operator==(const ResMultiArrayContainer& other) const noexcept{
            return std::equal(data.data(), 
                              data.data() + data.num_elements(),
                              other.data.data(),
                              [](double a, double b){
                                  return std::abs(a - b) < 1e-6;
                              });
        }

        const size_t nR, nr, nc;
    private:
        data_t data;
    };

    inline bool operator==(const ResVectorContainer& a, const ResVectorContainer& b) noexcept{
        return ResMultiArrayContainer(a) == ResMultiArrayContainer(b);
    }

    inline bool operator==(const ResVectorContainer& a, const ResMultiArrayContainer& b) noexcept{
        return ResMultiArrayContainer(a) == b;
    }
    
    inline bool operator==(const ResMultiArrayContainer& a, const ResVectorContainer& b) noexcept{
        return a == ResMultiArrayContainer(b);
    }

    inline void print_nonzero(const ResMultiArrayContainer& a){
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

    inline void print_nonzero(const ResVectorContainer& a){
        print_nonzero(ResMultiArrayContainer(a));
    }
};

#endif
