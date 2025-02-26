#ifndef SROTHMAN_EECS_RES4_RESULT_H
#define SROTHMAN_EECS_RES4_RESULT_H

#include "usings.h"
#include "Res3Axes.h"

#include "SRothman/SimonTools/src/histutil.h"

#include <vector>
#include <numeric>
#include <algorithm>
#include <boost/multi_array.hpp>

namespace EEC{

    class Res3VectorContainer{
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

        Res3VectorContainer(
                const size_t nR, 
                const size_t nr, 
                const size_t nc)  noexcept :
            nR(nR), nr(nr), nc(nc){

            data.reserve(100);
        }

        Res3VectorContainer() noexcept :
            Res3VectorContainer(0, 0, 0) {}

        void fill(unsigned iR, unsigned ir, 
                  unsigned ic, double wt) noexcept{
            data.emplace_back(iR, ir, ic, wt);
        }

        const data_t& get_data() const noexcept{
            return data;
        }

        Res3VectorContainer& operator+=(const Res3VectorContainer& other) noexcept{
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


    class Res3MultiArrayContainer{
    public:
        using data_t = multi_array<double, 3>;

        Res3MultiArrayContainer() noexcept :
            Res3MultiArrayContainer(0, 0, 0) {}

        Res3MultiArrayContainer(
                const size_t nR, 
                const size_t nr, 
                const size_t nc) noexcept : 
            nR(nR), nr(nr), nc(nc){

            data.resize(boost::extents[nR][nr][nc]);
            std::fill(data.data(), 
                      data.data() + data.num_elements(), 0);
        }

        Res3MultiArrayContainer(
                const Res3VectorContainer& other) noexcept:
            Res3MultiArrayContainer(other.nR, 
                                    other.nr, 
                                    other.nc){

            for (const auto& e : other.get_data()){
                fill(e.iR, e.ir, e.ic, e.wt);
            }
        }

        void fill(unsigned iR, unsigned ir, 
                  unsigned ic, double wt) noexcept {
            //printf("Res3MultiArrayContainer::fill(%u, %u, %u, %g)\n", iR, ir, ic, wt);
            //fflush(stdout);
            data[iR][ir][ic] += wt;
        }

        const data_t& get_data() const noexcept{
            return data;
        }

        Res3MultiArrayContainer& operator+=(const Res3MultiArrayContainer& other) noexcept{
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

        bool operator==(const Res3MultiArrayContainer& other) const noexcept{
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

    inline bool operator==(const Res3VectorContainer& a, const Res3VectorContainer& b) noexcept{
        return Res3MultiArrayContainer(a) == Res3MultiArrayContainer(b);
    }

    inline bool operator==(const Res3VectorContainer& a, const Res3MultiArrayContainer& b) noexcept{
        return Res3MultiArrayContainer(a) == b;
    }
    
    inline bool operator==(const Res3MultiArrayContainer& a, const Res3VectorContainer& b) noexcept{
        return a == Res3MultiArrayContainer(b);
    }

    inline void print_nonzero(const Res3MultiArrayContainer& a){
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

    inline void print_nonzero(const Res3VectorContainer& a){
        print_nonzero(Res3MultiArrayContainer(a));
    }

    template <class Container>
    class Res3Result{
    public:
        Res3Result() : 
            data(),
            pt_denom_set(false),
            pt_denom(-1) {}

        Res3Result(
                const size_t nR, 
                const size_t nr, 
                const size_t nc) noexcept :
            data(nR, nr, nc),
            pt_denom_set(false),
            pt_denom(-1) {}

        Res3Result(const Res3Axes& axes) noexcept :
                Res3Result(simon::AXextent(axes.R),
                           simon::AXextent(axes.r),
                           simon::AXextent(axes.c)) {}
        
        template <typename CALCULATOR>
        Res3Result(const CALCULATOR& calculator) noexcept :
            Res3Result(calculator.get_axes()) {}

        void fill(
                unsigned iR, unsigned ir,
                unsigned ic, double wt) noexcept{
            data.fill(iR, ir, ic, wt);
        }

        void set_pt_denom(double pt_denom_) {
            if (pt_denom_set){
                throw std::runtime_error("pt_denom already set");
            } else {
                pt_denom = pt_denom_;
                pt_denom_set = true;
            }
        }

        double get_pt_denom() const {
            if (!pt_denom_set){
                throw std::runtime_error("pt_denom not set");
            } else {
                return pt_denom;
            }
        }

        const Container& get_data() const noexcept{
            return data;
        }

        double total_weight() const noexcept{
            return data.total_weight();
        }

        template <class OtherContainer>
        bool operator==(const Res3Result<OtherContainer>& other) const noexcept{
            return data == other.get_data();
        }

    private:
        Container data;
        
        bool pt_denom_set;
        double pt_denom;
    };

    using Res3Result_MultiArray = Res3Result<Res3MultiArrayContainer>;
    using Res3Result_Vector = Res3Result<Res3VectorContainer>;
};

#endif
