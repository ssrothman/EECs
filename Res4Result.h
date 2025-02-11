#ifndef SROTHMAN_EECS_RES4_RESULT_H
#define SROTHMAN_EECS_RES4_RESULT_H

#include "usings.h"
#include "Res4Axes.h"
#include "Res4Calculator.h"

#include "SRothman/SimonTools/src/histutil.h"

#include <vector>
#include <numeric>
#include <algorithm>
#include <boost/multi_array.hpp>

namespace EEC{

    class Res4VectorContainer{
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

        Res4VectorContainer(
                const size_t nR, 
                const size_t nr, 
                const size_t nc)  noexcept :
            nR(nR), nr(nr), nc(nc){

            data.reserve(100);
        }

        Res4VectorContainer() noexcept :
            Res4VectorContainer(0, 0, 0) {}

        void fill(unsigned iR, unsigned ir, 
                  unsigned ic, double wt) noexcept{
            data.emplace_back(iR, ir, ic, wt);
        }

        const data_t& get_data() const noexcept{
            return data;
        }

        Res4VectorContainer& operator+=(const Res4VectorContainer& other) noexcept{
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


    class Res4MultiArrayContainer{
    public:
        using data_t = multi_array<double, 3>;

        Res4MultiArrayContainer() noexcept :
            Res4MultiArrayContainer(0, 0, 0) {}

        Res4MultiArrayContainer(
                const size_t nR, 
                const size_t nr, 
                const size_t nc) noexcept : 
            nR(nR), nr(nr), nc(nc){

            data.resize(boost::extents[nR][nr][nc]);
            std::fill(data.data(), 
                      data.data() + data.num_elements(), 0);
        }

        Res4MultiArrayContainer(
                const Res4VectorContainer& other) noexcept:
            Res4MultiArrayContainer(other.nR, 
                                    other.nr, 
                                    other.nc){

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

        Res4MultiArrayContainer& operator+=(const Res4MultiArrayContainer& other) noexcept{
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

        bool operator==(const Res4MultiArrayContainer& other) const noexcept{
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

    inline bool operator==(const Res4VectorContainer& a, const Res4VectorContainer& b) noexcept{
        return Res4MultiArrayContainer(a) == Res4MultiArrayContainer(b);
    }

    inline bool operator==(const Res4VectorContainer& a, const Res4MultiArrayContainer& b) noexcept{
        return Res4MultiArrayContainer(a) == b;
    }
    
    inline bool operator==(const Res4MultiArrayContainer& a, const Res4VectorContainer& b) noexcept{
        return a == Res4MultiArrayContainer(b);
    }

    inline void print_nonzero(const Res4MultiArrayContainer& a){
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

    inline void print_nonzero(const Res4VectorContainer& a){
        print_nonzero(Res4MultiArrayContainer(a));
    }

    template <class Container>
    class Res4Result{
    public:
        Res4Result() : 
            dipole(),
            tee(),
            triangle(),
            pt_denom_set(false),
            pt_denom(-1) {}

        Res4Result(
                const size_t nR, 
                const size_t nr_dipole,
                const size_t nc_dipole,
                const size_t nr_tee, 
                const size_t nc_tee,
                const size_t nr_triangle, 
                const size_t nc_triangle) noexcept :
            dipole(nR, nr_dipole, nc_dipole),
            tee(nR, nr_tee, nc_tee),
            triangle(nR, nr_triangle, nc_triangle),
            pt_denom_set(false),
            pt_denom(-1) {}

        Res4Result(const Res4Axes& axes) noexcept :
                Res4Result(simon::AXextent(axes.R),
                           simon::AXextent(axes.r_dipole),
                           simon::AXextent(axes.c_dipole),
                           simon::AXextent(axes.r_tee),
                           simon::AXextent(axes.c_tee),
                           simon::AXextent(axes.r_triangle),
                           simon::AXextent(axes.c_triangle)) {}
        
        Res4Result(const Res4Calculator& calculator) noexcept :
            Res4Result(calculator.get_axes()) {}

        Res4Result(const Res4TransferCalculator& calculator) noexcept :
            Res4Result(calculator.get_axes_gen()) {}

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
        bool operator==(const Res4Result<OtherContainer>& other) const noexcept{
            return dipole == other.get_dipole() &&
                   tee == other.get_tee() &&
                   triangle == other.get_triangle();
        }

    private:
        Container dipole;
        Container tee;
        Container triangle;
        
        bool pt_denom_set;
        double pt_denom;
    };

    using Res4Result_MultiArray = Res4Result<Res4MultiArrayContainer>;
    using Res4Result_Vector = Res4Result<Res4VectorContainer>;
};

#endif
