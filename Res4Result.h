#ifndef SROTHMAN_EECS_RES4_RESULT_H
#define SROTHMAN_EECS_RES4_RESULT_H

#include "usings.h"
#include "Res4Axes.h"
#include "ResResultContainers.h"

#include "SRothman/SimonTools/src/histutil.h"

#include <vector>
#include <numeric>
#include <boost/multi_array.hpp>

namespace EEC{
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
        
        template <typename CALCULATOR>
        Res4Result(const CALCULATOR& calculator) noexcept :
            Res4Result(calculator.get_axes()) {}

        template <typename OtherContainer>
        Res4Result(const OtherContainer& dipole, 
                   const OtherContainer& tee, 
                   const OtherContainer& triangle) noexcept :
            dipole(dipole),
            tee(tee),
            triangle(triangle),
            pt_denom_set(false),
            pt_denom(-1) {}

        template <typename OtherContainer>
        Res4Result(OtherContainer&& dipole, 
                   OtherContainer&& tee, 
                   OtherContainer&& triangle) noexcept :
            dipole(std::move(dipole)),
            tee(std::move(tee)),
            triangle(std::move(triangle)),
            pt_denom_set(false),
            pt_denom(-1) {}

        template <typename OtherContainer>
        Res4Result(const Res4Result<OtherContainer>& other) noexcept :
                dipole(other.get_dipole()),
                tee(other.get_tee()),
                triangle(other.get_triangle()){
            
            if (other.is_pt_denom_set()){
                pt_denom = other.get_pt_denom();
                pt_denom_set = true;
            } else {
                pt_denom_set = false;
                pt_denom = -1;
            }
        }

        template <typename OtherContainer>
        Res4Result(Res4Result<OtherContainer>&& other) noexcept :
                dipole(std::move(other.get_dipole())),
                tee(std::move(other.get_tee())),
                triangle(std::move(other.get_triangle())){

            if (other.is_pt_denom_set()){
                pt_denom = other.get_pt_denom();
                pt_denom_set = true;
            } else {
                pt_denom_set = false;
                pt_denom = -1;
            }
        }

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

        bool is_pt_denom_set() const noexcept{
            return pt_denom_set;
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

        template <class OtherContainer>
        bool operator!=(const Res4Result<OtherContainer>& other) const noexcept{
            return !(*this == other);
        }

        template <class OtherContainer>
        Res4Result<Container>& operator+=(const Res4Result<OtherContainer>& other) noexcept{
            dipole += other.get_dipole();
            tee += other.get_tee();
            triangle += other.get_triangle();
            return *this;
        }

        template <class OtherContainer>
        Res4Result<Container> operator+(const Res4Result<OtherContainer>& other) const noexcept{
            Res4Result<Container> result(*this);
            result += other;
            return result;
        }

        template <class OtherContainer>
        Res4Result<Container>& operator-=(const Res4Result<OtherContainer>& other) noexcept{
            if constexpr(std::is_same_v<Container, ResMultiArrayContainer>){
                dipole -= other.get_dipole();
                tee -= other.get_tee();
                triangle -= other.get_triangle();
            } else {
                //we only want a valid operator 
                //for ResMultiArrayContainer
                //this is probably the wrong way 
                //to do this,
                //but it should fail to compile
                //if you try to use this operator
                //with a ResVectorContainer
                static_assert(std::is_same_v<Container, ResMultiArrayContainer>);
            }
            return *this;
        }

        template <class OtherContainer>
        Res4Result<ResMultiArrayContainer> operator-(const Res4Result<OtherContainer>& other) const noexcept{
            Res4Result<ResMultiArrayContainer> result(*this);
            result -= other;
            return result;
        }

    private:
        Container dipole;
        Container tee;
        Container triangle;
        
        bool pt_denom_set;
        double pt_denom;
    };

    using Res4Result_MultiArray = Res4Result<ResMultiArrayContainer>;
    using Res4Result_Vector = Res4Result<ResVectorContainer<unsigned>>;
    using Res4Result_Unbinned = Res4Result<ResVectorContainer<double>>;
};

#endif
