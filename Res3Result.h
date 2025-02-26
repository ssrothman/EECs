#ifndef SROTHMAN_EECS_RES4_RESULT_H
#define SROTHMAN_EECS_RES4_RESULT_H

#include "usings.h"
#include "Res3Axes.h"
#include "ResResultContainers.h"

#include "SRothman/SimonTools/src/histutil.h"

#include <vector>
#include <numeric>
#include <algorithm>
#include <boost/multi_array.hpp>

namespace EEC{
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

    using Res3Result_MultiArray = Res3Result<ResMultiArrayContainer>;
    using Res3Result_Vector = Res3Result<ResVectorContainer>;
};

#endif
