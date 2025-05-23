#ifndef SROTHMAN_EECS_PROJ_RESULT_H
#define SROTHMAN_EECS_PROJ_RESULT_H

#include "usings.h"
#include "ProjAxes.h"
#include "ProjResultContainers.h"

#include "SRothman/SimonTools/src/histutil.h"

#include <vector>
#include <numeric>
#include <algorithm>
#include <array>
#include <boost/multi_array.hpp>

namespace EEC{
    template <class Container>
    class ProjResult{
    public:
        ProjResult() noexcept :
            data(),
            pt_denom_set(false),
            pt_denom(-1) {}

        using T = typename Container::TYPE;
        constexpr static bool SHOULD_BIN = Container::SHOULD_BIN;
        constexpr static bool IS_ARRAY = Container::IS_ARRAY;

        ProjResult(const unsigned nR) noexcept :
            data({{nR, nR, nR, nR, nR}}),
            pt_denom_set(false),
            pt_denom(-1) {}

        ProjResult(const ProjAxes& axes) noexcept :
            ProjResult(simon::AXextent(axes.R)) {}

        template <typename CALCULATOR>
        ProjResult(const CALCULATOR& calculator) noexcept :
            ProjResult(calculator.get_axes()) {}

        ProjResult(const std::array<Container, 5>& data) noexcept :
            data(data),
            pt_denom_set(false),
            pt_denom(-1) {}

        ProjResult(std::array<Container, 5>&& data) noexcept :
            data(std::move(data)),
            pt_denom_set(false),
            pt_denom(-1) {}

        ProjResult(const Container& data0,
                   const Container& data1,
                   const Container& data2,
                   const Container& data3,
                   const Container& data4) noexcept :
            data({data0, data1, data2, data3, data4}),
            pt_denom_set(false),
            pt_denom(-1) {}

        ProjResult(Container&& data0,
                   Container&& data1,
                   Container&& data2,
                   Container&& data3,
                   Container&& data4) noexcept :
            data({std::move(data0), std::move(data1), std::move(data2),
                    std::move(data3), std::move(data4)}),
            pt_denom_set(false),
            pt_denom(-1) {}

        template <class OtherContainer>
        ProjResult(const ProjResult<OtherContainer>& other) noexcept :
            data({other.get_data()[0], other.get_data()[1],
                  other.get_data()[2], other.get_data()[3],
                  other.get_data()[4]}),
            pt_denom_set(other.is_pt_denom_set()),
            pt_denom(other.get_pt_denom_noexcept()) {}

        template <class OtherContainer>
        ProjResult(ProjResult<OtherContainer>&& other) noexcept :
            data({std::move(other.get_data()[0]),
                  std::move(other.get_data()[1]),
                  std::move(other.get_data()[2]),
                  std::move(other.get_data()[3]),
                  std::move(other.get_data()[4])}),
            pt_denom_set(other.is_pt_denom_set()),
            pt_denom(other.get_pt_denom_noexcept()) {}

        template <unsigned ORDER>
        void fill(T iR, double wt) noexcept {
            data[ORDER-2].fill(iR, wt);
        }

        void set_pt_denom(double pt_denom_) noexcept {
            pt_denom_set = true;
            pt_denom = pt_denom_;
        }

        double get_pt_denom() const {
            if (!pt_denom_set) {
                throw std::runtime_error("pt_denom not set");
            }
            return pt_denom;
        }

        double get_pt_denom_noexcept() const noexcept {
            return pt_denom;
        }

        bool is_pt_denom_set() const noexcept {
            return pt_denom_set;
        }

        const std::array<Container, 5>& get_data() const noexcept {
            return data;
        }

        double total_weight(unsigned order) const noexcept {
            return data[order-2].total_weight();
        }

        template <class OtherContainer>
        bool operator==(const ProjResult<OtherContainer>& other) const noexcept {
            return data[0] == other.get_data()[0] && 
                   data[1] == other.get_data()[1] &&
                   data[2] == other.get_data()[2] &&
                   data[3] == other.get_data()[3] &&
                   data[4] == other.get_data()[4];
        }

        template <class OtherContainer>
        bool operator!=(const ProjResult<OtherContainer>& other) const noexcept {
            return !(*this == other);
        }

        template <class OtherContainer>
        ProjResult<Container>& operator+=(const ProjResult<OtherContainer>& other) noexcept {
            data[0] += other.get_data()[0];
            data[1] += other.get_data()[1];
            data[2] += other.get_data()[2];
            data[3] += other.get_data()[3];
            data[4] += other.get_data()[4];
            return *this;
        }

        template <class OtherContainer>
        ProjResult<Container>& operator-=(const ProjResult<OtherContainer>& other) noexcept {
            data[0] -= other.get_data()[0];
            data[1] -= other.get_data()[1];
            data[2] -= other.get_data()[2];
            data[3] -= other.get_data()[3];
            data[4] -= other.get_data()[4];
            return *this;
        }

        template <class OtherContainer>
        ProjResult<ProjArrayContainer> operator+(const ProjResult<OtherContainer>& other) const noexcept {
            ProjResult<ProjArrayContainer> result(*this);
            result += other;
            return result;
        }

        template <class OtherContainer>
        ProjResult<ProjArrayContainer> operator-(const ProjResult<OtherContainer>& other) const noexcept {
            ProjResult<ProjArrayContainer> result(*this);
            result -= other;
            return result;
        }
    private:
        std::array<Container, 5> data;
        bool pt_denom_set;
        double pt_denom;
    };

    using ProjResult_Array = ProjResult<ProjArrayContainer>;
    using ProjResult_Vector = ProjResult<ProjVectorContainer<unsigned>>;
    using ProjResult_Unbinned = ProjResult<ProjVectorContainer<double>>;
};

#endif
