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
