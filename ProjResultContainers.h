#ifndef SROTHMAN_EEC_PROJ_CONTAINERS_H
#define SROTHMAN_EEC_PROJ_CONTAINERS_H

#include "usings.h"
#include <vector>
#include <numeric>
#include <boost/multi_array.hpp>

namespace EEC{
    template <typename T>
    class ProjVectorContainer;

    class ProjArrayContainer{
    public:
        using data_t = std::vector<double>;
        using TYPE = unsigned;

        constexpr static bool SHOULD_BIN = true;
        constexpr static bool IS_ARRAY = true;
            
        ProjArrayContainer() noexcept :
            nR(0),
            data() {}

        ProjArrayContainer(const size_t nR) noexcept :
            nR(nR),
            data(nR, 0.0) {}

        ProjArrayContainer(ProjArrayContainer&& other) noexcept :
            nR(other.nR),
            data(std::move(other.data)) {}

        ProjArrayContainer(const ProjArrayContainer& other) noexcept :
            nR(other.nR),
            data(other.data) {}

        ProjArrayContainer(const ProjVectorContainer<unsigned>& other) noexcept;

        void fill(unsigned iR, double wt) noexcept{
            data[iR] += wt;
        }

        const data_t& get_data() const noexcept{
            return data;
        }

        double total_weight() const noexcept {
            return std::accumulate(data.begin(), data.end(), 0.0);
        }

        ProjArrayContainer& operator+=(const ProjArrayContainer& other) noexcept{
            std::transform(data.begin(), data.end(),
                          other.data.begin(), data.begin(),
                          std::plus<double>());
            return *this;
        }

        ProjArrayContainer& operator+=(const ProjVectorContainer<unsigned> other) noexcept;

        template <class OtherContainer>
        ProjArrayContainer operator+(const OtherContainer& other) const noexcept{
            ProjArrayContainer result(*this);
            result += other;
            return result;
        }

        ProjArrayContainer& operator-=(const ProjArrayContainer& other) noexcept{
            std::transform(data.begin(), data.end(),
                          other.data.begin(), data.begin(),
                          std::minus<double>());
            return *this;
        }

        ProjArrayContainer& operator-=(const ProjVectorContainer<unsigned>& other) noexcept;

        template <class OtherContainer>
        ProjArrayContainer operator-(const OtherContainer& other) const noexcept{
            ProjArrayContainer result(*this);
            result -= other;
            return result;
        }

        bool operator==(const ProjArrayContainer& other) const noexcept{
            return std::equal(data.begin(), data.end(),
                               other.data.begin(),
                               [](double a, double b){
                                   return std::abs(a - b) < 1e-6;
                               });
        }

        bool operator==(const ProjVectorContainer<unsigned>& other) const noexcept {
            return *this== ProjArrayContainer(other);
        }

        const size_t nR;
    private:
        data_t data;

    };

    template <typename T>
    class ProjVectorContainer{
    public:
        struct entry{
            T iR;
            double wt;
            entry(T iR, double wt) noexcept :
                iR(iR), wt(wt) {}
        };
        using data_t = std::vector<entry>;
        using TYPE = T;

        constexpr static bool SHOULD_BIN = std::is_same<T, unsigned>::value;
        constexpr static bool IS_ARRAY = false;

        ProjVectorContainer(const size_t nR) noexcept :
            nR(nR){
            data.reserve(100);
        }

        ProjVectorContainer() noexcept :
            ProjVectorContainer(0) {}

        ProjVectorContainer(ProjVectorContainer&& other) noexcept :
            nR(other.nR),
            data(std::move(other.data)) {}

        ProjVectorContainer(const ProjVectorContainer& other) noexcept :
            nR(other.nR),
            data(other.data) {}

        void fill(T iR, double wt) noexcept{
            data.emplace_back(iR, wt);
        }

        ProjVectorContainer<T>& operator+=(const ProjVectorContainer<T>& other) noexcept {
            data.insert(data.end(),
                        other.data.begin(),
                        other.data.end());
            return *this;
        }

        template <class OtherContainer>
        ProjVectorContainer<T> operator+(const OtherContainer& other) const noexcept {
            ProjVectorContainer<T> result(*this);
            result += other;
            return result;
        }

        template <class OtherContainer>
        ProjArrayContainer operator-(const OtherContainer& other) const noexcept {
            ProjArrayContainer result(*this);
            result -= other;
            return result;
        }

        bool operator==(const ProjVectorContainer& other) const noexcept{
            return ProjArrayContainer(*this) == ProjArrayContainer(other);
        }

        bool operator==(const ProjArrayContainer& other) const noexcept{
            return ProjArrayContainer(*this) == other;
        }

        const data_t& get_data() const noexcept{
            return data;
        }

        double total_weight() const noexcept {
            return std::accumulate(data.begin(), data.end(), 0.0,
                                   [](double sum, const entry& e) {
                                       return sum + e.wt;
                                   });
        }

        const size_t nR;
    private:
        data_t data;
    };
};

#endif
