#ifndef SROTHMAN_EEC_PROJ_TRANSFER_CONTAINERS_H
#define SROTHMAN_EEC_PROJ_TRANSFER_CONTAINERS_H

#include "usings.h"
#include <vector>
#include <numeric>
#include <boost/multi_array.hpp>

#include "ProjResultContainers.h"

namespace EEC{
    template <typename T>
    class ProjTransferVectorContainer;

    class ProjTransferArrayContainer{
    public:
        using data_t = multi_array<double, 2>;
        using TYPE = unsigned;

        constexpr static bool SHOULD_BIN = true;

        ProjTransferArrayContainer() noexcept :
            nR_reco(0), nR_gen(0),
            data() {
                printf("WARNING: CALLING DEFAULT CONSTRUCTOR FOR ProjTransferArrayContainer\n");
            }

        ProjTransferArrayContainer(const size_t nR_reco, 
                                   const size_t nR_gen) noexcept :
            nR_reco(nR_reco), nR_gen(nR_gen){
                data.resize(boost::extents[nR_reco][nR_gen]);
                std::fill(data.data(), data.data() + data.num_elements(), 0.0);
            }

        ProjTransferArrayContainer(ProjTransferArrayContainer&& other) noexcept :
            nR_reco(other.nR_reco), nR_gen(other.nR_gen){
        
            data.resize(boost::extents[nR_reco][nR_gen]);
            data = other.get_data();
        }

        ProjTransferArrayContainer(const ProjTransferArrayContainer& other) noexcept :
            nR_reco(other.nR_reco), nR_gen(other.nR_gen){

            data.resize(boost::extents[nR_reco][nR_gen]);
            data = other.get_data();
        }

        ProjTransferArrayContainer(const ProjTransferVectorContainer<unsigned>& other) noexcept;

        void fill(unsigned iR_reco, unsigned iR_gen, double wt_reco, [[maybe_unused]] double wt_gen) noexcept{
            data[iR_reco][iR_gen] += wt_reco;
        }

        const data_t& get_data() const noexcept{
            return data;
        }

        double total_weight_reco() const noexcept {
            return std::accumulate(data.data(), data.data()+data.num_elements(), 0.0);  
        }

        double total_weight_gen() const noexcept {
            return std::accumulate(data.data(), data.data()+data.num_elements(), 0.0);  
        }

        ProjTransferArrayContainer& operator+=(const ProjTransferArrayContainer& other) noexcept{
            std::transform(data.data(), data.data()+data.num_elements(),
                          other.data.data(), data.data(),
                          std::plus<double>());
            return *this;
        }

        ProjTransferArrayContainer& operator+=(const ProjTransferVectorContainer<unsigned> other) noexcept;

        template <class OtherContainer>
        ProjTransferArrayContainer operator+(const OtherContainer& other) const noexcept{
            ProjTransferArrayContainer result(*this);
            result += other;
            return result;
        }

        ProjTransferArrayContainer operator-=(const ProjTransferArrayContainer& other) noexcept{
            std::transform(data.data(), data.data()+data.num_elements(),
                          other.data.data(), data.data(),
                          std::minus<double>());
            return *this;
        }

        ProjTransferArrayContainer operator-=(const ProjTransferVectorContainer<unsigned>& other) noexcept;

        template <class OtherContainer>
        ProjTransferArrayContainer operator-(const OtherContainer& other) const noexcept{
            ProjTransferArrayContainer result(*this);
            result -= other;
            return result;
        }

        bool operator==(const ProjTransferArrayContainer& other) const noexcept{
            return std::equal(data.data(), data.data()+data.num_elements(),
                               other.data.data(),
                               [](double a, double b){
                                   return std::abs(a - b) < 1e-6;
                               });
        }

        bool operator==(const ProjTransferVectorContainer<unsigned>& other) const noexcept {
            return *this== ProjTransferArrayContainer(other);
        }

        ProjArrayContainer get_sum_over_gen() const noexcept{
            ProjArrayContainer sum(nR_reco);
            for (size_t iR_reco = 0; iR_reco < nR_reco; ++iR_reco){
                for (size_t iR_gen = 0; iR_gen < nR_gen; ++iR_gen){
                    sum.fill(iR_reco, data[iR_reco][iR_gen]);
                }
            }
            return sum;
        }

        ProjArrayContainer get_sum_over_reco() const noexcept{
            ProjArrayContainer sum(nR_gen);
            for (size_t iR_reco = 0; iR_reco < nR_reco; ++iR_reco){
                for (size_t iR_gen = 0; iR_gen < nR_gen; ++iR_gen){
                    sum.fill(iR_gen, data[iR_reco][iR_gen]);
                }
            }
            return sum;
        }

        const size_t nR_reco, nR_gen;
    private:
        data_t data;

    };

    template <typename T>
    class ProjTransferVectorContainer{
    public:
        struct entry{
            T iR_reco, iR_gen;
            double wt_reco, wt_gen;
            entry(T iR_reco, T iR_gen, double wt_reco, double wt_gen) noexcept :
                iR_reco(iR_reco), iR_gen(iR_gen), wt_reco(wt_reco), wt_gen(wt_gen) {}
        };
        using data_t = std::vector<entry>;
        using TYPE = T;

        constexpr static bool SHOULD_BIN = std::is_same<T, unsigned>::value;

        ProjTransferVectorContainer(const size_t nR_reco,
                                   const size_t nR_gen) noexcept :
            nR_reco(nR_reco), nR_gen(nR_gen){
            data.reserve(100);
        }

        ProjTransferVectorContainer() noexcept :
            nR_reco(0), nR_gen(0),
            data() {
                printf("WARNING: CALLING DEFAULT CONSTRUCTOR FOR ProjTransferVectorContainer\n");
            }

        ProjTransferVectorContainer(ProjTransferVectorContainer&& other) noexcept :
            nR_reco(other.nR_reco), nR_gen(other.nR_gen),
            data(std::move(other.data)) {}

        ProjTransferVectorContainer(const ProjTransferVectorContainer& other) noexcept :
            nR_reco(other.nR_reco), nR_gen(other.nR_gen),
            data(other.data) {}

        void fill(T iR_reco, T iR_gen, double wt_reco, double wt_gen) noexcept{
            data.emplace_back(iR_reco, iR_gen, wt_reco, wt_gen);
        }

        ProjTransferVectorContainer<T>& operator+=(const ProjTransferVectorContainer<T>& other) noexcept {
            data.insert(data.end(),
                        other.data.begin(),
                        other.data.end());
            return *this;
        }

        template <class OtherContainer>
        ProjTransferVectorContainer<T> operator+(const OtherContainer& other) const noexcept {
            ProjTransferVectorContainer<T> result(*this);
            result += other;
            return result;
        }

        template <class OtherContainer>
        ProjTransferArrayContainer operator-(const OtherContainer& other) const noexcept {
            ProjTransferArrayContainer result(*this);
            result -= other;
            return result;
        }

        bool operator==(const ProjTransferVectorContainer& other) const noexcept{
            return ProjTransferArrayContainer(*this) == ProjTransferArrayContainer(other);
        }

        bool operator==(const ProjTransferArrayContainer& other) const noexcept{
            return ProjTransferArrayContainer(*this) == other;
        }

        const data_t& get_data() const noexcept{
            return data;
        }

        double total_weight_reco() const noexcept {
            return std::accumulate(data.begin(), data.end(), 0.0,
                                   [](double sum, const entry& e) {
                                       return sum + e.wt_reco;
                                   });
        }

        double total_weight_gen() const noexcept {
            return std::accumulate(data.begin(), data.end(), 0.0,
                                   [](double sum, const entry& e) {
                                       return sum + e.wt_gen;
                                   });
        }

        ProjArrayContainer get_sum_over_gen() const noexcept{
            ProjArrayContainer sum(nR_reco);
            for (const auto& entry : data){
                sum.fill(entry.iR_reco, entry.wt_reco);
            }
            return sum;
        }

        ProjArrayContainer get_sum_over_reco() const noexcept{
            ProjArrayContainer sum(nR_gen);
            for (const auto& entry : data){
                sum.fill(entry.iR_gen, entry.wt_reco);
            }
            return sum;
        }

        const size_t nR_reco, nR_gen;
    private:
        data_t data;
    };
};

#endif
