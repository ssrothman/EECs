#ifndef EEC_ACCU_H
#define EEC_ACCU_H

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/sum_kahan.hpp>
#include <vector>

template <typename value_type=double, 
          typename stat_tag=boost::accumulators::tag::sum_kahan,
          bool bounds_check=false>
class EECaccu1d{
    public:
        EECaccu1d():
            accu_() {}

        EECaccu1d(unsigned nOrder, size_t nDR):
            accu_(nOrder, 
                std::vector<boost::accumulators::accumulator_set<
                value_type, boost::accumulators::stats<stat_tag>>>(nDR))
        {}

        EECaccu1d(unsigned nOrder, 
                std::vector<size_t> nDRs):
            accu_(nOrder)
        {
            for(size_t iDR = 0; iDR < nDRs.size(); ++iDR){
                accu_[iDR].resize(nDRs[iDR]);
            }
        }

        void accumulate(unsigned order, size_t iDR, value_type val) {
            return getAccu_(order, iDR)(val);
        }

        value_type get(unsigned order, size_t iDR) const{
            return boost::accumulators::extract_result<stat_tag>(getAccu_(order, iDR));
        }

        std::vector<value_type> data(unsigned order) const{
            std::vector<value_type> out(nDR(order));
            for(size_t iDR = 0; iDR < nDR(order); ++iDR){
                out[iDR] = get(order, iDR);
            }
            return out;
        }

        inline size_t nDR(unsigned order) const {
            return accu_.at(order-2).size();
        }

        inline size_t nOrder() const {
            return accu_.size();
        }

    private:
        std::vector<std::vector<boost::accumulators::accumulator_set<
            value_type, boost::accumulators::stats<
                boost::accumulators::tag::sum_kahan>>>> accu_;

        boost::accumulators::accumulator_set<value_type, 
            boost::accumulators::stats<
                boost::accumulators::tag::sum_kahan
            >
        >& getAccu_(unsigned order, size_t iDR){
            if constexpr(bounds_check){
                return accu_.at(order-2).at(iDR);
            } else {
                return accu_[order-2][iDR];
            }
        }

        const boost::accumulators::accumulator_set<value_type, 
            boost::accumulators::stats<
                boost::accumulators::tag::sum_kahan
            >
        >& getAccu_(unsigned order, size_t iDR) const{
            if constexpr(bounds_check){
                return accu_.at(order-2).at(iDR);
            } else {
                return accu_[order-2][iDR];
            }
        }
};

template <typename value_type=double, 
          typename stat_tag=boost::accumulators::tag::sum_kahan,
          bool bounds_check=false>
class EECaccu2d{
    public:
        EECaccu2d():
            accu_() {}

        EECaccu2d(unsigned nOrder, size_t nDR, size_t nPart){
            accu_.resize(nOrder);
            for(size_t iOrder = 0; iOrder < nOrder; ++iOrder){
                accu_[iOrder].resize(nDR);
                for(size_t iDR = 0; iDR < nDR; ++iDR){
                    accu_[iOrder][iDR].resize(nPart);
                }
            }
        }

        EECaccu2d(unsigned nOrder, std::vector<size_t> nDR, std::vector<size_t> nPart){
            accu_.resize(nOrder);
            for(size_t iOrder = 0; iOrder < nOrder; ++iOrder){
                accu_[iOrder].resize(nDR[iOrder]);
                for(size_t iDR = 0; iDR < nDR[iOrder]; ++iDR){
                    accu_[iOrder][iDR].resize(nPart[iOrder]);
                }
            }
        }

        void accumulate(unsigned order, size_t iDR, size_t iPart, value_type val) {
            return getAccu_(order, iDR, iPart)(val);
        }

        value_type get(unsigned order, size_t iDR, size_t iPart) const{
            return boost::accumulators::extract_result<stat_tag>(getAccu_(order, iDR, iPart));
        }

        arma::Mat<value_type> data(unsigned order) const{
            arma::Mat<value_type> out(nDR(order), nPart(order));
            for(size_t iDR = 0; iDR < nDR(order); ++iDR){
                for(size_t iPart = 0; iPart < nPart(order); ++iPart){
                    out.at(iDR, iPart) = get(order, iDR, iPart);
                }
            }
            return out;
        }

        inline size_t nDR(unsigned order) const {
            return accu_.at(order-2).size();
        }

        inline size_t nPart(unsigned order) const {
            return accu_.at(order-2).at(0).size();
        }

        inline size_t nOrder() const {
            return accu_.size();
        }

    private:
        std::vector<std::vector<std::vector<boost::accumulators::accumulator_set<
            value_type, boost::accumulators::stats<
                boost::accumulators::tag::sum_kahan
            >>>>> accu_;

        boost::accumulators::accumulator_set<value_type, 
            boost::accumulators::stats<
                boost::accumulators::tag::sum_kahan
            >
        >& getAccu_(unsigned order, size_t iDR, size_t iPart){
            if constexpr(bounds_check){
                return accu_.at(order-2).at(iDR).at(iPart);
            } else {
                return accu_[order-2][iDR][iPart];
            }
        }

        const boost::accumulators::accumulator_set<value_type, 
            boost::accumulators::stats<
                boost::accumulators::tag::sum_kahan
            >
        >& getAccu_(unsigned order, size_t iDR, size_t iPart) const{
            if constexpr(bounds_check){
                return accu_.at(order-2).at(iDR).at(iPart);
            } else {
                return accu_[order-2][iDR][iPart];
            }
        }

};

typedef EECaccu1d<double, boost::accumulators::tag::sum_kahan, true> accu1d;
typedef EECaccu2d<double, boost::accumulators::tag::sum_kahan, true> accu2d;

#endif

