#ifndef SROTHMAN_EECS_RESOLVEDCOORDS_H
#define SROTHMAN_EECS_RESOLVEDCOORDS_H

#include "SRothman/SimonTools/src/jets.h"
#include "jetinfo.h"
#include "adj.h"

#include <vector>
#include <memory>
#include <string>

#include <boost/histogram.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/sum_kahan.hpp>

#include "FWCore/ParameterSet/interface/ParameterSet.h"

using uvec = std::vector<unsigned>;
using stats_tag = boost::accumulators::tag::sum;
using axis_t = boost::histogram::axis::variable<double>;
using axisvec = std::vector<axis_t>;
using accu_t = boost::accumulators::accumulator_set<
    double, boost::accumulators::stats<stats_tag>>;

class EECindexer{
    public:
        virtual uvec getIndex(const uvec& ord,
                              const uvec& comp) const = 0;

        uvec getIndex(const uvec& ord) const{
            uvec ones(ord.size(), 1);
            return getIndex(ord, ones);
        }

        EECindexer() = default;
        virtual ~EECindexer() {};
};

struct trianglespec{
    double RMoRL; //ratio of RM to RL
    double RSoRL; //ratio of RS to RL
    double tol; //proportional tolerance

    trianglespec() = default;

    trianglespec(double RMoRL, double RSoRL, double tol) :
        RMoRL(RMoRL), RSoRL(RSoRL), tol(tol) {}

    trianglespec(const edm::ParameterSet& conf){
        RMoRL = conf.getParameter<double>("RMoRL");
        RSoRL = conf.getParameter<double>("RSoRL");
        tol = conf.getParameter<double>("tol");
    }
};

class EECtransferAccumulator{
    public:
        EECtransferAccumulator() = default;
        ~EECtransferAccumulator() = default;

        EECtransferAccumulator(unsigned Naccu, const axisvec& axes);
        EECtransferAccumulator(unsigned Naccu, const axis_t& axis);

        void setupMaxDRIndexers(const jet& jReco, const jet& jGen);
        void setupSortedDRIndexers(const jet& jReco, const jet& jGen);
        void setupSortedDRIndexers(const jet& jReco, const jet& jGen, 
                const unsigned Nax);
        void setupThirdOrderIndexers(const jet& jReco, const jet& jGen);
        void setupFourthOrderIndexers(const jet& jReco, const jet& jGen,
                const struct trianglespec& spec);
        void setupIndexers(const edm::ParameterSet& conf,
                const jet& jReco, const jet& jGen);

        void accumulate(const uvec& ordReco, const uvec& ordGen,
                const uvec& compReco, const uvec& compGen,
                const std::vector<double>& weights);
        void accumulate(const uvec& ordGen, const uvec& compGen,
                const adjacency& adj, const arma::mat& ptrans,
                double nextwt, unsigned iAcc);

        double get(unsigned iAcc, const uvec& ord) const;
        ArbitraryMatrix<double> data(unsigned iAcc) const;
    protected:
        std::unique_ptr<EECindexer> indexerReco_;
        std::unique_ptr<EECindexer> indexerGen_;

        axisvec axes_;

        std::vector<unsigned> sizes_;
        std::vector<ArbitraryMatrix<accu_t>> accus_;

        void fillSizes(const axisvec& axes);
};

class EECweightAccumulator{
    public:
        EECweightAccumulator() = default;
        ~EECweightAccumulator() = default;

        EECweightAccumulator(unsigned Naccu, const axisvec& axes);
        EECweightAccumulator(unsigned Naccu, const axis_t& axis);

        void setupMaxDRIndexer(const jet& j);
        void setupSortedDRIndexer(const jet& j);
        void setupSortedDRIndexer(const jet& j, const unsigned Nax);
        void setupThirdOrderIndexer(const jet& j);
        void setupFourthOrderIndexer(const jet& j,
                const struct trianglespec& spec);
        void setupIndexer(const edm::ParameterSet& conf,
                const jet& j);      

        void accumulate(const uvec& ord, const uvec& comp, 
                const std::vector<double>& weights,
                unsigned offset=0,
                const uvec * const input_idx=nullptr,
                uvec* return_idx=nullptr);
        void accumulate(const uvec& ord, const uvec& comp,
                double wt, unsigned iAcc);
        void accumulate(const uvec& idx, 
                const std::vector<double>& wts,
                unsigned offset=0);

        double get(unsigned iAcc, const uvec& ord) const;
        ArbitraryMatrix<double> data(unsigned iAcc) const;
    protected:
        std::unique_ptr<EECindexer> indexer_;

        axisvec axes_;

        std::vector<unsigned> sizes_;
        std::vector<ArbitraryMatrix<accu_t>> accus_;

        void fillSizes(const axisvec& axes);

};

#endif
