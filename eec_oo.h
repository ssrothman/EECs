#ifndef EECs_EEC_OO_H
#define EECs_EEC_OO_H

#include <vector>
#include <memory>

#include <armadillo>

#include "SRothman/SimonTools/src/jets.h"

#include "SRothman/SimonTools/src/util.h"
#include "SRothman/SimonTools/src/deltaR.h"
#include "SRothman/SimonTools/src/vecND.h"
#include "SRothman/SimonTools/src/combinatorics.h"
#include <boost/histogram.hpp>
#include <boost/accumulators/accumulators.hpp>

#include "eecaccu.h"

#include "jetinfo.h"
#include "compositions.h"
#include "maxDR.h"
#include "adj.h"
#include "EECaccumulation.h"


class EECCalculator{
public:
    using uvec = std::vector<unsigned>;
    using axis_t = boost::histogram::axis::variable<double>;
    using axisptr = std::shared_ptr<axis_t>;
    using axisvec = std::vector<axis_t>;
    using wtaccptr = std::shared_ptr<EECweightAccumulator>;
    using transaccptr = std::shared_ptr<EECtransferAccumulator>;

    EECCalculator(int verbose=0);

    void setupProjected(const jet& j1, const unsigned maxOrder,
                        const axisptr& RLaxis);
    void addPU(const std::vector<bool>& PU);
    void addTransfer(const jet& j1, 
                     const jet& j2, 
                     const arma::mat& rawmat);
    void enableRes3(const axisptr& xi3axis, 
                    const axisptr& phi3axis);
    void enableRes4(const axisptr& RM4axis,
                    const axisptr& phi4axis);

    void run();

    const std::vector<double> getproj(unsigned order) const;
    const std::vector<double> getproj_PU(unsigned order) const;
    const std::vector<double> getres3() const;
    const std::vector<double> getres4() const;
    const std::vector<double> getres3_PU() const;
    const std::vector<double> getres4_PU() const;
    const arma::mat getTransferproj(unsigned order) const;
    const arma::mat getTransferres3() const;
    const arma::mat getTransferres4() const;

    bool hasRun() const;
    unsigned getMaxOrder() const;

    void initialize();
private:
    arma::mat makePtrans(const jet& genjet, const jet& recojet,
                         const arma::mat& rawmat, bool normToRaw);

    void checkPU() const;
    void checkRes3() const;
    void checkRes4() const;
    void checkOrder(unsigned order) const;
    void checkRan() const;
    void checkTransfer() const;

    void computeMwayContribution(unsigned M);
    void accumulateWt(const unsigned M,
                      const uvec& ord);
    void computePointAtZero();

    unsigned maxOrder_;
    bool doRes3_;
    bool doRes4_;

    //the jet to do
    struct jetinfo J1_;
    jet J1_J_;

    //optionally track the PU component
    std::vector<bool> PU_;

    //only present if also doing unfolding
    struct jetinfo J2_; 
    jet J2_J_;
    arma::mat ptrans_;
    adjacency adj_;
    bool doTrans_;

    //precomputed quantities
    std::vector<comp_t> comps_; //indexed by [order, M] -> vector<composition>

    //quantities to compute
    wtaccptr projwts_;
    wtaccptr res3wts_;
    wtaccptr res4wts_;

    //optional
    transaccptr projtrans_;
    transaccptr res3trans_;
    transaccptr res4trans_;

    wtaccptr projwts_PU_;
    wtaccptr res3wts_PU_;
    wtaccptr res4wts_PU_;

    bool ran_;

    int verbose_;

    bool doPU_;

    axisptr RLaxis_, xi3axis_, phi3axis_, RM4axis_, phi4axis_;
};

#endif
