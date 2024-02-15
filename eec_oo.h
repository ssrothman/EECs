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

#include "compositions.h"
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

    /*
     * enum for how to normalize the jet pTs
     * RAWPT: normalize to raw pT = sum of all constituents in the jet
     *            including those that fail any particle selections
     *            (e.g. charged-only, fromPV, etc)
     * CORRPT: normalize to JEC-corrected pT
     * SUMPT: normalize to sum pT(constituents)
     *            including only those that pass particle selections
     */
    enum normType {
        RAWPT, 
        CORRPT,
        SUMPT, 
    };

    EECCalculator(int verbose=0);

    //set up simple projected calculation
    void setupProjected(const jet& j1, const unsigned maxOrder,
                        const axisptr& RLaxis, 
                        const normType& norm);
    //setup third-order resolved calculation
    void enableRes3(const axisptr& xi3axis, 
                    const axisptr& phi3axis);
    //setup fourth-order resolved calculation
    void enableRes4(const axisptr& RM4axis,
                    const axisptr& phi4axis,
                    struct trianglespec trispec);
    //add parallel calculation of unmatched contribution
    void addPU(const std::vector<bool>& PU);
    //add parallel calculation of transfer matrix
    void addTransfer(const jet& j2, 
                     const arma::mat& rawmat,
                     const normType& norm);

    //setup all internal variables and prepare for calculation
    void initialize();

    //run the calculation
    void run();

    //get results
    const std::vector<double> getproj(unsigned order) const;
    const std::vector<double> getproj_PU(unsigned order) const;
    const std::vector<double> getres3() const;
    const std::vector<double> getres4() const;
    const std::vector<double> getres3_PU() const;
    const std::vector<double> getres4_PU() const;
    const arma::mat getTransferproj(unsigned order) const;
    const arma::mat getTransferres3() const;
    const arma::mat getTransferres4() const;

    //very simple getters
    bool hasRun() const;
    unsigned getMaxOrder() const;

private:
    //renormalize particle transfer matrix appropriately 
    arma::mat makePtrans(const jet& genjet, const jet& recojet,
                         const arma::mat& rawmat, 
                         const normType& norm);
    //normalize jet pTs
    std::vector<double> normalizePt(const jet& j, const normType& norm);
    double getNormFact(const jet& j, const normType& norm);

    //throw errors if you try to access results that don't exist
    void checkPU() const;
    void checkRes3() const;
    void checkRes4() const;
    void checkOrder(unsigned order) const;
    void checkRan() const;
    void checkTransfer() const;

    //actually run the EEC calculation
    void computeMwayContribution(unsigned M);
    void accumulateWt(const unsigned M,
                      const uvec& ord);
    void computePointAtZero();

    //parameters
    unsigned maxOrder_;
    bool doRes3_;
    bool doRes4_;
    bool doTrans_;
    bool doPU_;
    int verbose_;

    //binning axes
    axisptr RLaxis_, xi3axis_, phi3axis_, RM4axis_, phi4axis_;
    struct trianglespec trispec_; //triangle for fourth order resolved

    //the jet to do
    jet J1_;
    std::vector<double> J1E_;

    //optionally track the PU component
    std::vector<bool> PU_;

    //reco jet for unfolding onto
    jet J2_;
    std::vector<double> J2E_;

    //particle-level transfer matrix
    //and adjacency list
    arma::mat ptrans_;
    adjacency adj_;

    //precomputed quantities
    std::vector<comp_t> comps_; //indexed by [order, M] -> vector<composition>

    //quantities to compute
    wtaccptr projwts_;
    wtaccptr res3wts_;
    wtaccptr res4wts_;

    wtaccptr projwts_PU_;
    wtaccptr res3wts_PU_;
    wtaccptr res4wts_PU_;

    transaccptr projtrans_;
    transaccptr res3trans_;
    transaccptr res4trans_;

    //did we actually run?
    bool ran_;
};

#endif
