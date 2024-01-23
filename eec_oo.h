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


enum EECKind{
    PROJECTED=0,
    RESOLVED=1,
};

template <enum EECKind K=PROJECTED, bool nonIRC=false>
class EECCalculator{
public:
    EECCalculator() :
        maxOrder_(0), J1_(), PU_(), J2_(),
        ptrans_(), adj_(), doTrans_(false),
        comps_(), ran_(false), verbose_(0), doPU_(false) {}
    
    EECCalculator(int verbose) : 
        maxOrder_(0), J1_(), PU_(), J2_(), 
        ptrans_(), adj_(), doTrans_(false),
        comps_(), ran_(false), verbose_(verbose), doPU_(false) {
            if(verbose_){
                printf("made empty EECcalculator\n");
            }
        }

    template<typename Axis>
    void setup(const jet& j1, const unsigned maxOrder, 
               const Axis& ax, bool normToRaw){
        if(verbose_){
            printf("making plain EEC calculator\n");
        }
        checkNonIRC<false>();
        maxOrder_ = maxOrder;
        J1_ = jetinfo(j1, ax, normToRaw);
        comps_ = getCompositions(maxOrder_);
        nDRbins_ = ax.size()+2;
        binAtZero_ = ax.index(0.0) + 1;
        initialize();
    }

    template<typename Axis>
    void setup(const jet& j1, const unsigned maxOrder,
               const std::vector<bool>& PU,
               const Axis& ax, bool normToRaw){
        if(verbose_){
            printf("making plain EEC calculator with PU\n");
        }
        checkNonIRC<false>();
        maxOrder_ = maxOrder;
        J1_ = jetinfo(j1, ax, normToRaw);
        comps_ = getCompositions(maxOrder_);
        PU_ = PU;
        nDRbins_ = ax.size()+2;
        binAtZero_ = ax.index(0.0) + 1;
        initialize();
    }

    template <typename Axis>
    void setup(const jet& j1, const unsigned maxOrder,
               const std::vector<bool>& PU,
               const unsigned p1, const unsigned p2,
               const Axis& ax, bool normToRaw){
        if(verbose_){
            printf("making nonIRC calculator with PU\n");
        }
        checkNonIRC<true>();
        maxOrder_ = maxOrder;
        J1_ = jetinfo(j1, ax, normToRaw);
        PU_ = PU;
        comps_ = getCustomComps(p1, p2);
        p1_ = p1;
        p2_ = p2;
        nDRbins_ = ax.size()+2;
        binAtZero_ = ax.index(0.0) + 1;
        initialize();
    }

    template <typename Axis>
    void setup(const jet& j1, const unsigned maxOrder,
               const unsigned p1, const unsigned p2,
               const Axis& ax, bool normToRaw){
        if(verbose_){
            printf("making nonIRC calculator\n");
        }
        checkNonIRC<true>();
        maxOrder_ = maxOrder;
        J1_ = jetinfo(j1, ax, normToRaw);
        comps_ = getCustomComps(p1, p2);
        p1_ = p1;
        p2_ = p2;
        nDRbins_ = ax.size()+2;
        binAtZero_ = ax.index(0.0) + 1;
        initialize();
    }

    template <typename Axis>
    void setup(const jet& j1, const unsigned maxOrder,
               const arma::mat& rawmat, const jet& j2,
               const Axis& ax, bool normToRaw){
        if(verbose_){
            printf("making calculator with transfer\n");
        }
        ptrans_ = makePtrans(j1, j2, rawmat, normToRaw);
        adj_ = adjacency(ptrans_);
        J2_ = jetinfo(j2, ax, normToRaw);
        doTrans_ = true;
        nDRbins_ = ax.size()+2;
        binAtZero_ = ax.index(0.0) + 1;
        setup(j1, maxOrder, ax, normToRaw);
    }

    template <typename Axis>
    void setup(const jet& j1, const unsigned maxOrder,
               const std::vector<bool>& PU,
               const arma::mat& rawmat, const jet& j2,
               const Axis& ax, bool normToRaw){
        if(verbose_){
            printf("making calculator with transfer\n");
        }
        ptrans_ = makePtrans(j1, j2, rawmat, normToRaw);
        adj_ = adjacency(ptrans_);
        J2_ = jetinfo(j2, ax, normToRaw);
        doTrans_ = true;
        nDRbins_ = ax.size()+2;
        binAtZero_ = ax.index(0.0) + 1;
        setup(j1, maxOrder, PU, ax, normToRaw);
    }

    template <typename Axis>
    void setup(const jet& j1, const unsigned maxOrder,
               const arma::mat& rawmat, const jet& j2,
               const unsigned p1, const unsigned p2,
               const Axis& ax, bool normToRaw){
        if(verbose_){
            printf("making nonIRC calculator with transfer\n");
        }
        ptrans_ = makePtrans(j1, j2, rawmat, normToRaw);
        adj_ = adjacency(ptrans_);
        J2_ = jetinfo(j2, ax, normToRaw);
        doTrans_ = true;
        nDRbins_ = ax.size()+2;
        binAtZero_ = ax.index(0.0) + 1;
        setup(j1, maxOrder, p1, p2, ax, normToRaw);
    }

    template <typename Axis>
    void setup(const jet& j1, const unsigned maxOrder,
               const std::vector<bool>& PU,
               const arma::mat& rawmat, const jet& j2,
               const unsigned p1, const unsigned p2,
               const Axis& ax, bool normToRaw){
        if(verbose_){
            printf("making nonIRC calculator with transfer\n");
        }
        ptrans_ = makePtrans(j1, j2, rawmat, normToRaw);
        adj_ = adjacency(ptrans_);
        J2_ = jetinfo(j2, ax, normToRaw);
        doTrans_ = true;
        nDRbins_ = ax.size()+2;
        binAtZero_ = ax.index(0.0) + 1;
        setup(j1, maxOrder, PU, p1, p2, ax, normToRaw);
    }

    void setVerbosity(int verbose){
        verbose_ = verbose;
    }

    void run(){
        if(verbose_>1){
            printf("top of run\n");
        }
        computePointAtZero();
        if(verbose_>1){
            printf("computed point at zero\n");
        }
        for(unsigned M=2; M<=maxOrder_; ++M){
            if(M > J1_.nPart){
                continue;
            } 
            computeMwayContribution(M);
            if(verbose_>1){
                printf("computed %d-way contribution\n", M);
            }
        }

        ran_=true;
    }

    const std::vector<double> getwts(unsigned order) const{
        checkRan();
        checkOrder(order);
        
        return wts_.data(order);
    }

    const std::vector<double> getwts_PU(unsigned order) const{
        checkRan();
        checkOrder(order);
        
        return wts_PU_.data(order);
    }

    const arma::mat getCov(unsigned order) const{
        checkRan();
        checkOrder(order);

        arma::mat cov(cov_.nDR(order), cov_.nPart(order));
        double factor = std::sqrt((J1_.nPart-2.0)/(2.0*J1_.nPart));
        double normfact;
        for(unsigned order=2; order<=maxOrder_; ++order){
            for(unsigned iPart=0; iPart<J1_.nPart; ++iPart){
                double base = 1/(1-J1_.Es[iPart]);
                int power; 
                if constexpr(nonIRC){
                    power = comps_[order-2][0][0].composition[0];
                } else {
                    power = order;
                }
                normfact = intPow(base, power);
                for(size_t iDR=0; iDR<wts_.nDR(order); ++iDR){
                    double contrib = -normfact * cov_.get(order, iDR, iPart);
                    double actual = (normfact-1) * wts_.get(order, iDR);
                    cov.at(iDR, iPart) = factor * (contrib + actual);
                }
            }
        }
        return cov;
    }

    const arma::mat getTransfer(unsigned order) const {
        checkRan();
        checkTransfer();

        if constexpr(nonIRC){
            throw std::logic_error("need to pass PU calculator for transfer normalization");
        }

        return transfer_.data(order);
    }

    const arma::mat getTransfer(unsigned order, 
                 const EECCalculator<K, true>& reco) const{
        if constexpr(!nonIRC){
            throw std::logic_error("passing PU calculator only necessary for nonIRC correlators");
        }

        arma::vec recoweights(reco.getwts_PU(order));
        arma::mat trans = transfer_.data(order);
        arma::vec transsum = arma::sum(trans, 1);
        transsum.replace(0, 1);

        trans = trans.each_col() % (recoweights / transsum);

        return trans;
    }

    bool hasRun() const{
        return ran_;
    }

    unsigned getMaxOrder() const {
        return maxOrder_;
    }

    unsigned getP1() const {
        return p1_;
    }

    unsigned getP2() const {
        return p2_;
    }

private:
    arma::mat makePtrans(const jet& genjet, const jet& recojet,
                         const arma::mat& rawmat, bool normToRaw) {
        arma::mat ans(rawmat);

        arma::vec genpt = genjet.ptvec();
        arma::vec recopt = recojet.ptvec();
        if(normToRaw){
            genpt/=genjet.sumpt;
            recopt/=recojet.sumpt;
        } else {
            genpt/=genjet.pt;
            recopt/=recojet.pt;
        }

        arma::vec predpt = ans * genpt;

        for(unsigned iGen=0; iGen < genjet.nPart; ++iGen){
            for(unsigned iReco=0; iReco < recojet.nPart; ++iReco){
                if (predpt[iReco] > 0){
                    ans.at(iReco, iGen) *= recopt[iReco]/predpt[iReco];
                }
            }
        }

        if (verbose_ > 1){
            printf("ptrans:\n");
            std::cout << ans << std::endl;
            printf("GEN\n");
            std::cout << genpt.t() << std::endl;
            printf("RECO\n");
            std::cout << recopt.t() << std::endl;
            printf("PRED\n");
            std::cout << (ans * genpt).t() << std::endl;
        }

        return ans;
    }

    void initialize() {
        if(verbose_){
            printf("top of initialize\n");
            printf("nDRbins = %lu\n", nDRbins_);
        }
        ran_ = false;

        if(PU_.size()){
            doPU_ = true;
        } else {
            doPU_ = false;
        }

        std::vector<size_t> nDRs;
        std::vector<size_t> nParts;
        for(unsigned order=2; order<=maxOrder_; ++order){
            if constexpr(K==EECKind::PROJECTED){
                nDRs.emplace_back(nDRbins_);
            } else {
                //there are (order choose 2) axes
                printf("for order %d, there are %lu axes\n", order, choose(order, 2));
                printf("for a total of %lu DR bins\n", intPow(nDRbins_, choose(order, 2)));
                nDRs.emplace_back(intPow(nDRbins_, choose(order, 2)));
            }
            nParts.emplace_back(J1_.nPart);
        }


        unsigned nOrder = maxOrder_-1;
        wts_ = accu1d(nOrder, nDRs);
        cov_ = accu2d(nOrder, nDRs, nParts);
        if(doPU_){
            wts_PU_ = accu1d(nOrder, nDRs);
        }

        if(doTrans_){
            transfer_ = accu2d(nOrder, nDRs, nDRs);
        }

        if(verbose_){
            printf("end of initialize\n");
        }
    }

    void checkResolved() const{
        static_assert(K == EECKind::RESOLVED, "Asking for resolved quantities on projected EEC");
    }

    template <bool haveNonIRC>
    void checkNonIRC(){
        static_assert(haveNonIRC == nonIRC, "Can only call setup with nonIRC for nonIRC calculators");
    }

    void checkOrder(unsigned order) const{
        if(order >maxOrder_){
            throw std::logic_error("Error! Asking for results of higher order than computed");
        } 
        if(order <2){
            throw std::logic_error("Observable not defined or computed for order < 2");
        }
    }

    void checkRan() const{
        if(!ran_){
            throw std::logic_error("Error! Asking for EEC results before running EEC calculation");
        }
    }

    void checkTransfer() const {
        if(!doTrans_){
            throw std::logic_error("Error! Asking for EEC transfer results when no second jet was supplied!");
        }
    }

    void computeMwayContribution(unsigned M){
        std::vector<unsigned> ord = ord0_nodiag(M);
        size_t ordidx = 0;
        do{
            if(verbose_>2){
                printf("accumulating weight for ");
                printOrd(ord);
                printf("\n");
                fflush(stdout);
            }
            accumulateWt(M, ord, ordidx);
            ++ordidx;
        } while (iterate_nodiag(M, ord, J1_.nPart));
    }

    void accumulateWt(const unsigned M,
                      std::vector<unsigned>& ord,
                      const size_t ordidx){ //ord actually is const I promise
        size_t dRidx;
        if constexpr(K==EECKind::PROJECTED){
            dRidx = getMaxDR(J1_, ord, false, binAtZero_);
            if(verbose_>2){
                printf("got dRidx = %lu\n", dRidx);
                fflush(stdout);
            }
        }

        for(unsigned order=M; order<=maxOrder_; ++order){//for each order
            const comp_t& thiscomps = comps_.at(order-2);
            for(const comp& c : thiscomps[M-1]){//for each composition
                if(verbose_>2){
                    printf("doing composition ");
                    printOrd(c.composition);
                    printf(" with factor %u\n", c.factor);
                    fflush(stdout);
                }

                if constexpr(K==EECKind::RESOLVED){
                    dRidx = getResolvedDR(J1_, ord, c.composition, binAtZero_);
                    if(verbose_>2){
                        printf("got dRidx = %lu\n", dRidx);
                        fflush(stdout);
                    }
                }

                double nextWt = c.factor;
                bool hasPU=false;
                for(unsigned i=0; i<M; ++i){
                    nextWt *= intPow(J1_.Es.at(ord.at(i)), c.composition.at(i));
                    if (doPU_ && PU_[ord.at(i)]){
                        hasPU = true;
                    }
                }
                //accumulate weight
                wts_.accumulate(order, dRidx, nextWt);
                if(verbose_>2){
                    printf("got wt = %f\n", nextWt);
                    fflush(stdout);
                }
                if (hasPU){
                    wts_PU_.accumulate(order, dRidx, nextWt);
                }
                
                //accumulate covariance
                for(unsigned i=0; i<M; ++i){
                    cov_.accumulate(order, dRidx, ord.at(i), nextWt);
                }
                if(verbose_>2){
                    printf("accumulated covariance\n");
                    fflush(stdout);
                }

                //accumulate transfer matrix
                if(doTrans_){
                    std::vector<unsigned> ord_J1;
                    if constexpr(nonIRC){
                        ord_J1 = ord;
                    } else {
                        ord_J1.reserve(order);
                        for(unsigned i=0; i<M; ++i){
                            ord_J1.insert(ord_J1.end(), c.composition.at(i), ord.at(i));
                        }
                    }

                    accumulateTransfer(ord_J1, dRidx, nextWt, order, 
                                       c.composition);
                    if(verbose_>2){
                        printf("accumulated transfer\n");
                        fflush(stdout);
                    }
                }
            }
        }
    }

    void computePointAtZero(){
        for(unsigned order=2; order<=maxOrder_; ++order){
            for(unsigned i=0; i<J1_.nPart; ++i){
                double nextwt;
                double rawwt;
                if constexpr(nonIRC){
                    rawwt = intPow(J1_.Es.at(i),
                        comps_.at(order-2)[0][0].composition[0]);
                    nextwt = rawwt * comps_.at(order-2)[0][0].factor;
                } else {
                    nextwt = intPow(J1_.Es.at(i), order);
                }

                //accumulate weight
                wts_.accumulate(order, binAtZero_, nextwt);
                if (doPU_ && PU_.at(i)){
                    wts_PU_.accumulate(order, 0, nextwt);
                }
                
                //accumulate covariance
                cov_.accumulate(order, binAtZero_, i, nextwt); 

                //accumulate transfer matrix
                if(doTrans_){
                    std::vector<unsigned> ord(order, i);
                    if constexpr(nonIRC){
                        for(unsigned M=2; M<=order; ++M){
                            for(const comp& c : comps_.at(order-2)[M-1]){
                                accumulateTransfer(ord, binAtZero_, 
                                                   rawwt, 
                                                   order, c.composition);

                            }
                        }
                    } else {
                        accumulateTransfer(ord, binAtZero_, nextwt, order);
                    }
                }
            }
        }
    }

    void accumulateTransfer(const std::vector<unsigned>& ord_J1,
                            const size_t dRidx_J1,
                            const double nextWt,
                            const unsigned order,
                            const std::vector<unsigned>& nircWts = {}){

        std::vector<unsigned> nadj(order);
        for(unsigned i=0; i<order; ++i){
            nadj.at(i) = adj_.data.at(ord_J1.at(i)).size();
            if(nadj.at(i)==0){
                return; //break out early if there are no neighbors
            }
        }

        std::vector<unsigned> ord_iter = ord0_full(order);
        do{
            double tfact = 1.0;
            std::vector<unsigned> ord_J2(order);
            for(unsigned i=0; i<order; ++i){
                ord_J2.at(i) = adj_.data.at(ord_J1.at(i))[ord_iter.at(i)];
                if constexpr(nonIRC){
                    tfact *= intPow(ptrans_.at(ord_J2.at(i), ord_J1.at(i)), nircWts.at(i));
                } else {
                    tfact *= ptrans_.at(ord_J2.at(i), ord_J1.at(i));
                }
            }
            size_t dRidx_J2;
            if constexpr(K==EECKind::PROJECTED){
                dRidx_J2 = getMaxDR(J2_, ord_J2, true, binAtZero_);
            } else if constexpr(K==EECKind::RESOLVED){
                dRidx_J2 = getResolvedDR(J2_, ord_J2, 
                                        std::vector<unsigned>(order, 1), 
                                        binAtZero_);
            }
            if(verbose_>1 && nextWt*tfact > 0.0){
                printf("transfer[%u] %lu -> %lu += %0.2g\n",
                        order, dRidx_J2, dRidx_J1, 
                        nextWt * tfact);
            }
            transfer_.accumulate(order, dRidx_J2, dRidx_J1, nextWt * tfact);
        } while (iterate_awkward(nadj, ord_iter));
    }

    unsigned maxOrder_;

    size_t nDRbins_;
    unsigned binAtZero_;

    //the jet to do
    struct jetinfo J1_;

    //optionally track the PU component
    std::vector<bool> PU_;

    //only present if also doing unfolding
    struct jetinfo J2_; 
    arma::mat ptrans_;
    adjacency adj_;
    bool doTrans_;

    //precomputed quantities
    std::vector<comp_t> comps_; //indexed by [order, M] -> vector<composition>
    unsigned p1_, p2_;

    //quantities to compute
    accu1d wts_;      
    accu2d cov_;

    //optional
    accu2d transfer_;
    accu1d wts_PU_;

    bool ran_;

    int verbose_;

    bool doPU_;
};

typedef EECCalculator<EECKind::PROJECTED, false> ProjectedEECCalculator;
typedef EECCalculator<EECKind::RESOLVED, false> ResolvedEECCalculator;

typedef EECCalculator<EECKind::PROJECTED, true> NonIRCEECCalculator;


#endif
