#ifndef EECs_EEC_OO_H
#define EECs_EEC_OO_H

#include <vector>
#include <memory>

#include <armadillo>

#include "toyjets/common.h"

#include "simon_util_cpp/util.h"
#include "simon_util_cpp/deltaR.h"
#include "simon_util_cpp/vecND.h"
#include "simon_util_cpp/combinatorics.h"

#include "jetinfo.h"
#include "compositions.h"
#include "maxDR.h"
#include "adj.h"

enum EECKind{
    PROJECTED=0,
    RESOLVED=1,
};

//can be either nodiagvec or symvec
//for projected or resolved respectively
template <enum EECKind K=PROJECTED, bool nonIRC=false, bool doPU=false>
class EECCalculator{
public:
    EECCalculator() :
        maxOrder_(0), J1_(), PU_(), J2_(), ptrans_(), adj_(), doTrans_(false),
        comps_() {}

    void setup(const jet& j1, const unsigned maxOrder){
        checkPU<false>();
        checkNonIRC<false>();
        maxOrder_ = maxOrder;
        J1_ = jetinfo(j1);
        comps_ = getCompositions(maxOrder_);
        initialize();
    }

    void setup(const jet& j1, const unsigned maxOrder,
               const std::vector<bool>& PU){
        checkPU<true>();
        checkNonIRC<false>();
        maxOrder_ = maxOrder;
        J1_ = jetinfo(j1);
        PU_ = PU;
        comps_ = getCompositions(maxOrder_);
        initialize();
    }

    void setup(const jet& j1, const unsigned maxOrder,
               const std::vector<bool>& PU,
               const unsigned p1, const unsigned p2){
        checkPU<true>();
        checkNonIRC<true>();
        maxOrder_ = maxOrder;
        J1_ = jetinfo(j1);
        PU_ = PU;
        comps_ = getCustomComps(p1, p2);
        p1_ = p1;
        p2_ = p2;
        initialize();
    }

    void setup(const jet& j1, const unsigned maxOrder,
               const unsigned p1, const unsigned p2){
        checkPU<false>();
        checkNonIRC<true>();
        maxOrder_ = maxOrder;
        J1_ = jetinfo(j1);
        comps_ = getCustomComps(p1, p2);
        p1_ = p1;
        p2_ = p2;
        initialize();
    }

    void setup(const jet& j1, const unsigned maxOrder,
               const arma::mat& ptrans, const jet& j2){
        ptrans_ = arma::mat(ptrans);
        adj_ = adjacency(ptrans);
        J2_ = jetinfo(j2);
        doTrans_ = true;
        setup(j1, maxOrder);
    }

    void setup(const jet& j1, const unsigned maxOrder,
               const arma::mat& ptrans, const jet& j2,
               const unsigned p1, const unsigned p2){
        ptrans_ = arma::mat(ptrans);
        adj_ = adjacency(ptrans);
        J2_ = jetinfo(j2);
        doTrans_ = true;
        setup(j1, maxOrder, p1, p2);
    }

    void run(){
        computePointAtZero();
        for(unsigned M=2; M<=maxOrder_; ++M){
            if(M > J1_.nPart){
                continue;
            } 
            computeMwayContribution(M);
        }

        finalizeCovariance();

        ran_=true;
    }

    const std::vector<double>& getResolvedDRs(unsigned order,
                                              unsigned n) const {
        checkRan();
        checkResolved();

        return resolveddRs_[order-2][n];
    }

    std::vector<double> getdRs() const{
        checkRan();

        std::vector<double> result(J1_.dR2s.data());
        result.emplace_back(0.0);
        for(unsigned i=0; i<result.size(); ++i){
            result.at(i) = std::sqrt(result.at(i));
        }
        return result;
    }

    const std::vector<double>& getwts(unsigned order) const{
        checkRan();
        checkOrder(order);
        
        return wts_[order-2];
    }

    const std::vector<double>& getwts_noPU(unsigned order) const{
        checkRan();
        checkOrder(order);
        if constexpr(!doPU){
            throw std::logic_error("Asking for wts_noPU when PU was not run");
        }
        
        return wts_noPU_[order-2];
    }

    const arma::mat& getCov(unsigned order) const{
        checkRan();
        checkOrder(order);

        return cov_[order-2];
    }

    std::vector<double> getdRs_J2() const{
        checkRan();
        checkTransfer();

        std::vector<double> result(J2_.dR2s.data());
        result.emplace_back(0.0);
        for(unsigned i=0; i<result.size(); ++i){
            result.at(i) = std::sqrt(result.at(i));
        }
        return result;
    }

    const arma::mat& getTransfer(unsigned order) const {
        checkRan();
        checkTransfer();

        if constexpr(nonIRC){
            throw std::logic_error("need to pass PU calculator for transfer normalization");
        }

        return transfer_[order-2];
    }

    const arma::mat getTransfer(unsigned order, 
                         const EECCalculator<K, true, true>& reco) const{
        if constexpr(!nonIRC){
            throw std::logic_error("passing PU calculator only necessary for nonIRC correlators");
        }

        arma::vec recoweights = reco.getwts_noPU(order);
        arma::vec transsum = arma::sum(transfer_[order-2], 1);
        transsum.replace(0, 1);

        arma::mat result =transfer_[order-2].each_col()
                             % (recoweights / transsum);
        return result;
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
    void initialize() {
        ran_ = false;

        if constexpr(K==EECKind::RESOLVED){
            size_t acc=0;
            offsets_.emplace_back(0);
            for(unsigned i=2; i<=maxOrder_; ++i){
                acc += choose(J1_.nPart, i);
                offsets_.emplace_back(acc);
            }

            if(doTrans_){
                size_t acc_J2=0;
                offsets_J2_.emplace_back(0);
                for(unsigned i=2; i<=maxOrder_; ++i){
                    acc_J2 += choose(J2_.nPart, i);
                    offsets_J2_.emplace_back(acc_J2);
                }
            }
        }

        for(unsigned order=2; order<=maxOrder_; ++order){
            size_t nconfig=0;
            if constexpr(K==EECKind::PROJECTED){
                nconfig = choose(J1_.nPart, 2);
            } else if constexpr(K==EECKind::RESOLVED){
                nconfig = offsets_[order-1];
            }
            nconfig += 1; //point at zero
            wts_.emplace_back(nconfig, 0.0);
            if constexpr(K==EECKind::RESOLVED){
                std::vector<std::vector<double>> next(choose(order,2));
                for(unsigned i=0; i<choose(order,2); ++i){
                    next.at(i).resize(nconfig, 0.0);
                }
                resolveddRs_.push_back(std::move(next));
            }
            if constexpr(doPU){
                wts_noPU_.emplace_back(nconfig);
            }
            cov_.emplace_back(nconfig, J1_.nPart, arma::fill::zeros);
            if(doTrans_){
                size_t nconfig_J2=0;
                if constexpr(K==EECKind::PROJECTED){
                    nconfig_J2 = choose(J2_.nPart, 2);
                } else if constexpr(K==EECKind::RESOLVED){
                    for(unsigned i=2; i<=order; ++i){
                        nconfig_J2 += choose(J2_.nPart, i);
                    }
                }
                nconfig_J2+=1; //point at zero
                transfer_.emplace_back(nconfig_J2,
                                      nconfig,
                                      arma::fill::zeros);
            }
        }
    }

    void checkResolved() const{
        static_assert(K == EECKind::RESOLVED, "Asking for resolved quantities on projected EEC");
    }

    template <bool havePU>
    void checkPU(){
        static_assert(doPU == havePU, "Can only call setup with PU for PU calculators");
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

    void finalizeCovariance(){
        double normfact;
        for(unsigned order=2; order<=maxOrder_; ++order){
            for(unsigned iPart=0; iPart<J1_.nPart; ++iPart){
                if constexpr(nonIRC){
                    normfact = intPow(1/(1-J1_.Es.at(iPart)), 
                            comps_.at(order-2)[0][0].composition[0]);
                } else {
                    normfact = intPow(1/(1-J1_.Es.at(iPart)), order);
                }
                for(size_t iDR=0; iDR<wts_[order-2].size(); ++iDR){
                    double contrib = -normfact*cov_[order-2](iDR,iPart);
                    double actual = (normfact-1)* wts_[order-2].at(iDR);
                    cov_[order-2](iDR, iPart) = contrib + actual;
                }
            }
        }
    }

    void computeMwayContribution(unsigned M){
        std::vector<unsigned> ord = ord0_nodiag(M);
        size_t ordidx = 0;
        do{
            accumulateWt(M, ord, ordidx);
            ++ordidx;
        } while (iterate_nodiag(M, ord, J1_.nPart));
    }

    void accumulateWt(const unsigned M,
                      std::vector<unsigned>& ord,
                      const size_t ordidx){ //ord actually is const I promise
        size_t dRidx;
        if constexpr(K==EECKind::PROJECTED){
            dRidx = getMaxDR(J1_, ord, false);
        } else if constexpr(K==EECKind::RESOLVED){
            dRidx = ordidx + offsets_[M-2];
            if(M==2){
                for(unsigned order=2; order<=maxOrder_; ++order){
                    resolveddRs_.at(order-2).at(0).at(dRidx) = std::sqrt(J1_.dR2s.at(ord));
                }
            } else if(M==3){
                std::vector<double> dRs(3);
                dRs.at(0) = std::sqrt(J1_.dR2s.at(std::vector<unsigned>({ord.at(0), ord.at(1)})));
                dRs.at(1) = std::sqrt(J1_.dR2s.at(std::vector<unsigned>({ord.at(0), ord.at(2)})));
                dRs.at(2) = std::sqrt(J1_.dR2s.at(std::vector<unsigned>({ord.at(1), ord.at(2)})));
                std::sort(dRs.begin(), dRs.end());
                for(unsigned order=3; order<=maxOrder_; ++order){
                    resolveddRs_.at(order-2).at(0).at(dRidx) = dRs.at(2);
                    resolveddRs_.at(order-2).at(1).at(dRidx) = dRs.at(1);
                    resolveddRs_.at(order-2).at(2).at(dRidx) = dRs.at(0);
                }
            } else if(M==4){
                std::vector<double> dRs(6);
                dRs.at(0) = std::sqrt(J1_.dR2s.at(std::vector<unsigned>({ord.at(0), ord.at(1)})));
                dRs.at(1) = std::sqrt(J1_.dR2s.at(std::vector<unsigned>({ord.at(0), ord.at(2)})));
                dRs.at(2) = std::sqrt(J1_.dR2s.at(std::vector<unsigned>({ord.at(0), ord.at(3)})));
                dRs.at(3) = std::sqrt(J1_.dR2s.at(std::vector<unsigned>({ord.at(1), ord.at(2)})));
                dRs.at(4) = std::sqrt(J1_.dR2s.at(std::vector<unsigned>({ord.at(1), ord.at(3)})));
                dRs.at(5) = std::sqrt(J1_.dR2s.at(std::vector<unsigned>({ord.at(2), ord.at(3)})));
                std::sort(dRs.begin(), dRs.end());
                for(unsigned order=4; order<=maxOrder_; ++order){
                    resolveddRs_.at(order-2).at(0).at(dRidx) = dRs.at(5);
                    resolveddRs_.at(order-2).at(1).at(dRidx) = dRs.at(4);
                    resolveddRs_.at(order-2).at(2).at(dRidx) = dRs.at(3);
                    resolveddRs_.at(order-2).at(3).at(dRidx) = dRs.at(2);
                    resolveddRs_.at(order-2).at(4).at(dRidx) = dRs.at(1);
                    resolveddRs_.at(order-2).at(5).at(dRidx) = dRs.at(0);
                }
            }
        }

        for(unsigned order=M; order<=maxOrder_; ++order){//for each order
            const comp_t& thiscomps = comps_.at(order-2);
            for(const comp& c : thiscomps[M-1]){//for each composition
                double nextWt = c.factor;
                bool hasPU=false;
                for(unsigned i=0; i<M; ++i){
                    nextWt *= intPow(J1_.Es.at(ord.at(i)), c.composition.at(i));
                    if constexpr(doPU){
                        if(PU_[ord.at(i)]){
                            hasPU = true;
                        }
                    }
                }
                //accumulate weight
                wts_[order-2].at(dRidx) += nextWt;
                if constexpr(doPU){
                    if(!hasPU){
                        wts_noPU_[order-2].at(dRidx) += nextWt;
                    }
                } else {
                    if(hasPU){
                        printf("this will never happen!\n");
                    }
                }
                
                //accumulate covariance
                for(unsigned i=0; i<M; ++i){
                    cov_[order-2](dRidx, ord.at(i)) += nextWt;
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
                }
            }
        }
    }

    void computePointAtZero(){
        for(unsigned order=2; order<=maxOrder_; ++order){
            size_t iDR = wts_[order-2].size()-1;
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
                wts_[order-2].at(iDR) += nextwt;
                if constexpr(doPU){
                    if(!PU_.at(i)){
                        wts_noPU_[order-2].at(iDR) += nextwt;
                    }
                }
                
                //accumulate covariance
                cov_[order-2](iDR, i) += nextwt; 

                //accumulate transfer matrix
                if(doTrans_){
                    std::vector<unsigned> ord(order, i);
                    if constexpr(nonIRC){
                        for(unsigned M=2; M<=order; ++M){
                            for(const comp& c : comps_.at(order-2)[M-1]){
                                accumulateTransfer(ord, iDR, 
                                                   rawwt, 
                                                   order, c.composition);

                            }
                        }
                    } else {
                        accumulateTransfer(ord, iDR, nextwt, order);
                    }
                }
            }
        }
    }

    void accumulateTransfer(const std::vector<unsigned>& ord_J1,
                            const size_t dRidx_J1,
                            const double nextWt,
                            const unsigned order,
                            const std::vector<unsigned>& comp = {}){

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
                    tfact *= intPow(ptrans_.at(ord_J2.at(i), ord_J1.at(i)), comp.at(i));
                } else {
                    tfact *= ptrans_.at(ord_J2.at(i), ord_J1.at(i));
                }
            }
            size_t dRidx_J2;
            if constexpr(K==EECKind::PROJECTED){
                dRidx_J2 = getMaxDR(J2_, ord_J2, true);
            } else if constexpr(K==EECKind::RESOLVED){
                std::sort(ord_J2.begin(), ord_J2.end());
                auto end = std::unique(ord_J2.begin(), ord_J2.end());
                unsigned M_J2 = std::distance(ord_J2.begin(), end);
                if(M_J2==1){
                    dRidx_J2 = offsets_J2_[order-1];
                } else {
                    dRidx_J2 = getNodiagIdx(ord_J2, J2_.nPart, M_J2) 
                             + offsets_J2_[M_J2-2];
                }
            }
            transfer_[order-2].at(dRidx_J2, dRidx_J1) += nextWt * tfact;
        } while (iterate_awkward(nadj, ord_iter));
    }

    unsigned maxOrder_;

    //the jet to do
    struct jetinfo J1_;

    //optionally track the PU-free component
    std::vector<bool> PU_;

    //only present if also doing unfolding
    struct jetinfo J2_; 
    arma::mat ptrans_;
    adjacency adj_;
    bool doTrans_;

    //precomputed quantities
    std::vector<comp_t> comps_; //indexed by [order, M] -> vector<composition>
    unsigned p1_, p2_;

    std::vector<size_t> offsets_;    //offsets into wts array for M-way
    std::vector<size_t> offsets_J2_; //contributions to resolved EEC


    //quantities to compute
    std::vector<std::vector<double>> wts_;      
    std::vector<arma::mat> cov_;
    std::vector<arma::mat> transfer_;
                       
    std::vector<std::vector<double>> wts_noPU_; 

    std::vector<std::vector<std::vector<double>>> resolveddRs_;

    bool ran_;
};

typedef EECCalculator<EECKind::PROJECTED> ProjectedEECCalculator;
typedef EECCalculator<EECKind::RESOLVED> ResolvedEECCalculator;

template <bool PU>
using NonIRCEECCalculator = EECCalculator<EECKind::PROJECTED, true, PU>;

#endif
