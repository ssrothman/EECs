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
    EECCalculator(const jet& j1,
                  const unsigned maxOrder,
                  const std::vector<bool>& PU,
                  const std::shared_ptr<std::vector<comp_t>> customComps = nullptr):
        maxOrder(maxOrder), J1(j1),
        PU(PU),
        J2(), ptrans(nullptr), adj(nullptr), 
        comps(customComps ? customComps : getCompositions(maxOrder))
    {
        checkPU(true);
        checkNonIRC(customComps);
        initialize();
    }

    EECCalculator(const jet& j1,
                  const unsigned maxOrder,
                  const std::shared_ptr<std::vector<comp_t>> customComps = nullptr):
        maxOrder(maxOrder), J1(j1),
        J2(), ptrans(nullptr), adj(nullptr), 
        comps(customComps ? customComps : getCompositions(maxOrder))
    {
        checkPU(false);
        checkNonIRC(customComps);
        initialize();
    }

    EECCalculator(const jet& j1, 
                  const unsigned maxOrder,
                  const std::shared_ptr<const arma::mat> ptrans,
                  const std::shared_ptr<const jet> j2,
                  const std::shared_ptr<std::vector<comp_t>> customComps = nullptr):
        maxOrder(maxOrder),
        J1(j1), J2(j2), ptrans(ptrans),
        adj(std::make_shared<adjacency>(ptrans)), 
        comps(customComps ? customComps : getCompositions(maxOrder))
    {
        checkPU(false);
        checkNonIRC(customComps);
        initialize();
    }

    void run(){
        computePointAtZero();
        for(unsigned M=2; M<=maxOrder; ++M){
            if(M > J1.nPart){
                continue;
            } 
            computeMwayContribution(M);
        }

        finalizeCovariance();

        ran=true;
    }

    const std::vector<double>& getResolvedDRs(unsigned order,
                                              unsigned n) const {
        checkRan();
        checkResolved();

        return resolveddRs[order-2][n];
    }

    std::vector<double> getdRs() const{
        checkRan();

        std::vector<double> result(J1.dR2s->data());
        result.emplace_back(0.0);
        for(unsigned i=0; i<result.size(); ++i){
            result.at(i) = std::sqrt(result.at(i));
        }
        return result;
    }

    const std::vector<double>& getwts(unsigned order) const{
        checkRan();
        checkOrder(order);
        
        return wts[order-2];
    }

    const std::vector<double>& getwts_noPU(unsigned order) const{
        checkRan();
        checkOrder(order);
        if constexpr(!doPU){
            throw std::logic_error("Asking for wts_noPU when PU was not run");
        }
        
        return wts_noPU[order-2];
    }

    const arma::mat& getCov(unsigned order) const{
        checkRan();
        checkOrder(order);

        return cov[order-2];
    }

    std::vector<double> getdRs_J2() const{
        checkRan();
        checkTransfer();

        std::vector<double> result(J2.dR2s->data());
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

        return transfer[order-2];
    }

    const arma::mat getTransfer(unsigned order, 
                             const EECCalculator<K, true, true>& reco){
        if constexpr(!nonIRC){
            throw std::logic_error("passing PU calculator only necessary for nonIRC correlators");
        }

        arma::vec recoweights = reco.getwts_noPU(order);
        arma::vec transsum = arma::sum(transfer[order-2], 1);
        transsum.replace(0, 1);

        arma::mat result =transfer[order-2].each_col()
                             % (recoweights / transsum);
        return result;
    }

    bool hasRun() const{
        return ran;
    }

    unsigned maxOrder;
private:
    void initialize() {
        ran = false;

        if constexpr(K==EECKind::RESOLVED){
            size_t acc=0;
            offsets.emplace_back(0);
            for(unsigned i=2; i<=maxOrder; ++i){
                acc += choose(J1.nPart, i);
                offsets.emplace_back(acc);
            }

            if(ptrans){
                size_t acc_J2=0;
                offsets_J2.emplace_back(0);
                for(unsigned i=2; i<=maxOrder; ++i){
                    acc_J2 += choose(J2.nPart, i);
                    offsets_J2.emplace_back(acc_J2);
                }
            }
        }

        for(unsigned order=2; order<=maxOrder; ++order){
            size_t nconfig=0;
            if constexpr(K==EECKind::PROJECTED){
                nconfig = choose(J1.nPart, 2);
            } else if constexpr(K==EECKind::RESOLVED){
                nconfig = offsets[order-1];
            }
            nconfig += 1; //point at zero
            wts.emplace_back(nconfig, 0.0);
            if constexpr(K==EECKind::RESOLVED){
                std::vector<std::vector<double>> next(choose(order,2));
                for(unsigned i=0; i<choose(order,2); ++i){
                    next.at(i).resize(nconfig, 0.0);
                }
                resolveddRs.push_back(std::move(next));
            }
            if constexpr(doPU){
                wts_noPU.emplace_back(nconfig);
            }
            cov.emplace_back(nconfig, J1.nPart, arma::fill::zeros);
            if(ptrans){
                size_t nconfig_J2=0;
                if constexpr(K==EECKind::PROJECTED){
                    nconfig_J2 = choose(J2.nPart, 2);
                } else if constexpr(K==EECKind::RESOLVED){
                    for(unsigned i=2; i<=order; ++i){
                        nconfig_J2 += choose(J2.nPart, i);
                    }
                }
                nconfig_J2+=1; //point at zero
                transfer.emplace_back(nconfig_J2,
                                      nconfig,
                                      arma::fill::zeros);
            }
        }
    }

    void checkResolved() const{
        if constexpr(K != EECKind::RESOLVED){
            throw std::logic_error("Asking for resolved quantities on projected EEC");
        }
    }

    void checkPU(bool wantPU){
        if (doPU != wantPU){
            throw std::logic_error("doPU template argument needs to match whether you pass PU");
        }
    }

    void checkNonIRC(std::shared_ptr<std::vector<comp_t>> customComps) const{
        if(customComps){
            if constexpr(!nonIRC){
                throw std::logic_error("Custom comps require NonIRC");
            }
            if (maxOrder != 2){
                throw std::logic_error("Custom comps only supported for second order correlators");
            }
        } else {
            if constexpr(nonIRC){
                throw std::logic_error("Non IRC requires Custom comps");
            }
        }
    }

    void checkOrder(unsigned order) const{
        if(order >maxOrder){
            throw std::logic_error("Error! Asking for results of higher order than computed");
        } 
        if(order <2){
            throw std::logic_error("Observable not defined or computed for order < 2");
        }
    }

    void checkRan() const{
        if(!ran){
            throw std::logic_error("Error! Asking for EEC results before running EEC calculation");
        }
    }

    void checkTransfer() const {
        if(!ptrans){
            throw std::logic_error("Error! Asking for EEC transfer results when no second jet was supplied!");
        }
    }

    void finalizeCovariance(){
        double normfact;
        for(unsigned order=2; order<=maxOrder; ++order){
            for(unsigned iPart=0; iPart<J1.nPart; ++iPart){
                if constexpr(nonIRC){
                    normfact = intPow(1/(1-J1.Es->at(iPart)), 
                            comps->at(order-2)[0][0].composition[0]);
                } else {
                    normfact = intPow(1/(1-J1.Es->at(iPart)), order);
                }
                for(size_t iDR=0; iDR<wts[order-2].size(); ++iDR){
                    double contrib = -normfact*cov[order-2](iDR,iPart);
                    double actual = (normfact-1)* wts[order-2].at(iDR);
                    cov[order-2](iDR, iPart) = contrib + actual;
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
        } while (iterate_nodiag(M, ord, J1.nPart));
    }

    void accumulateWt(const unsigned M,
                      std::vector<unsigned>& ord,
                      const size_t ordidx){ //ord actually is const I promise
        size_t dRidx;
        if constexpr(K==EECKind::PROJECTED){
            dRidx = getMaxDR(J1, ord, false);
        } else if constexpr(K==EECKind::RESOLVED){
            dRidx = ordidx + offsets[M-2];
            if(M==2){
                resolveddRs.at(0).at(0).at(dRidx) = std::sqrt(J1.dR2s->at(ord));
                resolveddRs.at(1).at(0).at(dRidx) = std::sqrt(J1.dR2s->at(ord));
                resolveddRs.at(2).at(0).at(dRidx) = std::sqrt(J1.dR2s->at(ord));
            } else if(M==3){
                std::vector<double> dRs(3);
                dRs.at(0) = std::sqrt(J1.dR2s->at(std::vector<unsigned>({ord.at(0), ord.at(1)})));
                dRs.at(1) = std::sqrt(J1.dR2s->at(std::vector<unsigned>({ord.at(0), ord.at(2)})));
                dRs.at(2) = std::sqrt(J1.dR2s->at(std::vector<unsigned>({ord.at(1), ord.at(2)})));
                std::sort(dRs.begin(), dRs.end());
                resolveddRs.at(1).at(0).at(dRidx) = dRs.at(2);
                resolveddRs.at(1).at(1).at(dRidx) = dRs.at(1);
                resolveddRs.at(1).at(2).at(dRidx) = dRs.at(0);
                resolveddRs.at(2).at(0).at(dRidx) = dRs.at(2);
                resolveddRs.at(2).at(1).at(dRidx) = dRs.at(1);
                resolveddRs.at(2).at(2).at(dRidx) = dRs.at(0);
            } else if(M==4){
                std::vector<double> dRs(6);
                dRs.at(0) = std::sqrt(J1.dR2s->at(std::vector<unsigned>({ord.at(0), ord.at(1)})));
                dRs.at(1) = std::sqrt(J1.dR2s->at(std::vector<unsigned>({ord.at(0), ord.at(2)})));
                dRs.at(2) = std::sqrt(J1.dR2s->at(std::vector<unsigned>({ord.at(0), ord.at(3)})));
                dRs.at(3) = std::sqrt(J1.dR2s->at(std::vector<unsigned>({ord.at(1), ord.at(2)})));
                dRs.at(4) = std::sqrt(J1.dR2s->at(std::vector<unsigned>({ord.at(1), ord.at(3)})));
                dRs.at(5) = std::sqrt(J1.dR2s->at(std::vector<unsigned>({ord.at(2), ord.at(3)})));
                std::sort(dRs.begin(), dRs.end());
                resolveddRs.at(2).at(0).at(dRidx) = dRs.at(5);
                resolveddRs.at(2).at(1).at(dRidx) = dRs.at(4);
                resolveddRs.at(2).at(2).at(dRidx) = dRs.at(3);
                resolveddRs.at(2).at(3).at(dRidx) = dRs.at(2);
                resolveddRs.at(2).at(4).at(dRidx) = dRs.at(1);
                resolveddRs.at(2).at(5).at(dRidx) = dRs.at(0);
            }
        }

        for(unsigned order=M; order<=maxOrder; ++order){//for each order
            const comp_t& thiscomps = comps->at(order-2);
            for(const comp& c : thiscomps[M-1]){//for each composition
                double nextWt = c.factor;
                bool hasPU=false;
                for(unsigned i=0; i<M; ++i){
                    nextWt *= intPow(J1.Es->at(ord.at(i)), c.composition.at(i));
                    if constexpr(doPU){
                        if(PU[ord.at(i)]){
                            hasPU = true;
                        }
                    }
                }
                //accumulate weight
                wts[order-2].at(dRidx) += nextWt;
                if constexpr(doPU){
                    if(!hasPU){
                        wts_noPU[order-2].at(dRidx) += nextWt;
                    }
                } else {
                    if(hasPU){
                        printf("this will never happen!\n");
                    }
                }
                
                //accumulate covariance
                for(unsigned i=0; i<M; ++i){
                    cov[order-2](dRidx, ord.at(i)) += nextWt;
                }

                //accumulate transfer matrix
                if(ptrans){
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
        for(unsigned order=2; order<=maxOrder; ++order){
            size_t iDR = wts[order-2].size()-1;
            for(unsigned i=0; i<J1.nPart; ++i){
                double nextwt;
                double rawwt;
                if constexpr(nonIRC){
                    rawwt = intPow(J1.Es->at(i),
                        comps->at(order-2)[0][0].composition[0]);
                    nextwt = rawwt * comps->at(order-2)[0][0].factor;
                } else {
                    nextwt = intPow(J1.Es->at(i), order);
                }

                //accumulate weight
                wts[order-2].at(iDR) += nextwt;
                if constexpr(doPU){
                    if(!PU.at(i)){
                        wts_noPU[order-2].at(iDR) += nextwt;
                    }
                }
                
                //accumulate covariance
                cov[order-2](iDR, i) += nextwt; 

                //accumulate transfer matrix
                if(ptrans){
                    std::vector<unsigned> ord(order, i);
                    if constexpr(nonIRC){
                        for(unsigned M=2; M<=order; ++M){
                            for(const comp& c : comps->at(order-2)[M-1]){
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
            nadj.at(i) = adj->data.at(ord_J1.at(i)).size();
            if(nadj.at(i)==0){
                return; //break out early if there are no neighbors
            }
        }

        std::vector<unsigned> ord_iter = ord0_full(order);
        do{
            double tfact = 1.0;
            std::vector<unsigned> ord_J2(order);
            for(unsigned i=0; i<order; ++i){
                ord_J2.at(i) = adj->data.at(ord_J1.at(i))[ord_iter.at(i)];
                if constexpr(nonIRC){
                    tfact *= intPow(ptrans->at(ord_J2.at(i), ord_J1.at(i)),
                                    comp.at(i));
                } else {
                    tfact *= ptrans->at(ord_J2.at(i), ord_J1.at(i));
                }
            }
            size_t dRidx_J2;
            if constexpr(K==EECKind::PROJECTED){
                dRidx_J2 = getMaxDR(J2, ord_J2, true);
            } else if constexpr(K==EECKind::RESOLVED){
                std::sort(ord_J2.begin(), ord_J2.end());
                auto end = std::unique(ord_J2.begin(), ord_J2.end());
                unsigned M_J2 = std::distance(ord_J2.begin(), end);
                if(M_J2==1){
                    dRidx_J2 = offsets_J2[order-1];
                } else {
                    dRidx_J2 = getNodiagIdx(ord_J2, J2.nPart, M_J2) 
                             + offsets_J2[M_J2-2];
                }
            }
            transfer[order-2].at(dRidx_J2, dRidx_J1) += nextWt * tfact;
        } while (iterate_awkward(nadj, ord_iter));
    }

    //the jet to do
    const struct jetinfo J1;

    //optionally track the PU-free component
    const std::vector<bool> PU;

    //only present if also doing unfolding
    const struct jetinfo J2; 
    const std::shared_ptr<const arma::mat> ptrans;
    const std::shared_ptr<adjacency> adj;

    //precomputed quantities
    const std::shared_ptr<std::vector<comp_t>> comps; //indexed by [order, M] -> vector<composition>

    std::vector<size_t> offsets;    //offsets into wts array for M-way
    std::vector<size_t> offsets_J2; //contributions to resolved EEC


    //quantities to compute
    std::vector<std::vector<double>> wts;      
    std::vector<arma::mat> cov;
    std::vector<arma::mat> transfer;
                       
    std::vector<std::vector<double>> wts_noPU; 

    std::vector<std::vector<std::vector<double>>> resolveddRs;

    bool ran;
};

typedef EECCalculator<EECKind::PROJECTED> ProjectedEECCalculator;
typedef EECCalculator<EECKind::RESOLVED> ResolvedEECCalculator;

template <bool PU>
using NonIRCEECCalculator = EECCalculator<EECKind::PROJECTED, true, PU>;

#endif
