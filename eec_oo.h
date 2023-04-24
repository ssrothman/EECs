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

/*
 * 3 possible invocations
 *
 * 1) EECCalculator(jet, order) will just calculate the EEC, no frills
 * 2) EECCalculator(jet, order, ptrans, jet) will calculate the EEC and a transfer matrix
 * 3) EECCalculator(jet, order, nullptr, nullptr, customComps) will calculate a nonIRC EEC
 *
 * note that EECCalculator(jet, order, ptrans, jet, customComps) is NOT supported
 */

enum EECKind{
    PROJECTED=0,
    RESOLVED=1
};

//can be either nodiagvec or symvec
//for projected or resolved respectively
template <enum EECKind K=PROJECTED>
class EECCalculator{
public:
    EECCalculator(const std::shared_ptr<const jet> j1,
                  const unsigned maxOrder,
                  const std::shared_ptr<std::vector<comp_t>> customComps = nullptr):
        maxOrder(maxOrder), J1(j1),
        J2(nullptr), ptrans(nullptr), adj(nullptr), 
        comps(customComps ? customComps : getCompositions(maxOrder))
    {
        initialize();
    }

    EECCalculator(const std::shared_ptr<const jet> j1, 
                  const unsigned maxOrder,
                  const std::shared_ptr<const arma::mat> ptrans,
                  const std::shared_ptr<const jet> j2):
        maxOrder(maxOrder),
        J1(j1), J2(j2), ptrans(ptrans),
        adj(std::make_unique<adjacency>(ptrans)), 
        comps(getCompositions(maxOrder))
    {
        initialize();
    }

    void run(){
        printf("top of run\n");
        computePointAtZero();
        printf("did point at zero\n");
        for(unsigned M=2; M<=maxOrder; ++M){
            if(M > J1.nPart){
                printf("skipping %u\n", M);
                continue;
            } 
            computeMwayContribution(M);
            printf("did %u-way\n",M);
        }

        finalizeCovariance();
        printf("finalized cov\n");

        ran=true;
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

        return transfer[order-2];
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
        printf("initialized\n");
        printf("offsets:\n");
        printOrd(offsets);
        printf("\n");
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
                normfact = intPow(1/(1-J1.Es->at(iPart)), order);
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
        }

        for(unsigned order=M; order<=maxOrder; ++order){//for each order
            const comp_t& thiscomps = comps->at(order-2);
            for(const comp& c : thiscomps[M-1]){//for each composition
                double nextWt = c.factor;
                for(unsigned i=0; i<M; ++i){
                    nextWt *= intPow(J1.Es->at(ord[i]), c.composition[i]);
                }
                //accumulate weight
                wts[order-2].at(dRidx) += nextWt;
                
                //accumulate covariance
                for(unsigned i=0; i<M; ++i){
                    cov[order-2](dRidx, ord[i]) += nextWt;
                }

                //accumulate transfer matrix
                if(ptrans){
                    std::vector<unsigned> ord_J1;
                    ord_J1.reserve(order);
                    for(unsigned i=0; i<M; ++i){
                        ord_J1.insert(ord_J1.end(), c.composition[i], ord[i]);
                    }

                    accumulateTransfer(ord_J1, dRidx, nextWt, order);
                }
            }
        }
    }

    void computePointAtZero(){
        for(unsigned order=2; order<=maxOrder; ++order){
            size_t iDR = wts[order-2].size()-1;
            for(unsigned i=0; i<J1.nPart; ++i){
                double nextwt = intPow(J1.Es->at(i), order);
                nextwt *= comps->at(order-2)[0][0].factor; //this is probably always 1...
                //accumulate weight
                wts[order-2].at(iDR) += nextwt;
                
                //accumulate covariance
                cov[order-2](iDR, i) += nextwt; 

                //accumulate transfer matrix
                if(ptrans){
                    std::vector<unsigned> ord(order, i);
                    accumulateTransfer(ord, J1.dR2s->size(), nextwt, order);
                }
            }
        }
    }

    void accumulateTransfer(const std::vector<unsigned>& ord_J1,
                            const unsigned dRidx_J1,
                            const double nextWt,
                            const unsigned order){

        std::vector<unsigned> nadj(order);
        for(unsigned i=0; i<order; ++i){
            nadj[i] = adj->data.at(ord_J1[i]).size();
            if(nadj[i]==0){
                return; //break out early if there are no neighbors
            }
        }

        std::vector<unsigned> ord_iter = ord0_full(order);
        do{
            double tfact = 1.0;
            std::vector<unsigned> ord_J2(order);
            for(unsigned i=0; i<order; ++i){
                ord_J2[i] = adj->data.at(ord_J1[i])[ord_iter[i]];
                tfact *= ptrans->at(ord_J2[i], ord_J1[i]);
            }
            size_t dRidx_J2 = getMaxDR(J2, ord_J2, true);
            transfer[order-2].at(dRidx_J2, dRidx_J1) += nextWt * tfact;
        } while (iterate_awkward(nadj, ord_iter));
    }

    //the jet to do
    const struct jetinfo J1;

    //only present if also doing unfolding
    const struct jetinfo J2; 
    const std::shared_ptr<const arma::mat> ptrans;
    const std::unique_ptr<adjacency> adj;

    //precomputed quantities
    const std::shared_ptr<std::vector<comp_t>> comps; //indexed by [order, M] -> vector<composition>
                                                      //
    std::vector<size_t> offsets; //offsets into wts array for M-way
                                 //contributions to resolved EEC

    //quantities to compute
    std::vector<std::vector<double>> wts; //shape of weights vec different for projected vs resolved
    std::vector<arma::mat> cov;
    std::vector<arma::mat> transfer;
                       
    bool ran;
};

typedef EECCalculator<EECKind::PROJECTED> ProjectedEECCalculator;
typedef EECCalculator<EECKind::RESOLVED> ResolvedEECCalculator;

#endif
