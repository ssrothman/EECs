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

class EECCalculator{
public:
    EECCalculator(const std::shared_ptr<const jet> j1, 
                  const unsigned order,
                  const std::shared_ptr<const arma::mat> ptrans=nullptr,
                  const std::shared_ptr<const jet> j2=nullptr,
                  const std::shared_ptr<comp_t> customComps = nullptr): 
        order(order), J1(j1), 
        J2(j2), ptrans(ptrans),
        adj(std::make_unique<adjacency>(ptrans)), 
        comps(customComps ? customComps : getCompositions(order)),
        powOrder(comps->at(0)[0].composition[0]),
        wts(J1.nPart, 2u),
        cov(wts.size()+1, J1.nPart, arma::fill::zeros),
        ptAtZero(0),
        transfer(j2 ? J2.dR2s->size()+1 : 0,
                 j2 ? J1.dR2s->size()+1 : 0, 
                 arma::fill::zeros),
        ran(false) {}

    void run(){
        computePointAtZero();
        for(unsigned M=2; M<=order; ++M){
            computeMwayContribution(M);
        }

        finalizeCovariance();

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

    std::vector<double> getwts() const{
        checkRan();
        
        std::vector<double> result(wts.data());
        result.emplace_back(ptAtZero);
        return result;
    }

    const arma::mat& getCov() const{
        checkRan();

        return cov;
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

    const arma::mat& getTransfer() const {
        checkRan();
        checkTransfer();

        return transfer;
    }

    //info about the calculation
    const unsigned order;

private:
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
        for(unsigned iPart=0; iPart<J1.nPart; ++iPart){
            normfact = intPow(1/(1-J1.Es->at(iPart)), powOrder);
            for(unsigned iDR=0; iDR<J1.dR2s->size(); ++iDR){
                double contrib = -normfact * cov(iDR, iPart);
                double actual = (normfact-1) * wts.at(iDR);
                cov(iDR, iPart) = contrib + actual;
            }
            //special case for point at zero
            unsigned iDR = J1.dR2s->size();
            cov(iDR, iPart) = -normfact * cov(iDR, iPart)
                                  +(normfact-1) * ptAtZero;
        }
    }

    void computeMwayContribution(unsigned M){
        std::vector<unsigned> ord = ord0_nodiag(M);
        bool loop;
        do{
            accumulateWt(M, ord);
        } while (iterate_nodiag(M, ord, J1.nPart));
    }

    void accumulateWt(const unsigned M,
                      const std::vector<unsigned>& ord){
        unsigned dRidx = getMaxDR(J1, ord);

        for(const comp& c : comps->at(M-1)){//for each composition
            double nextWt = c.factor;
            for(unsigned i=0; i<M; ++i){
                nextWt *= intPow(J1.Es->at(ord[i]), c.composition[i]);
            }
            //accumulate weight
            wts.at(dRidx) += nextWt;
            
            //accumulate covariance
            for(unsigned i=0; i<M; ++i){
                cov(dRidx, ord[i]) += nextWt;
            }

            //accumulate transfer matrix
            if(ptrans){
                std::vector<unsigned> ord_J1;
                ord_J1.reserve(powOrder);
                for(unsigned i=0; i<M; ++i){
                    ord_J1.insert(ord_J1.end(), c.composition[i], ord[i]);
                }

                accumulateTransfer(ord_J1, dRidx, nextWt);
            }
        }
    }

    void computePointAtZero(){
        for(unsigned i=0; i<J1.nPart; ++i){
            double nextwt = intPow(J1.Es->at(i), powOrder);
            nextwt *= comps->at(0)[0].factor;
            printf("F = %u\n", comps->at(0)[0].factor);

            //accumulate weight
            ptAtZero += nextwt;
            
            //accumulate covariance
            cov(J1.dR2s->size(), i) += nextwt; 

            //accumulate transfer matrix
            if(ptrans){
                std::vector<unsigned> ord(powOrder, i);
                accumulateTransfer(ord, J1.dR2s->size(), nextwt);
            }
        }
    }

    void accumulateTransfer(const std::vector<unsigned>& ord_J1,
                            const unsigned dRidx_J1,
                            const double nextWt){

        std::vector<unsigned> nadj(powOrder);
        for(unsigned i=0; i<powOrder; ++i){
            nadj[i] = adj->data.at(ord_J1[i]).size();
            if(nadj[i]==0){
                return; //break out early if there are no neighbors
            }
        }

        printf("Transfering ");
        printOrd(ord_J1);
        printf(":\n");
        std::vector<unsigned> ord_iter = ord0_full(powOrder);
        do{
            float tfact = 1.0f;
            std::vector<unsigned> ord_J2(powOrder);
            for(unsigned i=0; i<powOrder; ++i){
                ord_J2[i] = adj->data.at(ord_J1[i])[ord_iter[i]];
                tfact *= ptrans->at(ord_J2[i], ord_J1[i]);
            }
            printf("  -> ");
            printOrd(ord_J2);
            printf(" [which is ");
            printOrd(ord_iter);
            printf("]\n");
            printf(": (%0.5f) * (%0.5f) \n\t\t\t\t= %0.5f\n", tfact, nextWt, tfact*nextWt);
            unsigned dRidx_J2 = getMaxDR(J2, ord_J2);
            transfer.at(dRidx_J2, dRidx_J1) += nextWt * tfact;
        } while (iterate_awkward(nadj, ord_iter));
    }


    //the jet to do
    const struct jetinfo J1;

    //only present if also doing unfolding
    const struct jetinfo J2; 
    const std::shared_ptr<const arma::mat> ptrans;
    const std::unique_ptr<adjacency> adj;

    //precomputed quantities
    const std::shared_ptr<comp_t> comps;
    const unsigned powOrder; //for customComps

    //quantities to compute
    vecND::nodiagvec wts;
    arma::mat cov;
    float ptAtZero;
    arma::mat transfer; //only if doing unfolding
    
                       
    bool ran;
};

#endif

