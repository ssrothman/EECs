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

    void setupProjected(const jet& j1, const unsigned maxOrder,
                        const axisptr& RLaxis){
        RLaxis_ = RLaxis;
        maxOrder_ = maxOrder;
        //TODO: normToRaw
        J1_ = jetinfo(j1, *RLaxis, true);
        J1_J_ = j1;
        comps_ = getCompositions(maxOrder_);

        doPU_ = false;
        doTrans_= false;
        doRes3_ = false;
        doRes4_ = false;
    }
    
    void addPU(const std::vector<bool>& PU){
        PU_ = PU;
        doPU_ = true;
    }

    void addTransfer(const jet& j1, 
                     const jet& j2, 
                     const arma::mat& rawmat){
        //TODO: normToRaw
        ptrans_ = makePtrans(j1, j2, rawmat, true);
        adj_ = adjacency(ptrans_);
        //TODO: normToRaw
        J2_ = jetinfo(j2, *RLaxis_, true);
        J2_J_ = j2;
        doTrans_ = true;
    }

    void enableRes3(const axisptr& xi3axis, 
                    const axisptr& phi3axis){
        doRes3_ = true;
        xi3axis_ = xi3axis;
        phi3axis_ = phi3axis;
    }

    void enableRes4(const axisptr& RM4axis,
                    const axisptr& phi4axis){
        doRes4_ = true;
        RM4axis_ = RM4axis;
        phi4axis_ = phi4axis;
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

    const std::vector<double> getproj(unsigned order) const{
        checkRan();
        checkOrder(order);
        
        return projwts_->data(order-2).data();
    }

    const std::vector<double> getproj_PU(unsigned order) const{
        checkRan();
        checkOrder(order);
        checkPU();
        
        return projwts_PU_->data(order-2).data();
    }

    const std::vector<double> getres3() const{
        checkRan();
        checkOrder(3);
        checkRes3();

        return res3wts_->data(0).data();
    }

    const std::vector<double> getres4() const{
        checkRan();
        checkOrder(4);
        checkRes4();

        return res4wts_->data(0).data();
    }

    const std::vector<double> getres3_PU() const{
        checkRan();
        checkOrder(3);
        checkRes3();
        checkPU();

        return res3wts_PU_->data(0).data();
    }

    const std::vector<double> getres4_PU() const{
        checkRan();
        checkOrder(4);
        checkRes4();
        checkPU();

        return res4wts_PU_->data(0).data();
    }

    const arma::mat getTransferproj(unsigned order) const {
        checkRan();
        checkTransfer();
        checkOrder(order);

        ArbitraryMatrix<double> AM = projtrans_->data(order-2);
        uvec dims = AM.dims();
        unsigned ndim = dims.size()/2;
        unsigned size = 1;
        for(unsigned i=0; i<ndim; ++i){
            size *= dims[i];
        }
        const std::vector<double>& data = AM.data();
        arma::mat result(size, size, arma::fill::none);
        unsigned k=0;
        for(unsigned i=0; i<size; ++i){
            for(unsigned j=0; j<size; ++j){
                result(i, j) = data[k++];
            }
        }
        return result;
        //TODO: check if this is correct
    }

    const arma::mat getTransferres3() const {
        checkRan();
        checkTransfer();
        checkOrder(3);
        checkRes3();

        ArbitraryMatrix<double> AM = res3trans_->data(0);
        uvec dims = AM.dims();
        unsigned ndim = dims.size()/2;
        unsigned size = 1;
        for(unsigned i=0; i<ndim; ++i){
            size *= dims[i];
        }
        const std::vector<double>& data = AM.data();
        arma::mat result(size, size, arma::fill::none);
        unsigned k=0;
        for(unsigned i=0; i<size; ++i){
            for(unsigned j=0; j<size; ++j){
                result(i, j) = data[k++];
            }
        }
        return result;
        //TODO: check if this is correct
    }

    const arma::mat getTransferres4() const {
        checkRan();
        checkTransfer();
        checkOrder(4);
        checkRes4();

        ArbitraryMatrix<double> AM = res4trans_->data(0);
        uvec dims = AM.dims();
        unsigned ndim = dims.size()/2;
        unsigned size = 1;
        for(unsigned i=0; i<ndim; ++i){
            size *= dims[i];
        }
        const std::vector<double>& data = AM.data();
        arma::mat result(size, size, arma::fill::none);
        unsigned k=0;
        for(unsigned i=0; i<size; ++i){
            for(unsigned j=0; j<size; ++j){
                result(i, j) = data[k++];
            }
        }
        return result;
        //TODO: check if this is correct
    }

    bool hasRun() const{
        return ran_;
    }

    unsigned getMaxOrder() const {
        return maxOrder_;
    }

    void initialize(){
        if(verbose_){
            printf("top of initialize\n");
        }
        ran_ = false;

        if(PU_.size()){
            doPU_ = true;
        } else {
            doPU_ = false;
        }

        unsigned Nacc = maxOrder_ - 1;
        projwts_ = std::make_shared<EECweightAccumulator>(Nacc, *RLaxis_);
        projwts_->setupMaxDRIndexer(J1_J_);
        
        res3wts_ = std::make_shared<EECweightAccumulator>(1, axisvec({*RLaxis_, *xi3axis_, 
                                                                        *phi3axis_}));
        res3wts_->setupThirdOrderIndexer(J1_J_);

        res4wts_ = std::make_shared<EECweightAccumulator>(1, axisvec({*RLaxis_, *RM4axis_, 
                                                                        *phi4axis_}));
        res4wts_->setupFourthOrderIndexer(J1_J_, 
                trianglespec(0,0,0));

        if(doPU_){
            projwts_PU_ = std::make_shared<EECweightAccumulator>(Nacc, *RLaxis_);
            projwts_PU_->setupMaxDRIndexer(J1_J_);

            res3wts_PU_ = std::make_shared<EECweightAccumulator>(1, axisvec({*RLaxis_, 
                                                                            *xi3axis_, 
                                                                            *phi3axis_}));
            res3wts_PU_->setupThirdOrderIndexer(J1_J_);

            res4wts_PU_ = std::make_shared<EECweightAccumulator>(1, axisvec({*RLaxis_, 
                                                                            *RM4axis_, 
                                                                        *phi4axis_}));
            res4wts_PU_->setupFourthOrderIndexer(J1_J_, 
                    trianglespec(0,0,0));
        }

        if(doTrans_){
            projtrans_ = std::make_shared<EECtransferAccumulator>(Nacc, *RLaxis_);
            projtrans_->setupMaxDRIndexers(J2_J_, J1_J_);

            res3trans_ = std::make_shared<EECtransferAccumulator>(1, axisvec({*RLaxis_,
                                                                            *xi3axis_,
                                                                            *phi3axis_}));
            res3trans_->setupThirdOrderIndexers(J2_J_, J1_J_);

            res4trans_ = std::make_shared<EECtransferAccumulator>(1, axisvec({*RLaxis_, 
                                                                            *RM4axis_,
                                                                        *phi4axis_}));
            res4trans_->setupFourthOrderIndexers(J2_J_, J1_J_, 
                    trianglespec(0,0,0));
        }

        if(verbose_){
            printf("end of initialize\n");
        }
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

    void checkPU() const{
        if(!doPU_){
            throw std::logic_error("Error! Asking for PU results without doing PU calculation");
        }
    }

    void checkRes3() const{
        if(!doRes3_){
            throw std::logic_error("Error! Asking for res3 results without doing res3 calculation");
        }
    }

    void checkRes4() const{
        if(!doRes4_){
            throw std::logic_error("Error! Asking for res4 resoluts without doing res4 calculation");
        }
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
        uvec ord = ord0_nodiag(M);
        do{
            if(verbose_>2){
                printf("accumulating weight for ");
                printOrd(ord);
                printf("\n");
                fflush(stdout);
            }
            accumulateWt(M, ord);
        } while (iterate_nodiag(M, ord, J1_.nPart));
    }

    void accumulateWt(const unsigned M,
                      const uvec& ord){

        //wait to accumulate projected weights till end of loop
        //we can do this becase the maxDR binning is indep of 
        //composition
        //
        //Note that this is not the case for the transfer fills
        //or for the resolved correlators
        std::vector<double> wts, wts_PU;
        wts.resize(maxOrder_-M+1, 0);
        if(doPU_){
            wts_PU.resize(maxOrder_-M+1, 0);
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

                double nextWt = c.factor;
                bool hasPU=false;
                for(unsigned i=0; i<M; ++i){
                    unsigned iP = ord.at(i);
                    nextWt *= intPow(J1_.Es.at(iP), c.composition.at(i));
                    if (doPU_ && PU_[ord.at(i)]){
                        hasPU = true;
                    }
                }
                wts[order-M] += nextWt;
                if (hasPU){
                    wts_PU[order-M] += nextWt;
                }

                if (doTrans_){
                    projtrans_->accumulate(ord, c.composition, 
                            adj_, ptrans_, nextWt, order-2);
                }
                if(doRes3_ && order==3){
                    res3wts_->accumulate(ord, c.composition, {nextWt});
                    if(hasPU){
                        res3wts_PU_->accumulate(ord, c.composition, 
                                {nextWt});
                    }
                    if (doTrans_){
                        res3trans_->accumulate(ord, c.composition,
                                adj_, ptrans_, nextWt, order-2);
                    }
                } else if(doRes4_ && order==4){
                    res4wts_->accumulate(ord, c.composition, {nextWt});
                    if(hasPU){
                        res4wts_PU_->accumulate(ord, c.composition,
                                {nextWt});
                    }
                    if (doTrans_){
                        res4trans_->accumulate(ord, c.composition,
                                adj_, ptrans_, nextWt, order-2);
                    }
                }
            }//end for each composition
        }//end for each order
        projwts_->accumulate(ord, uvec(M, 1), wts, M-2); 
        if(doPU_){
            projwts_PU_->accumulate(ord, uvec(M, 1), wts_PU, M-2);
        }
    }

    void computePointAtZero(){
        for(unsigned i=0; i<J1_.nPart; ++i){//for each particle
            std::vector<double> wts, wts_PU;
            wts.reserve(maxOrder_-1);
            if(doPU_){
                wts_PU.reserve(maxOrder_-1);
            }
            for(unsigned order=2; order<=maxOrder_; ++order){//for each order
                double nextwt;
                nextwt = intPow(J1_.Es.at(i), order);

                //accumulate weight
                wts.push_back(nextwt); 
                if (doPU_ && PU_.at(i)){
                    wts_PU.push_back(nextwt);
                }
            }//end for each order
            projwts_->accumulate(uvec({i}), uvec({1}), wts);
            if(doRes3_){
                res3wts_->accumulate(uvec({i}), uvec({1}), 
                                     std::vector<double>({wts[1]}));
            }
            if(doRes4_){
                res4wts_->accumulate(uvec({i}), uvec({1}), 
                                     std::vector<double>({wts[2]}));
            }
            if(doPU_){
                projwts_PU_->accumulate(uvec({i}), uvec({1}), wts_PU);
                if(doRes3_){
                    res3wts_PU_->accumulate(uvec({i}), uvec({1}),
                                            std::vector<double>({wts_PU[1]}));
                }
                if(doRes4_){
                    res4wts_PU_->accumulate(uvec({i}), uvec({1}), 
                                            std::vector<double>({wts_PU[2]}));
                }
            }
            if(doTrans_){
                for(unsigned order=2; order<=maxOrder_; ++order){
                    projtrans_->accumulate(uvec({i}), uvec({order}), adj_, ptrans_,
                                           wts[order-2], order-2);
                }
                if(doRes3_){
                    res3trans_->accumulate(uvec({i}), uvec({3}), adj_, ptrans_,
                                          wts[1], 0);
                }
                if(doRes4_){
                    res4trans_->accumulate(uvec({i}), uvec({4}), adj_, ptrans_,
                                          wts[2], 0);
                }
            }
        }//end for each particle
    }

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
    unsigned p1_, p2_;

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
