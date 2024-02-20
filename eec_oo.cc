#include "eec_oo.h"

EECCalculator::EECCalculator(int verbose) {
    maxOrder_ = 0;
    ran_ = false;
    doPU_ = false;
    doTrans_ = false;
    doRes3_ = false;
    doRes4_ = false;
    verbose_ = verbose;
    if(verbose_>1){
        printf("made empty calculator\n");
    }

    duration_total_ = 0;
    duration_transferproj_ = 0;
    duration_res3_ = 0;
    duration_res4_ = 0;
    duration_transferres3_ = 0;
    duration_transferres4_ = 0;
}

void EECCalculator::setupProjected(const jet& j1, 
                                   const unsigned maxOrder,
                                   const axisptr& RLaxis,
                                   const normType& norm){
    RLaxis_ = RLaxis;
    maxOrder_ = maxOrder;
    
    J1_ = j1;
    J1E_ = normalizePt(j1, norm);
    comps_ = getCompositions(maxOrder_);
    if(verbose_>1){
        printf("set up projected\n");
    }
}

void EECCalculator::addPU(const std::vector<bool>& PU){
    PU_ = PU;
    doPU_ = true;
    if(verbose_>1){
        printf("added PU\n");
    }
}

void EECCalculator::addTransfer(
                 const jet& j2, 
                 const arma::mat& rawmat,
                 const normType& norm){
    ptrans_ = makePtrans(J1_, j2, rawmat, norm);
    adj_ = adjacency(ptrans_);

    J2_ = j2;
    J2E_ = normalizePt(j2, norm);

    doTrans_ = true;
    if(verbose_>1){
        printf("added transfer\n");
    }
}

void EECCalculator::enableRes3(const axisptr& xi3axis, 
                const axisptr& phi3axis){
    doRes3_ = true;
    xi3axis_ = xi3axis;
    phi3axis_ = phi3axis;
    if(verbose_>1){
        printf("enabled res3\n");
    }
}

void EECCalculator::enableRes4(const axisptr& RM4axis,
                               const axisptr& phi4axis,
                               struct trianglespec trispec){
    doRes4_ = true;
    RM4axis_ = RM4axis;
    phi4axis_ = phi4axis;
    trispec_ = trispec;
    if(verbose_>1){
        printf("enabled res4\n");
    }
}

void EECCalculator::run(){
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

    if(verbose_){
        printf("TIME STATISTICS\n");
        long proj = duration_total_ - duration_transferproj_ - duration_res3_ - duration_res4_ - duration_transferres3_ - duration_transferres4_;
        printf("total: %0.3f\n", std::log10(duration_total_));
        printf("\tproj:         %0.3f\n", std::log10(proj));
        printf("\ttransferproj: %0.3f\n", std::log10(duration_transferproj_));
        printf("\tres3:         %0.3f\n", std::log10(duration_res3_));
        printf("\ttransferres3: %0.3f\n", std::log10(duration_transferres3_));
        printf("\tres4:         %0.3f\n", std::log10(duration_res4_));
        printf("\ttransferres4: %0.3f\n", std::log10(duration_transferres4_));
    }
    ran_=true;
}

const std::vector<double> EECCalculator::getproj(unsigned order) const{
    checkRan();
    checkOrder(order);
    
    return projwts_->data(order-2).data();
}

const std::vector<double> EECCalculator::getproj_PU(unsigned order) const{
    checkRan();
    checkOrder(order);
    checkPU();
    
    return projwts_PU_->data(order-2).data();
}

const std::vector<double> EECCalculator::getres3() const{
    checkRan();
    checkOrder(3);
    checkRes3();

    return res3wts_->data(0).data();
}

const std::vector<double> EECCalculator::getres4() const{
    checkRan();
    checkOrder(4);
    checkRes4();

    return res4wts_->data(0).data();
}

const std::vector<double> EECCalculator::getres3_PU() const{
    checkRan();
    checkOrder(3);
    checkRes3();
    checkPU();

    return res3wts_PU_->data(0).data();
}

const std::vector<double> EECCalculator::getres4_PU() const{
    checkRan();
    checkOrder(4);
    checkRes4();
    checkPU();

    return res4wts_PU_->data(0).data();
}

const arma::mat EECCalculator::getTransferproj(unsigned order) const {
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

const arma::mat EECCalculator::getTransferres3() const {
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

const arma::mat EECCalculator::getTransferres4() const {
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

bool EECCalculator::hasRun() const{
    return ran_;
}

unsigned EECCalculator::getMaxOrder() const {
    return maxOrder_;
}

void EECCalculator::initialize(){
    if(verbose_){
        printf("top of initialize\n");
    }

    //double sum=0;
    //printf("particle momenta:\n");
    //for(const auto& p : J1E_){
    //    printf("%0.5g\n", p);
    //    sum += p;
    //}
    //printf("sumpt: %0.5g\n", sum);

    ran_ = false;

    unsigned Nacc = maxOrder_ - 1;
    projwts_ = std::make_shared<EECweightAccumulator>(Nacc, *RLaxis_);
    projwts_->setupMaxDRIndexer(J1_);
    if(verbose_>2){
        printf("initialized projwts\n");
    }
    
    if(doRes3_){
        res3wts_ = std::make_shared<EECweightAccumulator>(1, axisvec({*RLaxis_, *xi3axis_, 
                                                                        *phi3axis_}));
        res3wts_->setupThirdOrderIndexer(J1_);
    }
    if(verbose_>2){
        printf("initialized res3wts\n");
    }

    if(doRes4_){
        res4wts_ = std::make_shared<EECweightAccumulator>(1, axisvec({*RLaxis_, *RM4axis_, 
                                                                        *phi4axis_}));
        res4wts_->setupFourthOrderIndexer(J1_, trispec_);
    }
    if(verbose_>2){
        printf("initialized res4\n");
    }

    if(doPU_){
        projwts_PU_ = std::make_shared<EECweightAccumulator>(Nacc, *RLaxis_);
        projwts_PU_->setupMaxDRIndexer(J1_);

        if(doRes3_){
            res3wts_PU_ = std::make_shared<EECweightAccumulator>(1, axisvec({*RLaxis_, 
                                                                            *xi3axis_, 
                                                                            *phi3axis_}));
            res3wts_PU_->setupThirdOrderIndexer(J1_);
        }

        if(doRes4_){
            res4wts_PU_ = std::make_shared<EECweightAccumulator>(1, axisvec({*RLaxis_, 
                                                                            *RM4axis_, 
                                                                        *phi4axis_}));
            res4wts_PU_->setupFourthOrderIndexer(J1_, trispec_);
        }
    }
    if(verbose_>2){
        printf("initialized PUwts\n");
    }

    if(doTrans_){
        projtrans_ = std::make_shared<EECtransferAccumulator>(Nacc, *RLaxis_);
        projtrans_->setupMaxDRIndexers(J2_, J1_);

        if(doRes3_){
            res3trans_ = std::make_shared<EECtransferAccumulator>(1, axisvec({*RLaxis_,
                                                                            *xi3axis_,
                                                                            *phi3axis_}));
            res3trans_->setupThirdOrderIndexers(J2_, J1_);
        }

        if(doRes4_){
            res4trans_ = std::make_shared<EECtransferAccumulator>(1, axisvec({*RLaxis_, 
                                                                            *RM4axis_,
                                                                        *phi4axis_}));
            res4trans_->setupFourthOrderIndexers(J2_, J1_, trispec_);
        }
    }
    if(verbose_>2){
        printf("initialized transfer\n");
    }

    if(verbose_){
        printf("end of initialize\n");
    }
}

arma::mat EECCalculator::makePtrans(const jet& genjet, 
                                    const jet& recojet,
                                    const arma::mat& rawmat, 
                                    const normType& norm) {
    arma::mat ans(rawmat);

    arma::vec genpt = genjet.ptvec();
    arma::vec recopt = recojet.ptvec();
    double normfactorG = getNormFact(genjet, norm);
    double normfactorR = getNormFact(recojet, norm);
    genpt/=normfactorG;
    recopt/=normfactorR;

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

void EECCalculator::checkPU() const{
    if(!doPU_){
        throw std::logic_error("Error! Asking for PU results without doing PU calculation");
    }
}

void EECCalculator::checkRes3() const{
    if(!doRes3_){
        throw std::logic_error("Error! Asking for res3 results without doing res3 calculation");
    }
}

void EECCalculator::checkRes4() const{
    if(!doRes4_){
        throw std::logic_error("Error! Asking for res4 resoluts without doing res4 calculation");
    }
}

void EECCalculator::checkOrder(unsigned order) const{
    if(order >maxOrder_){
        throw std::logic_error("Error! Asking for results of higher order than computed");
    } 
    if(order <2){
        throw std::logic_error("Observable not defined or computed for order < 2");
    }
}

void EECCalculator::checkRan() const{
    if(!ran_){
        throw std::logic_error("Error! Asking for EEC results before running EEC calculation");
    }
}

void EECCalculator::checkTransfer() const {
    if(!doTrans_){
        throw std::logic_error("Error! Asking for EEC transfer results when no second jet was supplied!");
    }
}

void EECCalculator::computeMwayContribution(unsigned M){
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

void EECCalculator::accumulateWt(const unsigned M,
                  const uvec& ord){

    auto start = std::chrono::high_resolution_clock::now();
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
                nextWt *= intPow(J1E_.at(iP), c.composition.at(i));
                if (doPU_ && PU_[ord.at(i)]){
                    hasPU = true;
                }
            }
            wts[order-M] += nextWt;
            //if(order==5){
            //    printf("Weight for ");
            //    printOrd(ord);
            //    printf(" is %0.5g\n", nextWt);
            //}
            if (hasPU){
                wts_PU[order-M] += nextWt;
            }


            if (doTrans_){
                auto beforetrans = std::chrono::high_resolution_clock::now();
                projtrans_->accumulate(ord, c.composition, 
                        adj_, ptrans_, nextWt, order-2);
                auto aftertrans = std::chrono::high_resolution_clock::now();
                duration_transferproj_ += std::chrono::duration_cast<std::chrono::nanoseconds>(aftertrans-beforetrans).count();
            }
            if(doRes3_ && order==3){
                auto beforeres3 = std::chrono::high_resolution_clock::now();
                res3wts_->accumulate(ord, c.composition, {nextWt});
                if(hasPU){
                    res3wts_PU_->accumulate(ord, c.composition, 
                            {nextWt});
                }
                auto afterres3 = std::chrono::high_resolution_clock::now();
                duration_res3_ += std::chrono::duration_cast<std::chrono::nanoseconds>(afterres3-beforeres3).count();
                if (doTrans_){
                    auto beforeres3trans = std::chrono::high_resolution_clock::now();
                    res3trans_->accumulate(ord, c.composition,
                            adj_, ptrans_, nextWt, 0);
                    auto afterres3trans = std::chrono::high_resolution_clock::now();
                    duration_transferres3_ += std::chrono::duration_cast<std::chrono::nanoseconds>(afterres3trans-beforeres3trans).count();
                }
            } 

            if(doRes4_ && order==4){
                auto beforeres4 = std::chrono::high_resolution_clock::now();
                res4wts_->accumulate(ord, c.composition, {nextWt});
                if(hasPU){
                    res4wts_PU_->accumulate(ord, c.composition,
                            {nextWt});
                }
                auto afterres4 = std::chrono::high_resolution_clock::now();
                duration_res4_ += std::chrono::duration_cast<std::chrono::nanoseconds>(afterres4-beforeres4).count();
                if (doTrans_){
                    auto beforeres4trans = std::chrono::high_resolution_clock::now();
                    res4trans_->accumulate(ord, c.composition,
                            adj_, ptrans_, nextWt, 0);
                    auto afterres4trans = std::chrono::high_resolution_clock::now();
                    duration_transferres4_ += std::chrono::duration_cast<std::chrono::nanoseconds>(afterres4trans-beforeres4trans).count();
                }
            }
        }//end for each composition
    }//end for each order
    //if(M<=5){
    //    printf("wtsvec: \n");
    //    for(auto w : wts){
    //        printf("\t%0.5g \n", w);
    //    }
    //}
    projwts_->accumulate(ord, uvec(M, 1), wts, M-2); 
    if(doPU_){
        projwts_PU_->accumulate(ord, uvec(M, 1), wts_PU, M-2);
    }
    auto end = std::chrono::high_resolution_clock::now();
    duration_total_ += std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count();
}

void EECCalculator::computePointAtZero(){
    auto start = std::chrono::high_resolution_clock::now();
    for(unsigned i=0; i<J1_.nPart; ++i){//for each particle
        std::vector<double> wts, wts_PU;
        wts.reserve(maxOrder_-1);
        if(doPU_){
            wts_PU.reserve(maxOrder_-1);
        }
        for(unsigned order=2; order<=maxOrder_; ++order){//for each order
            double nextwt;
            nextwt = intPow(J1E_.at(i), order);

            //accumulate weight
            wts.push_back(nextwt); 
            if (doPU_ && PU_.at(i)){
                wts_PU.push_back(nextwt);
            }
        }//end for each order
        projwts_->accumulate(uvec({i}), uvec({1}), wts);
        if(doPU_){
            projwts_PU_->accumulate(uvec({i}), uvec({1}), wts_PU);
        }

        if(doRes3_){
            auto beforeres3 = std::chrono::high_resolution_clock::now();
            res3wts_->accumulate(uvec({i}), uvec({3}), 
                                 std::vector<double>({wts[1]}));
            if(doRes3_){
                res3wts_PU_->accumulate(uvec({i}), uvec({3}),
                                        std::vector<double>({wts_PU[1]}));
            }
            auto afterres3 = std::chrono::high_resolution_clock::now();
            duration_res3_ += std::chrono::duration_cast<std::chrono::nanoseconds>(afterres3-beforeres3).count();
        }
        if(doRes4_){
            auto beforeres4 = std::chrono::high_resolution_clock::now();
            res4wts_->accumulate(uvec({i}), uvec({4}), 
                                 std::vector<double>({wts[2]}));
            if(doRes4_){
                res4wts_PU_->accumulate(uvec({i}), uvec({4}), 
                                        std::vector<double>({wts_PU[2]}));
            }
            auto afterres4 = std::chrono::high_resolution_clock::now();
            duration_res4_ += std::chrono::duration_cast<std::chrono::nanoseconds>(afterres4-beforeres4).count();
        }
        if(doTrans_){
            auto beforetrans = std::chrono::high_resolution_clock::now();
            for(unsigned order=2; order<=maxOrder_; ++order){
                projtrans_->accumulate(uvec({i}), uvec({order}), adj_, ptrans_,
                                       wts[order-2], order-2);
            }
            auto aftertrans = std::chrono::high_resolution_clock::now();
            duration_transferproj_ = std::chrono::duration_cast<std::chrono::nanoseconds>(aftertrans-beforetrans).count();
            if(doRes3_){
                auto beforetrans3 = std::chrono::high_resolution_clock::now();
                res3trans_->accumulate(uvec({i}), uvec({3}), adj_, ptrans_,
                                      wts[1], 0);
                auto aftertrans3 = std::chrono::high_resolution_clock::now();
                duration_transferres3_ = std::chrono::duration_cast<std::chrono::nanoseconds>(aftertrans3-beforetrans3).count();

            }
            if(doRes4_){
                auto beforetrans4 = std::chrono::high_resolution_clock::now();
                res4trans_->accumulate(uvec({i}), uvec({4}), adj_, ptrans_,
                                      wts[2], 0);
                auto aftertrans4 = std::chrono::high_resolution_clock::now();
                duration_transferres4_ = std::chrono::duration_cast<std::chrono::nanoseconds>(aftertrans4-beforetrans4).count();
            }
        }
    }//end for each particle
    auto end = std::chrono::high_resolution_clock::now();
    duration_total_ += std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count();
}

std::vector<double> EECCalculator::normalizePt(const jet& j,
                                               const normType& norm){
    std::vector<double> ans;
    ans.reserve(j.nPart);
    double normFact = getNormFact(j, norm);
    for(unsigned i=0; i<j.nPart; ++i){
        ans.push_back(j.particles.at(i).pt/normFact);
    }
    return ans;
}

double EECCalculator::getNormFact(const jet& j, const normType& norm){
    switch (norm){
        case RAWPT:
            return j.rawpt;
        case SUMPT:
            return j.sumpt;
        case CORRPT:
            return j.pt;
        default:
            throw std::invalid_argument("Invalid normType");
    }
}
