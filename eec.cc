#include "eec.h"
#include "simon_util_cpp/util.h"
#include "simon_util_cpp/deltaR.h"

projectedEEC packageResults(const constdata& cd,
                            computedata& ans){
    using fvec=std::vector<float>;
    auto dRs = std::make_unique<fvec>(cd.dR2s->data());
    dRs->emplace_back(0.0);
    
    auto wts = std::make_unique<fvec>(ans.wts->data());
    wts->emplace_back(ans.ptAtZero);

    if (cd.dR2s_o){
        auto dRs_o = std::make_unique<fvec>(cd.dR2s_o->data());
        dRs_o->emplace_back(0.0);

        return projectedEEC(cd.order, cd.nPart,
                            std::move(dRs),
                            std::move(wts),
                            std::move(dRs_o),
                            std::move(ans.cov),
                            std::move(ans.transfer));

    } else {
        return projectedEEC(cd.order, cd.nPart,
                            std::move(dRs),
                            std::move(wts),
                            nullptr,
                            std::move(ans.cov),
                            nullptr);
    }
}

void finalizeCovariance(const constdata& cd,
                        computedata& ans){
    float normfact;
    for(unsigned iPart=0; iPart<cd.nPart; ++iPart){
        normfact = intPow<float>(1/(1-cd.Es->at(iPart)), cd.order);
        for(unsigned iDR=0; iDR<cd.dR2s->size(); ++iDR){
            float rescaled_contrib = -normfact * (*ans.cov)(iDR, iPart);
            float rescaled_actual = (normfact-1) * ans.wts->at(iDR);
            (*ans.cov)(iDR, iPart) = rescaled_contrib + rescaled_actual;
        }
        unsigned iDR = cd.dR2s->size();
        (*ans.cov)(iDR, iPart) = -normfact * (*ans.cov)(iDR, iPart)
                              +(normfact-1) * ans.ptAtZero;
    }
}

void addTransfer(const constdata& cd,
                 const std::vector<unsigned>& ord_t,
                 const unsigned dRidx,
                 const float nextWt,
                 computedata& ans){

    std::vector<unsigned> ord_o = ord0_full(cd.order);

    do{
        float tfact = 1;
        for(unsigned i=0; i<cd.order && tfact>0; ++i){
            tfact *= cd.ptrans->at(ord_o[i], ord_t[i]);
        }
        if(tfact==0){
            continue;
        }
        unsigned dRidx_o = getMaxDR(*cd.dR2s_o, ord_o);
        ans.transfer->at(dRidx_o, dRidx) += nextWt * tfact;
    } while(iterate_full(cd.order, ord_o, cd.nPart_o));
}


projectedEEC doProjected(const std::shared_ptr<const jet> j,
                         const unsigned order,
                         const std::shared_ptr<const arma::fmat> ptrans,
                         const std::shared_ptr<const jet> j_o){
    constdata cd = getConstdata(j, order, j_o, ptrans);
    computedata ans = initializeComputedata(cd);
    
    pointAtZero(cd, ans);

    for(unsigned M=2; M<=cd.order; ++M){
        MwayContribution(cd, M, ans);
    }

    finalizeCovariance(cd, ans);

    return packageResults(cd, ans);            
}

void pointAtZero(const constdata& cd, computedata& ans){
    float result = 0;
    for(unsigned i=0; i<cd.nPart; ++i){
        float nextwt = intPow<float>(cd.Es->at(i), cd.order);
        result += nextwt;
        //accumulate covariance
        (*ans.cov)(cd.dR2s->size(), i) += nextwt; 

        if(cd.ptrans){
            std::vector<unsigned> ord_t(cd.order, i);
            addTransfer(cd, ord_t, cd.dR2s->size(), nextwt, ans);
        }
    }
    ans.ptAtZero = result;

}

void MwayContribution(const constdata& cd,
                      const unsigned M,
                      computedata& ans){

    std::vector<unsigned> ord = ord0_nodiag(M);
    bool loop;
    do{
        addWt(cd, M, ord, ans);
    } while (iterate_nodiag(M, ord, cd.nPart));
}

void addWt(const constdata& cd,
           const unsigned M,
           const std::vector<unsigned>& ord,
           computedata& ans){
    unsigned dRidx = getMaxDR(*cd.dR2s, ord);


    for(const comp& c : cd.comps->at(M-1)){//for each composition
        double nextWt = c.factor;
        for(unsigned i=0; i<M; ++i){
            nextWt *= intPow<float>(cd.Es->at(ord[i]), c.composition[i]);
        }
        ans.wts->at(dRidx) += nextWt;
        
        //accumulate covariance
        for(unsigned i=0; i<M; ++i){
            (*ans.cov)(dRidx, ord[i]) += nextWt;
        }

        if(cd.ptrans){
            std::vector<unsigned> ord_t;
            ord_t.reserve(cd.order);
            for(unsigned i=0; i<M; ++i){
                ord_t.insert(ord_t.end(), c.composition[i], ord[i]);
            }

            addTransfer(cd, ord_t, dRidx, nextWt, ans);
        }
    }
}

unsigned getMaxDR(const vecND::nodiagvec& dR2s,
                  const std::vector<unsigned> ord){
    if(uniform(ord)){
        return dR2s.size();
    }
    float maxDR=0;
    unsigned maxidx=-1;

    unsigned idx;
    float dR;
    const unsigned M = ord.size();
    std::vector<unsigned> ord2 = ord0_nodiag(2);
    do{
        std::vector<unsigned> ord_test(2);
        ord_test[0] = ord[ord2[0]];
        ord_test[1] = ord[ord2[1]];
        dR = dR2s.at(ord_test, &idx);
        if(dR > maxDR){
            maxDR = dR;
            maxidx = idx;
        }
    } while(iterate_nodiag(2u, ord2, M));
    return maxidx;
}

void getDR2(vecND::nodiagvec& result,
            const std::vector<float>& eta,
            const std::vector<float>& phi){
    //check that arguments make sense
    if(result.dim()!=2){
        throw std::logic_error("getDR() needs a 2-dimensional vecND");
    }
    if(result.nPart()!=eta.size()){
        throw std::logic_error("result.nPart != eta.size()");
    }
    if(eta.size()!=phi.size()){
        throw std::logic_error("eta.size() != phi.size()");
    }

    std::vector<unsigned> ord = result.ord0();
    bool loop;
    for(size_t i=0, loop=true; loop; loop=result.iterate(ord), ++i){//iterate over pairs
        result.at(i) = dR2(eta[ord[0]], phi[ord[0]],
                           eta[ord[1]], phi[ord[1]]);
    }
}

void normalizePt(const std::vector<float> &pts,
                 std::vector<float>& Eout){
    float sumPt=0;
    for(const float& pt: pts){
        sumPt += pt;
    }

    Eout.clear();
    Eout.reserve(pts.size());
    
    for(const float& pt: pts){
        Eout.push_back(pt/sumPt);
    }
}

constdata getConstdata(const std::shared_ptr<const jet> j,
                       const unsigned order,
                       const std::shared_ptr<const jet> j_o,
                       const std::shared_ptr<const arma::fmat> ptrans){
    unsigned nPart = j->pt.size();
    
    auto Es = std::make_unique<std::vector<float>>();
    normalizePt(j->pt, *Es);

    auto dR2s = std::make_unique<vecND::nodiagvec>(nPart, 2u);
    getDR2(*dR2s, j->eta, j->phi);

    auto comps = std::make_unique<comp_t>();
    fillCompositions(order, *comps); 

    if(j_o){
        unsigned nPart_o = j_o->pt.size();
        auto dR2s_o = std::make_unique<vecND::nodiagvec>(nPart_o, 2u);
        getDR2(*dR2s_o, j_o->eta, j_o->phi);
        return constdata(order, nPart,
                         std::move(Es),
                         std::move(dR2s),
                         std::move(comps),
                         std::move(ptrans),
                         nPart_o,
                         std::move(dR2s_o));
    } else {
        return constdata(order, nPart, 
                         std::move(Es), 
                         std::move(dR2s), 
                         std::move(comps));
    }
}

computedata initializeComputedata(const constdata& cd){
    auto wts = std::make_unique<vecND::nodiagvec>(cd.nPart, 2u);
    auto cov = std::make_unique<arma::fmat>(wts->size()+1, cd.nPart,
                                            arma::fill::zeros);

    if(cd.ptrans){
        unsigned size_o = choose(cd.nPart_o, 2u);
        auto transfer  = std::make_unique<arma::fmat>(size_o+1,
                                                     wts->size()+1,
                                                     arma::fill::zeros);
        return computedata(std::move(wts), 
                           std::move(cov),
                           std::move(transfer)); 
    } else {
        return computedata(std::move(wts),
                           std::move(cov));
    }
}
