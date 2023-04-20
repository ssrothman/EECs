#include "eec.h"
#include "simon_util_cpp/util.h"
#include "simon_util_cpp/deltaR.h"

projectedEEC packageResults(const constdata& cd,
                            computedata& ans){
    using fvec=std::vector<f_t>;
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
    f_t normfact;
    for(unsigned iPart=0; iPart<cd.nPart; ++iPart){
        normfact = intPow<f_t>(1/(1-cd.Es->at(iPart)), cd.order);
        for(unsigned iDR=0; iDR<cd.dR2s->size(); ++iDR){
            f_t rescaled_contrib = -normfact * (*ans.cov)(iDR, iPart);
            f_t rescaled_actual = (normfact-1) * ans.wts->at(iDR);
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
                 const f_t nextWt,
                 computedata& ans){

    std::vector<unsigned> ord_o = ord0_full(cd.order);

    std::vector<unsigned> nadj(cd.order);
    for(unsigned i=0; i<cd.order; ++i){
        nadj[i] = cd.adj->at(ord_t[i]).size();
        if(nadj[i]==0){
            return; //break out early if there are no neighbors
        }
    }

    do{
        float tfact = 1.0f;
        std::vector<unsigned> ord_test(cd.order);
        for(unsigned i=0; i<cd.order; ++i){
            ord_test[i] = cd.adj->at(ord_t[i])[ord_o[i]];
            tfact *= cd.ptrans->at(ord_test[i], ord_t[i]);
        }
        unsigned dRidx_o = getMaxDR(*cd.dR2s_o, ord_test);
        ans.transfer->at(dRidx_o, dRidx) += nextWt * tfact;
    } while (iterate_awkward(nadj, ord_o));
}


projectedEEC doProjected(const std::shared_ptr<const jet> j,
                         const unsigned order,
                         const std::shared_ptr<const mat_t> ptrans,
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
    f_t result = 0;
    for(unsigned i=0; i<cd.nPart; ++i){
        f_t nextwt = intPow<f_t>(cd.Es->at(i), cd.order);
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
            nextWt *= intPow<f_t>(cd.Es->at(ord[i]), c.composition[i]);
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

unsigned getMaxDR(const nodiagvec& dR2s,
                  const std::vector<unsigned> ord){
    if(uniform(ord)){
        return dR2s.size();
    }
    f_t maxDR=0;
    unsigned maxidx=-1;

    unsigned idx;
    f_t dR;
    const unsigned M = ord.size();
    std::vector<unsigned> ord2 = ord0_nodiag(2);
    do{
        std::vector<unsigned> ord_test(2);
        ord_test[0] = ord[ord2[0]];
        ord_test[1] = ord[ord2[1]];
        if(uniform(ord_test)){
            continue;
        }
        dR = dR2s.at(ord_test, &idx);
        if(dR > maxDR){
            maxDR = dR;
            maxidx = idx;
        }
    } while(iterate_nodiag(2u, ord2, M));
    return maxidx;
}

void getDR2(nodiagvec& result,
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
                 std::vector<f_t>& Eout){
    f_t sumPt=0;
    for(const f_t& pt: pts){
        sumPt += pt;
    }

    Eout.clear();
    Eout.reserve(pts.size());
    
    for(const f_t& pt: pts){
        Eout.push_back(pt/sumPt);
    }
}

void getAdj(const mat_t& ptrans, adjacency_t& adj){
    adj.clear();
    adj.resize(ptrans.n_cols);
    for(unsigned i=0; i<ptrans.n_cols; ++i){
        for(unsigned j=0; j<ptrans.n_rows; ++j){
            if(ptrans(j, i)>0){
                adj[i].emplace_back(j);
            }
        }
    }
}

constdata getConstdata(const std::shared_ptr<const jet> j,
                       const unsigned order,
                       const std::shared_ptr<const jet> j_o,
                       const std::shared_ptr<const mat_t> ptrans){
    unsigned nPart = j->pt.size();
    
    auto Es = std::make_unique<std::vector<f_t>>();
    normalizePt(j->pt, *Es);

    auto dR2s = std::make_unique<nodiagvec>(nPart, 2u);
    getDR2(*dR2s, j->eta, j->phi);

    auto comps = std::make_unique<comp_t>();
    fillCompositions(order, *comps); 

    if(j_o){
        unsigned nPart_o = j_o->pt.size();
        auto dR2s_o = std::make_unique<nodiagvec>(nPart_o, 2u);
        getDR2(*dR2s_o, j_o->eta, j_o->phi);

        auto adj = std::make_unique<adjacency_t>();
        getAdj(*ptrans, *adj);

        return constdata(order, nPart,
                         std::move(Es),
                         std::move(dR2s),
                         std::move(comps),
                         std::move(ptrans),
                         std::move(adj),
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
    auto wts = std::make_unique<nodiagvec>(cd.nPart, 2u);
    auto cov = std::make_unique<mat_t>(wts->size()+1, cd.nPart,
                                            arma::fill::zeros);

    if(cd.ptrans){
        unsigned size_o = choose(cd.nPart_o, 2u);
        auto transfer  = std::make_unique<mat_t>(size_o+1,
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
