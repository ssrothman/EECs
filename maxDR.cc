#include "maxDR.h"
#include "SRothman/SimonTools/src/deltaR.h"

unsigned getResolvedDR(const jetinfo& J,
                       std::vector<unsigned>& ord,
                       const std::vector<unsigned>& comp){
    if(uniform(ord)){
        return 0;
    }
    std::vector<unsigned> ord_full;
    for(unsigned i=0; i<ord.size(); ++i){
        for(unsigned j=0; j<comp[i]; ++j){
            ord_full.push_back(ord[i]);
        }
    }

    unsigned M = ord_full.size();

    std::vector<unsigned> idxs;

    std::vector<unsigned> ord2 = ord0_nodiag(2);
    do{
        unsigned idx;
        std::vector<unsigned> ord_test(2);
        ord_test[0] = ord_full[ord2[0]];
        ord_test[1] = ord_full[ord2[1]];

        if(uniform(ord_test)){
            idx=0;
        } else {
            idx = J.dRidxs.at(ord_test);
        }
        idxs.push_back(idx);
    } while(iterate_nodiag(2u, ord2, M));

    std::sort(idxs.begin(), idxs.end());

    unsigned result=0;
    for(unsigned i=0; i<idxs.size(); ++i){
        result = (J.nPart * result) + idxs[i];
    }

    return result;
}

unsigned getMaxDR(const jetinfo& J,
                  std::vector<unsigned>& ord,
                  bool uniqify){
    unsigned M = ord.size();
    if(uniqify){
        std::sort(ord.begin(), ord.end());
        auto end = std::unique(ord.begin(), ord.end());
        M = std::distance(ord.begin(), end);
    }
    if(uniform(ord)){
        return 0;
    }

    double maxDR=0;
    unsigned maxidx=-1;

    double dR;
    std::vector<unsigned> ord2 = ord0_nodiag(2);
    do{
        std::vector<unsigned> ord_test(2);
        ord_test[0] = ord[ord2[0]];
        ord_test[1] = ord[ord2[1]];
        if(uniform(ord_test)){
            continue;
        }
        dR = J.dRs.at(ord_test);
        if(dR > maxDR){
            maxDR = dR;
            maxidx = J.dRidxs.at(ord_test);
        }
    } while(iterate_nodiag(2u, ord2, M));
    return maxidx;
}
