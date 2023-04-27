#include "maxDR.h"
#include "simon_util_cpp/deltaR.h"

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
        return J.dR2s.size();
    }
    double maxDR=0;
    unsigned maxidx=-1;

    unsigned idx;
    double dR;
    std::vector<unsigned> ord2 = ord0_nodiag(2);
    do{
        std::vector<unsigned> ord_test(2);
        ord_test[0] = ord[ord2[0]];
        ord_test[1] = ord[ord2[1]];
        if(uniform(ord_test)){
            continue;
        }
        dR = J.dR2s.at(ord_test, &idx);
        if(dR > maxDR){
            maxDR = dR;
            maxidx = idx;
        }
    } while(iterate_nodiag(2u, ord2, M));
    return maxidx;
}
