#include "jetinfo.h"

std::vector<double> normalizePt(const jet& j){
    std::vector<double> result;
    result.reserve(j.nPart);
    
    for(const particle& part: j.particles){
        result.emplace_back(part.pt/j.sumpt);
    }
    return result;
}

vecND::nodiagvec getdR2s(const jet& j){
    vecND::nodiagvec result(j.nPart, 2u);

    std::vector<unsigned> ord = result.ord0();
    unsigned i=0;
    do{
        const particle& p1 = j.particles[ord[0]];
        const particle& p2 = j.particles[ord[1]];
        result.at(i) = dR2(p1.eta, p1.phi,
                           p2.eta, p2.phi);
        ++i;
    } while(result.iterate(ord));
    return result;
}
