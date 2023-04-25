#include "jetinfo.h"

std::unique_ptr<std::vector<double>> normalizePt(const jet& j){
    auto result = std::make_unique<std::vector<double>>();
    result->reserve(j.nPart);
    
    for(const particle& part: j.particles){
        result->emplace_back(part.pt/j.sumpt);
    }
    return result;
}

std::unique_ptr<vecND::nodiagvec> getdR2s(const jet& j){
    auto result = std::make_unique<vecND::nodiagvec>(j.nPart, 2u);

    std::vector<unsigned> ord = result->ord0();
    unsigned i=0;
    do{
        const particle& p1 = j.particles[ord[0]];
        const particle& p2 = j.particles[ord[1]];
        result->at(i) = dR2(p1.eta, p1.phi,
                           p2.eta, p2.phi);
        ++i;
    } while(result->iterate(ord));
    return result;
}
