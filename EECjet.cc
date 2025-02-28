#include "EECjet.h"
#include "SRothman/SimonTools/src/histutil.h"
#include "SRothman/SimonTools/src/deltaR.h"

EEC::Singles::Singles(const simon::jet& J, const normType& nt){
    pt_denom = 1;
    switch(nt){
        case NONE:
            pt_denom = 1;
            break;
        case SUMPT:
            pt_denom = J.sumpt;
            break;
        case CORRPT:
            pt_denom = J.pt;
            break;
        case RAWPT:
            pt_denom = J.rawpt;
            break;
        default:
            throw std::invalid_argument("Invalid normType");
    }

    singles.clear();
    singles.resize(J.particles.size());

    for (unsigned i=0; i<J.particles.size(); ++i){
        const auto& p = J.particles[i];
        singles[i].eta = p.eta;
        singles[i].phi = p.phi;
        singles[i].pt = p.pt / pt_denom;
    }
}
