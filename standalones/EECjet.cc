#include "EECjet.h"
#include "SRothman/SimonTools/src/histutil.h"
#include "SRothman/SimonTools/src/deltaR.h"

standaloneEEC::EECjet::EECjet(const simon::jet& J, const normType& nt) {
   double pt_denominator = 1; 
    switch(nt){
        case NONE:
            pt_denominator = 1;
            break;
        case SUMPT:
            pt_denominator = J.sumpt;
            break;
        case CORRPT:
            pt_denominator = J.pt;
            break;
        case RAWPT:
            pt_denominator = J.rawpt;
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
        singles[i].E = p.pt / pt_denominator;
    }

    N = J.particles.size();
}

standaloneEEC::EECjet_precomputed::EECjet_precomputed(
        const simon::jet& J, const normType& nt) noexcept  :
    EECjet(J, nt) {

    pairs.resize(extents[J.nPart][J.nPart]);
    
    for (unsigned iPart=0; iPart < J.nPart; ++iPart){
        for(unsigned jPart=0; jPart < J.nPart; ++jPart){
            pairentry& entry = pairs[iPart][jPart];

            const oneentry& singleI = singles[iPart];
            const oneentry& singleJ = singles[jPart];

            entry.floatDR = simon::deltaR(singleI.eta, singleI.phi, 
                               singleJ.eta, singleJ.phi);
            entry.dphi = simon::deltaPhi(singleI.phi, singleJ.phi);
            entry.deta = simon::deltaEta(singleI.eta, singleJ.eta);
        }
    }
}
