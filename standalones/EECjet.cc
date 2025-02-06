#include "EECjet.h"
#include "SRothman/SimonTools/src/histutil.h"
#include "SRothman/SimonTools/src/deltaR.h"

standaloneEEC::EECjet::EECjet(const simon_jet& J, const normType& nt) {
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

    etas.clear();
    phis.clear();
    Es.clear();

    etas.reserve(J.particles.size());
    phis.reserve(J.particles.size());
    Es.reserve(J.particles.size());

    for (const auto& p : J.particles){
        etas.push_back(p.eta);
        phis.push_back(p.phi);
        Es.push_back(p.pt / pt_denominator);
    }

    N = J.particles.size();
}

standaloneEEC::EECjet_precomputed::EECjet_precomputed(
        const simon_jet& J, const normType& nt,
        const axis& ax) :
    EECjet(J, nt) {

    pairs.resize(extents[J.nPart][J.nPart]);
    
    for (unsigned iPart=0; iPart < J.nPart; ++iPart){
        for(unsigned jPart=0; jPart < J.nPart; ++jPart){
            pairentry& entry = pairs[iPart][jPart];

            entry.floatDR = deltaR(etas[iPart], phis[iPart], 
                               etas[jPart], phis[jPart]);
            entry.dphi = deltaPhi(phis[iPart], phis[jPart]);
            entry.deta = deltaEta(etas[iPart], etas[jPart]);
            entry.dRbin = getIndex(entry.floatDR, ax);
        }
    }
}
