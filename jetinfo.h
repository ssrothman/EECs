#ifndef EECs_JETINFO_H
#define EECs_JETINFO_H

#include <vector>
#include <memory>

#include <armadillo>

#include "SRothman/SimonTools/src/jets.h"
#include "SRothman/SimonTools/src/util.h"
#include "SRothman/SimonTools/src/deltaR.h"
#include "SRothman/SimonTools/src/vecND.h"
#include "SRothman/SimonTools/src/combinatorics.h"

std::vector<double> normalizePt(const jet& j, bool toRaw);

vecND::nodiagvec getdRs(const jet& j);

template <typename Axis>
vecND::inodiagvec getdRidxs(const vecND::nodiagvec& dRs, 
                            const Axis& axis){
    vecND::inodiagvec result(dRs.nPart(), dRs.dim());

    std::vector<unsigned> ord = result.ord0();
    unsigned i=0;
    do{
        result.at(i) = axis.index(dRs.at(i)) + 1; //+1 for underflow bin
        ++i;
    } while(result.iterate(ord));
    return result;
}     

struct jetinfo{
    unsigned nPart;
    std::vector<double> Es;
    vecND::nodiagvec dRs;
    vecND::inodiagvec dRidxs;

    jetinfo(): nPart(0), Es(), dRs(), dRidxs() {}

    jetinfo(unsigned nPart, 
            const std::vector<double>& Es,
            const vecND::nodiagvec& dRs, 
            const vecND::inodiagvec& dRidxs):
        nPart(nPart), Es(Es), dRs(dRs), 
        dRidxs(std::move(dRidxs)) {}

    jetinfo(unsigned nPart, 
            std::vector<double>&& Es,
            vecND::nodiagvec&& dRs, 
            vecND::inodiagvec&& dRidxs):
        nPart(nPart), 
        Es(std::move(Es)), 
        dRs(std::move(dRs)), 
        dRidxs(std::move(dRidxs)) {}

    template <typename Axis>
    jetinfo(const jet& j, const Axis& axis, bool normToRaw):
        nPart(j.nPart), 
        Es(normalizePt(j, normToRaw)),
        dRs(getdRs(j)),
        dRidxs(getdRidxs(dRs, axis)) {}

};

/*    jetinfo(): nPart(0), Es(), dRs() {}

    jetinfo(unsigned nPart,
            std::vector<double>&& Es,
            vecND::nodiagvec&& dRs,
            vecND::inodiagvec&& dRidxs):
        nPart(nPart), Es(std::move(Es)), dRs(std::move(dRs), 
                dRidxs(std::move(dRidxs)) {}

    };*/

#endif
