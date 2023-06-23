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

std::vector<double> normalizePt(const jet& j);

vecND::nodiagvec getdR2s(const jet& j);

struct jetinfo{
    unsigned nPart;
    std::vector<double> Es;
    vecND::nodiagvec dR2s;

    jetinfo(): nPart(0), Es(), dR2s() {}

    jetinfo(unsigned nPart,
            std::vector<double>&& Es,
            vecND::nodiagvec&& dR2s):
        nPart(nPart), Es(std::move(Es)), dR2s(std::move(dR2s)) {}

    jetinfo(const jet& j):
        nPart(j.nPart), 
        Es(normalizePt(j)),
        dR2s(getdR2s(j)) {}
};

#endif
