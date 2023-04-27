#ifndef EECs_JETINFO_H
#define EECs_JETINFO_H

#include <vector>
#include <memory>

#include <armadillo>

#include "toyjets/common.h"

#include "simon_util_cpp/util.h"
#include "simon_util_cpp/deltaR.h"
#include "simon_util_cpp/vecND.h"
#include "simon_util_cpp/combinatorics.h"

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
