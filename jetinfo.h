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

std::shared_ptr<std::vector<double>> normalizePt(const jet& j);

std::shared_ptr<vecND::nodiagvec> getdR2s(const jet& j);

struct jetinfo{
    unsigned nPart;
    std::shared_ptr<std::vector<double>> Es;
    std::shared_ptr<vecND::nodiagvec> dR2s;

    jetinfo():
        nPart(0), Es(nullptr), dR2s(nullptr) {
        }

    jetinfo(unsigned nPart,
            std::shared_ptr<std::vector<double>>&& Es,
            std::shared_ptr<vecND::nodiagvec>&& dR2s):
        nPart(nPart), Es(std::move(Es)), dR2s(std::move(dR2s)) {
        }

    jetinfo(const jet& j):
        nPart(j.nPart), 
        Es(normalizePt(j)),
        dR2s(getdR2s(j)) {
        }

    jetinfo(const std::shared_ptr<const jet>& j){
        if(j){
            nPart = j->nPart;
            Es = normalizePt(*j);
            dR2s = getdR2s(*j);
        } else {
            nPart = 0;
            Es = nullptr;
            dR2s = nullptr;
        }
    }
};

#endif
