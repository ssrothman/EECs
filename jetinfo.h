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

std::unique_ptr<std::vector<double>> normalizePt(const jet& j);

std::unique_ptr<vecND::nodiagvec> getdR2s(const jet& j);

struct jetinfo{
    unsigned nPart;
    std::unique_ptr<std::vector<double>> Es;
    std::unique_ptr<vecND::nodiagvec> dR2s;

    jetinfo():
        nPart(0), Es(nullptr), dR2s(nullptr) {
            printf("constructed empty jetinfo\n");
        }

    jetinfo(unsigned nPart,
            std::unique_ptr<std::vector<double>>&& Es,
            std::unique_ptr<vecND::nodiagvec>&& dR2s):
        nPart(nPart), Es(std::move(Es)), dR2s(std::move(dR2s)) {
            printf("constructed jetinfo from args\n");
        }

    jetinfo(const jet& j):
        nPart(j.nPart), 
        Es(normalizePt(j)),
        dR2s(getdR2s(j)) {
            printf("constructed jetinfo from jet\n");
        }

    jetinfo(const std::shared_ptr<const jet>& j){
        printf("constructing jetinfo from jet ptr\n");
        if(j){
            nPart = j->nPart;
            Es = normalizePt(*j);
            dR2s = getdR2s(*j);
        } else {
            nPart = 0;
            Es = nullptr;
            dR2s = nullptr;
        }
        printf("constructed jetinfo from jet ptr\n");
    }
};

#endif
