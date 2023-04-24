#ifndef EECs_COMPS_H
#define EECs_COMPS_H

#include "simon_util_cpp/combinatorics.h"
#include <memory>
#include <stdio.h>

inline std::shared_ptr<std::vector<comp_t>> getCompositions(unsigned order){
    printf("top of getCompositions\n");
    auto result = std::make_shared<std::vector<comp_t>>();
    result->resize(order-1);
    for(unsigned i=2; i<=order; ++i){
        fillCompositions(i, (*result)[i-2]);
    }
    printf("bottom of getCompositions\n");
    return result;
}

std::shared_ptr<std::vector<comp_t>> getCustomComps(unsigned p1, unsigned p2);

#endif
