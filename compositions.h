#ifndef EECs_COMPS_H
#define EECs_COMPS_H

#include "simon_util_cpp/combinatorics.h"
#include <memory>

inline std::shared_ptr<comp_t> getCompositions(unsigned order){
    auto result = std::make_shared<comp_t>();
    fillCompositions(order, *result);
    return result;
}

#endif
