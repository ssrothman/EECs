#ifndef EECs_COMPS_H
#define EECs_COMPS_H

#include "SRothman/SimonTools/src/combinatorics.h"
#include <memory>
#include <stdio.h>

std::vector<comp_t> getCompositions(unsigned order);

std::vector<comp_t> getCustomComps(unsigned p1, unsigned p2);

#endif
