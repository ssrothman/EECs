#ifndef EECs_MAXDR_H
#define EECs_MAXDR_H

#include "jetinfo.h"

unsigned getMaxDR(const jetinfo& J,
                  std::vector<unsigned>& ord,
                  bool uniqify = false);

unsigned getResolvedDR(const jetinfo& J,
                       std::vector<unsigned>& ord,
                       const std::vector<unsigned>& comp);

#endif
