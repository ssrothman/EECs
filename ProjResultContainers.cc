#include "ProjResultContainers.h"

EEC::ProjArrayContainer::ProjArrayContainer(const EEC::ProjVectorContainer<unsigned>& other) noexcept :
    nR(other.nR),
    data(nR, 0.0) {
        
        for (const auto& entry : other.get_data()) {
            data[entry.iR] += entry.wt;
        }
}

EEC::ProjArrayContainer& EEC::ProjArrayContainer::operator+=(const EEC::ProjVectorContainer<unsigned> other) noexcept {
    for (const auto& entry : other.get_data()) {
        data[entry.iR] += entry.wt;
    }
    return *this;
}

EEC::ProjArrayContainer& EEC::ProjArrayContainer::operator-=(const EEC::ProjVectorContainer<unsigned>& other) noexcept {
    for (const auto& entry : other.get_data()) {
        data[entry.iR] -= entry.wt;
    }
    return *this;
}
