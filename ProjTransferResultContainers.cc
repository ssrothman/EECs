#include "ProjTransferResultContainers.h"


EEC::ProjTransferArrayContainer::ProjTransferArrayContainer(const EEC::ProjTransferVectorContainer<unsigned>& other) noexcept :
    nR_reco(other.nR_reco), nR_gen(other.nR_gen){

    for (const auto& entry : other.get_data()) {
        data[entry.iR_reco][entry.iR_gen] += entry.wt;
    }
}

EEC::ProjTransferArrayContainer& EEC::ProjTransferArrayContainer::operator+=(const EEC::ProjTransferVectorContainer<unsigned> other) noexcept {
    for (const auto& entry : other.get_data()) {
        data[entry.iR_reco][entry.iR_gen] += entry.wt;
    }
    return *this;
}

EEC::ProjTransferArrayContainer EEC::ProjTransferArrayContainer::operator-=(const EEC::ProjTransferVectorContainer<unsigned>& other) noexcept{
    for (const auto& entry : other.get_data()) {
        data[entry.iR_reco][entry.iR_gen] -= entry.wt;
    }
    return *this;
}
