#ifndef SROTHMAN_EEC_RES4TRANSFERRESULT_FORWARD_H
#define SROTHMAN_EEC_RES4TRANSFERRESULT_FORWARD_H

#include "Res4Result_forward.h"

namespace EEC{
    template <class A>
    class Res4TransferResult;

    class Res4TransferVectorContainer;
    class Res4TransferMultiArrayContainer;

    using Res4TransferResult_Vector = Res4TransferResult<Res4TransferVectorContainer>;
    using Res4TransferResult_MultiArray = Res4TransferResult<Res4TransferMultiArrayContainer>;
};

#endif
