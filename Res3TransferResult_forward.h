#ifndef SROTHMAN_EEC_RES4TRANSFERRESULT_FORWARD_H
#define SROTHMAN_EEC_RES4TRANSFERRESULT_FORWARD_H

#include "Res3Result_forward.h"

namespace EEC{
    template <class A>
    class Res3TransferResult;

    class Res3TransferVectorContainer;
    class Res3TransferMultiArrayContainer;

    using Res3TransferResult_Vector = Res3TransferResult<Res3TransferVectorContainer>;
    using Res3TransferResult_MultiArray = Res3TransferResult<Res3TransferMultiArrayContainer>;
};

#endif
