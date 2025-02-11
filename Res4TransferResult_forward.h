#ifndef SROTHMAN_EEC_RES4TRANSFERRESULT_FORWARD_H
#define SROTHMAN_EEC_RES4TRANSFERRESULT_FORWARD_H

#include "Res4Result_forward.h"

namespace EEC{
    template <class A, class B>
    class Res4TransferResult;

    class Res4TransferVectorContainer;
    class Res4TransferMultiArrayContainer;

    using Res4TransferResult_Vector_Vector = Res4TransferResult<Res4TransferVectorContainer, Res4VectorContainer>;
    using Res4TransferResult_Vector_MultiArray = Res4TransferResult<Res4TransferVectorContainer, Res4MultiArrayContainer>;
    using Res4TransferResult_MultiArray_Vector = Res4TransferResult<Res4TransferMultiArrayContainer, Res4VectorContainer>;
    using Res4TransferResult_MultiArray_MultiArray = Res4TransferResult<Res4TransferMultiArrayContainer, Res4MultiArrayContainer>;
};

#endif
