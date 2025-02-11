#ifndef SROTHMAN_EEC_RES4RESULT_FORWARD_H
#define SROTHMAN_EEC_RES4RESULT_FORWARD_H

namespace EEC{
    template <class A>
    class Res4Result;

    class Res4VectorContainer;
    class Res4MultiArrayContainer;

    using Res4Result_Vector = Res4Result<Res4VectorContainer>;
    using Res4Result_MultiArray = Res4Result<Res4MultiArrayContainer>;
};

#endif
