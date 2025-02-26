#ifndef SROTHMAN_EEC_RES4RESULT_FORWARD_H
#define SROTHMAN_EEC_RES4RESULT_FORWARD_H

namespace EEC{
    template <class A>
    class Res3Result;

    class Res3VectorContainer;
    class Res3MultiArrayContainer;

    using Res3Result_Vector = Res3Result<Res3VectorContainer>;
    using Res3Result_MultiArray = Res3Result<Res3MultiArrayContainer>;
};

#endif
