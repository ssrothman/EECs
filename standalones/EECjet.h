#ifndef SROTHMAN_EECS_STANDALONE_EECJET_H
#define SROTHMAN_EECS_STANDALONE_EECJET_H

#include <boost/multi_array.hpp>
#include <vector>

#include "usings.h"
#include "SRothman/SimonTools/src/simon_jet.h"

namespace standaloneEEC{
    enum normType{
        NONE,
        SUMPT,
        CORRPT,
        RAWPT
    };

    class EECjet {
    public:
        struct oneentry{
            double E;
            double eta;
            double phi;
        };

        size_t N;
        vector<oneentry> singles;

        EECjet(const simon_jet& J, const normType& nt);
    };

    class EECjet_precomputed : public EECjet {
    public:
        struct pairentry {
            double deta;
            double dphi;
            double floatDR;
        };
        multi_array<pairentry, 2> pairs;

        EECjet_precomputed(const simon_jet& J, 
                           const normType& nt) noexcept ;
    };
};

#endif
