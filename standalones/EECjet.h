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
        size_t N;
        vector<double> etas, phis, Es;

        EECjet(const simon_jet& J, const normType& nt);
    };

    class EECjet_precomputed : public EECjet {
    public:
        struct pairentry {
            double floatDR;
            double dphi;
            double deta;
            unsigned dRbin;
        };
        multi_array<pairentry, 2> pairs;

        EECjet_precomputed(const simon_jet& J, 
                           const normType& nt,
                           const axis& ax);
    };
};

#endif
