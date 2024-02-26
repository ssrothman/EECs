#ifndef EECS_FASTSTRUCTS_H
#define EECS_FASTSTRUCTS_H

#include <boost/histogram.hpp>
#include <boost/multi_array.hpp>

namespace fastEEC{
    using namespace boost;
    using namespace std;

    using axis_t = histogram::axis::variable<double>;
    using axisptr = std::shared_ptr<axis_t>;

    using umat = multi_array<unsigned, 2>;
    using umatptr = std::shared_ptr<umat>;
    
    template <typename T>
    struct result{
        vector<T> wts2;
        vector<T> wts3;
        vector<T> wts4;
        vector<T> wts5;
        vector<T> wts6;

        vector<T> wts2_PU;
        vector<T> wts3_PU;
        vector<T> wts4_PU;
        vector<T> wts5_PU;
        vector<T> wts6_PU;

        multi_array<T, 2> transfer2;
        multi_array<T, 2> transfer3;
        multi_array<T, 2> transfer4;
        multi_array<T, 2> transfer5;
        multi_array<T, 2> transfer6;
    };

    enum normType {
        RAWPT, 
        CORRPT,
        SUMPT, 
    };

    template <typename T>
    struct transferInputs{
        umat dRs;
        adjacency adj;
        multi_array<T, 2> ptrans;
    };
};

#endif
