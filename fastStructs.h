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
        /*
         * Projected weights for orders 2-6
         */
        std::shared_ptr<vector<T>> wts2;
        std::shared_ptr<vector<T>> wts3;
        std::shared_ptr<vector<T>> wts4;
        std::shared_ptr<vector<T>> wts5;
        std::shared_ptr<vector<T>> wts6;

        std::shared_ptr<vector<T>> wts2_PU;
        std::shared_ptr<vector<T>> wts3_PU;
        std::shared_ptr<vector<T>> wts4_PU;
        std::shared_ptr<vector<T>> wts5_PU;
        std::shared_ptr<vector<T>> wts6_PU;

        std::shared_ptr<multi_array<T, 2>> transfer2;
        std::shared_ptr<multi_array<T, 2>> transfer3;
        std::shared_ptr<multi_array<T, 2>> transfer4;
        std::shared_ptr<multi_array<T, 2>> transfer5;
        std::shared_ptr<multi_array<T, 2>> transfer6;
        
        /*
         * shape [RL, xi, phi]
         *
         * where RL = longest side
         *       xi = shortest side/medium side
         *       phi = arcsin(1 - square(RL-RM)/square(RS))
         *
         * cf 2201.07800
         *    2205.02857
         */
        std::shared_ptr<multi_array<T, 3>> resolved3; 
        std::shared_ptr<multi_array<T, 3>> resolved3_PU;
        std::shared_ptr<multi_array<T, 6>> transfer_res3;

        /*
         * shape [shapeindex, Rl, r, theta]
         *
         * where shapeindex is:
         *      0: no special shape
         *         in this case r, theta are both zero
         *      1: dipole
         *         in this case theta is the cross angle, r is the short distance over the long distance
         *      2: tee
         *         in this case r is the short distance over the long distance , theta is the angle at the top of the T
         *      3: triangle
         *         in this case r is the distance from the top of the triangle to the point, over RL
         *         theta is the angle between that vertex and the point
         */
        std::shared_ptr<multi_array<T, 4>> resolved4_shapes; 
        std::shared_ptr<multi_array<T, 4>> resolved4_shapes_PU;
        std::shared_ptr<multi_array<T, 8>> transfer_res4_shapes;

        /*
         * Some fixed shapes for 4th order and 5th order
         * shape [shapeindex, RL]
         * where shapeindex is:
         *     0: no special shape
         *     1: square
         *     2: triangle (for fourth-order) or pentagon (for fifth-order)
         */
        std::shared_ptr<multi_array<T, 2>> resolved4_fixed;
        std::shared_ptr<multi_array<T, 2>> resolved4_fixed_PU;
        std::shared_ptr<multi_array<T, 4>> transfer_res4_fixed;

        //multi_array<T, 2> resolved5_fixed;
        //multi_array<T, 2> resolved5_fixed_PU;
        //multi_array<T, 4> transfer_res5_fixed;
    };

    enum normType {
        RAWPT, 
        CORRPT,
        SUMPT, 
    };

    template <typename T>
    struct resolvedInputs{
        multi_array<T, 2> floatDRs;
        std::vector<T> etas;
        std::vector<T> phis;

        axisptr coarseRL;
        axisptr xi;
        axisptr phi;

        axisptr r_dipole;
        axisptr ct_dipole;

        axisptr r_tee;
        axisptr ct_tee;

        axisptr r_triangle;
        axisptr ct_triangle;

        T shapetol;
    };

    template <typename T>
    struct transferInputs{
        umat dRs;

        adjacency adj;
        multi_array<T, 2> ptrans;

        resolvedInputs<T> rin;
    };
};

#endif
