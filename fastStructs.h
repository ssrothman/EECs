#ifndef EECS_FASTSTRUCTS_H
#define EECS_FASTSTRUCTS_H

#include "adj.h"

#include <boost/histogram.hpp>
#include <boost/multi_array.hpp>

#include "usings.h"
#include "util.h"

namespace fastEEC{
    
    template <typename T>
    struct result_t{
        /*
         * Projected weights for orders 2-6
         */
        std::array<std::shared_ptr<vector<T>>, 5> wts;
        std::array<std::shared_ptr<vector<T>>, 5> wts_PU;
        std::array<std::shared_ptr<multi_array<T, 2>>, 5> transfer_wts;

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

        result_t(const result_t&) = delete;
        result_t() = default;
    };

    template <typename T>
    struct jetDetails_t{
        multi_array<T, 2> floatDRs;
        multi_array<unsigned, 2> dRbins;

        std::vector<T> etas;
        std::vector<T> phis;
        std::vector<T> Es;

        jetDetails_t():
            floatDRs(extents[0][0]),
            dRbins(extents[0][0]),
            etas(0),
            phis(0),
            Es(0)
        {}

        jetDetails_t(const jet& J, const axisptr& ax, const normType nt):
            jetDetails_t()
        {
            getFloatDRs(floatDRs, J);
            getDRbins(dRbins, J, ax);
            getEtasPhis(etas, phis, J);
            getEs(Es, J, nt);
        }
    };

    struct res3axes_t{
        axisptr RL;
        axisptr xi;
        axisptr phi;
    };

    struct res4shapesAxes_t{
        axisptr RL;

        axisptr r_dipole;
        axisptr ct_dipole;

        axisptr r_tee;
        axisptr ct_tee;

        axisptr r_triangle;
        axisptr ct_triangle;

        float shapetol;
    };

    struct res4fixedAxes_t{
        axisptr RL;

        float shapetol;
    };

    template <typename T>
    struct transferInputs{
        jetDetails_t<T> recoJet;

        adjacency adj;
        multi_array<T, 2> ptrans;
    };
}

#endif
