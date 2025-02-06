#ifndef SROTHMAN_EECS_RES4_STANDALONE_H
#define SROTHMAN_EECS_RES4_STANDALONE_H

#include "standalone_structs.h"

namespace standaloneEEC{
    void res4_standalone_multi_array(
            res4_result_multi_array& res4,
            
            const simon_jet& J,
            const normType& nt,

            const axis& RL,
            const axis& r_dipole,
            const axis& c_dipole,
            const axis& r_tee,
            const axis& c_tee,
            const axis& r_triangle,
            const axis& c_triangle,

            double tolerance);

    void res4_standalone_vector(
            res4_result_vector& res4,
            
            const simon_jet& J,
            const normType& nt,

            const axis& RL,
            const axis& r_dipole,
            const axis& c_dipole,
            const axis& r_tee,
            const axis& c_tee,
            const axis& r_triangle,
            const axis& c_triangle,

            double tolerance);

    /*void res4_standalone_precompute(res4_result_multi_array& res4,
                         
                         const simon_jet& J,
                         const normType& nt,

                         const axis& RL,
                         const axis& r_dipole,
                         const axis& c_dipole,
                         const axis& r_tee,
                         const axis& c_tee,
                         const axis& r_triangle,
                         const axis& c_triangle,

                         double tolerance);*/
};

#endif
