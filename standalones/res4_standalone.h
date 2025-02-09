#ifndef SROTHMAN_EECS_RES4_STANDALONE_H
#define SROTHMAN_EECS_RES4_STANDALONE_H

#include "standalone_structs.h"

namespace standaloneEEC{
    /*
     * Functions to compute res4 result on a given jet
     * 
     * Two versions provided with different storage backends: 
     *    - multi array dense histogram
     *    - vector of sparse histogram entries
     *
     * Maybe a third version could be added with 
     * actual sparse matrices.....
     *
     * See the implmenetation file for some discussion of
     * the pros and cons of these methods
     *
     * All of the implementation details 
     * (including some template shenanagins)
     * are hidden in the implementation file
     * to avoid as much as possible recomplications
     * and to maintain a clean interface
     */
    void res4_standalone_multi_array(
            res4_result_multi_array& res4,
            
            const simon::jet& J,
            const normType& nt,

            const axis& RL,
            const axis& r_dipole,
            const axis& c_dipole,
            const axis& r_tee,
            const axis& c_tee,
            const axis& r_triangle,
            const axis& c_triangle,

            const double tolerance,
            const double tri_tolerance) noexcept ;

    void res4_standalone_multi_array_precomputed(
            res4_result_multi_array& res4,
            
            const simon::jet& J,
            const normType& nt,

            const axis& RL,
            const axis& r_dipole,
            const axis& c_dipole,
            const axis& r_tee,
            const axis& c_tee,
            const axis& r_triangle,
            const axis& c_triangle,

            const double tolerance,
            const double tri_tolerance) noexcept ;

    void res4_standalone_vector(
            res4_result_vector& res4,
            
            const simon::jet& J,
            const normType& nt,

            const axis& RL,
            const axis& r_dipole,
            const axis& c_dipole,
            const axis& r_tee,
            const axis& c_tee,
            const axis& r_triangle,
            const axis& c_triangle,

            const double tolerance,
            const double tri_tolerance) noexcept ;

    void res4_standalone_vector_precomputed(
            res4_result_vector& res4,
            
            const simon::jet& J,
            const normType& nt,

            const axis& RL,
            const axis& r_dipole,
            const axis& c_dipole,
            const axis& r_tee,
            const axis& c_tee,
            const axis& r_triangle,
            const axis& c_triangle,

            const double tolerance,
            const double tri_tolerance) noexcept ;

    void res4_standalone_transfer_multi_array(
            res4_result_multi_array& ans,
            res4_transfer_result_multi_array& transfer_ans,

            const simon::jet& J1,
            const simon::jet& J2,
            const normType& nt,

            const Eigen::MatrixXd& adjmat,

            const standaloneEEC::axis& ax_R,
            const standaloneEEC::axis& ax_r_dipole,
            const standaloneEEC::axis& ax_c_dipole,
            const standaloneEEC::axis& ax_r_tee,
            const standaloneEEC::axis& ax_c_tee,
            const standaloneEEC::axis& ax_r_triangle,
            const standaloneEEC::axis& ax_c_triangle,

            const double tolerance,
            const double tri_tolerance) noexcept ;

    void res4_standalone_transfer_vector(
            res4_result_multi_array& ans,
            res4_transfer_result_vector& transfer_ans,

            const simon::jet& J1,
            const simon::jet& J2,
            const normType& nt,

            const Eigen::MatrixXd& adjmat,

            const standaloneEEC::axis& ax_R,
            const standaloneEEC::axis& ax_r_dipole,
            const standaloneEEC::axis& ax_c_dipole,
            const standaloneEEC::axis& ax_r_tee,
            const standaloneEEC::axis& ax_c_tee,
            const standaloneEEC::axis& ax_r_triangle,
            const standaloneEEC::axis& ax_c_triangle,

            const double tolerance,
            const double tri_tolerance) noexcept ;
};

#endif
