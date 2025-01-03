#ifndef EECS_FAST_CLEAR_H
#define EECS_FAST_CLEAR_H

namespace fastEEC{
    template <typename T, bool doPU, bool doTransfer>
    void clear(result_t<T>& ans,
               axisptr ax,
               res3axes_t res3ax,
               res4shapesAxes_t res4ax) noexcept {

        unsigned NDR = AXextent(*ax);

        for(unsigned order=0; order < 5; ++order){
            ans.wts[order] = std::make_shared<vector<T>>(NDR, 0);

            if constexpr(doPU){
                ans.wts_PU[order] = std::make_shared<vector<T>>(NDR, 0);
            }

            if constexpr(doTransfer){
                ans.transfer_wts[order] = std::make_shared<multi_array<T, 2>>(
                        extents[NDR][NDR]
                );
                std::fill(ans.transfer_wts[order]->data(), 
                          ans.transfer_wts[order]->data() 
                            + ans.transfer_wts[order]->num_elements(),
                          0);
            }
        }

        if (res3ax.RL){
            unsigned NDR_coarse = AXextent(*res3ax.RL);
            unsigned Nxi = AXextent(*res3ax.xi);
            unsigned Nphi = AXextent(*res3ax.phi);

            ans.resolved3 = std::make_shared<multi_array<T, 3>>(
                    extents[NDR_coarse][Nxi][Nphi]
            );
            std::fill(ans.resolved3->data(), 
                      ans.resolved3->data() 
                        + ans.resolved3->num_elements(),
                      0);

            if constexpr(doPU){
                ans.resolved3_PU = std::make_shared<multi_array<T, 3>>(
                        extents[NDR_coarse][Nxi][Nphi]
                );

                std::fill(ans.resolved3_PU->data(), 
                          ans.resolved3_PU->data() 
                            + ans.resolved3_PU->num_elements(),
                          0);
            }

            if constexpr(doTransfer){
                ans.transfer_res3 = std::make_shared<multi_array<T, 6>>(
                        extents[NDR_coarse][Nxi][Nphi][NDR_coarse][Nxi][Nphi]
                );
                std::fill(ans.transfer_res3->data(), 
                          ans.transfer_res3->data() 
                            + ans.transfer_res3->num_elements(),
                          0);
            }
        }

        if (res4ax.RL){
            unsigned NRL_res4 = AXextent(*res4ax.RL);
            unsigned Nr_dipole_res4 = AXextent(*res4ax.r_dipole);
            unsigned Ntheta_dipole_res4 = AXextent(*res4ax.ct_dipole);
            unsigned Nr_tee_res4 = AXextent(*res4ax.r_tee);
            unsigned Ntheta_tee_res4 = AXextent(*res4ax.ct_tee);
            unsigned Nr_triangle_res4 = AXextent(*res4ax.r_triangle);
            unsigned Ntheta_triangle_res4 = AXextent(*res4ax.ct_triangle);
            unsigned Nr_minR = AXextent(*res4ax.r_minR);
            unsigned Nphi_minR = AXextent(*res4ax.phi_minR);
            unsigned Nphidiff_minR = AXextent(*res4ax.phidiff_minR);
            unsigned Ntheta_minR = AXextent(*res4ax.theta_minR);
            
            ans.resolved4_shapes = std::make_shared<res4shapes<T>>();

            ans.resolved4_shapes->setup(
                    NRL_res4, 
                    Nr_dipole_res4, Ntheta_dipole_res4,
                    Nr_tee_res4, Ntheta_tee_res4,
                    Nr_triangle_res4, Ntheta_triangle_res4,
                    Nr_minR, Nphi_minR,
                    Nphidiff_minR, Ntheta_minR
            );

            if constexpr(doPU){
                ans.resolved4_shapes_PU = std::make_shared<res4shapes<T>>();

                ans.resolved4_shapes_PU->setup(
                        NRL_res4,
                        Nr_dipole_res4, Ntheta_dipole_res4,
                        Nr_tee_res4, Ntheta_tee_res4,
                        Nr_triangle_res4, Ntheta_triangle_res4,
                        Nr_minR, Nphi_minR,
                        Nphidiff_minR, Ntheta_minR
                );
            }

            if constexpr (doTransfer){
                ans.transfer_res4_shapes = std::make_shared<res4shapes_transfer<T>>();

                ans.transfer_res4_shapes->setup(
                        NRL_res4,
                        Nr_dipole_res4, Ntheta_dipole_res4,
                        Nr_tee_res4, Ntheta_tee_res4,
                        Nr_triangle_res4, Ntheta_triangle_res4
                );
            }
        }
    }//end clear()
}

#endif
