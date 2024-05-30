#ifndef EECS_FAST_CLEAR_H
#define EECS_FAST_CLEAR_H

namespace fastEEC{
    template <typename T>
    void clear(result_t<T>& ans,
               axisptr ax,
               res3axes_t res3ax,
               res4shapesAxes_t res4ax,
               res4fixedAxes_t res4fixedax){

        unsigned NDR = histogram::axis::traits::extent(*ax);

        for(unsigned order=0; order < 5; ++order){
            ans.wts[order] = std::make_shared<vector<T>>(NDR, 0);
            ans.wts_PU[order] = std::make_shared<vector<T>>(NDR, 0);
            ans.transfer_wts[order] = std::make_shared<multi_array<T, 2>>(
                    extents[NDR][NDR]
            );
            std::fill(ans.transfer_wts[order]->data(), 
                      ans.transfer_wts[order]->data() 
                        + ans.transfer_wts[order]->num_elements(),
                      0);
        }

        if (res3ax.RL){
            unsigned NDR_coarse = histogram::axis::traits::extent(*res3ax.RL);
            unsigned Nxi = histogram::axis::traits::extent(*res3ax.xi);
            unsigned Nphi = histogram::axis::traits::extent(*res3ax.phi);

            ans.resolved3 = std::make_shared<multi_array<T, 3>>(
                    extents[NDR_coarse][Nxi][Nphi]
            );
            std::fill(ans.resolved3->data(), 
                      ans.resolved3->data() 
                        + ans.resolved3->num_elements(),
                      0);

            ans.resolved3_PU = std::make_shared<multi_array<T, 3>>(
                    extents[NDR_coarse][Nxi][Nphi]
            );
            std::fill(ans.resolved3_PU->data(), 
                      ans.resolved3_PU->data() 
                        + ans.resolved3_PU->num_elements(),
                      0);

            ans.transfer_res3 = std::make_shared<multi_array<T, 6>>(
                    extents[NDR_coarse][Nxi][Nphi][NDR_coarse][Nxi][Nphi]
            );
            std::fill(ans.transfer_res3->data(), 
                      ans.transfer_res3->data() 
                        + ans.transfer_res3->num_elements(),
                      0);
        }

        if (res4ax.RL){
            unsigned Nshape = 4;
            unsigned NRL = histogram::axis::traits::extent(*res4ax.RL);
            unsigned Nr = histogram::axis::traits::extent(*res4ax.r_dipole);
            unsigned Nct = histogram::axis::traits::extent(*res4ax.ct_dipole);

            ans.resolved4_shapes = std::make_shared<multi_array<T, 4>>(
                    extents[Nshape][NRL][Nr][Nct]
            );
            std::fill(ans.resolved4_shapes->data(), 
                      ans.resolved4_shapes->data() 
                        + ans.resolved4_shapes->num_elements(),
                      0);

            ans.resolved4_shapes_PU = std::make_shared<multi_array<T, 4>>(
                    extents[Nshape][NRL][Nr][Nct]
            );
            std::fill(ans.resolved4_shapes_PU->data(), 
                      ans.resolved4_shapes_PU->data() 
                        + ans.resolved4_shapes_PU->num_elements(),
                      0);

            ans.transfer_res4_shapes = std::make_shared<multi_array<T, 8>>(
                    extents[Nshape][NRL][Nr][Nct][Nshape][NRL][Nr][Nct]
            );
            std::fill(ans.transfer_res4_shapes->data(), 
                      ans.transfer_res4_shapes->data() 
                        + ans.transfer_res4_shapes->num_elements(),
                      0);
        }

        if (res4fixedax.RL){
            unsigned Nfixedshape = 3;
            unsigned NRL = histogram::axis::traits::extent(*res4fixedax.RL);

            ans.resolved4_fixed = std::make_shared<multi_array<T, 2>>(
                    extents[Nfixedshape][NRL]
            );
            std::fill(ans.resolved4_fixed->data(), 
                      ans.resolved4_fixed->data() 
                        + ans.resolved4_fixed->num_elements(),
                      0);

            ans.resolved4_fixed_PU = std::make_shared<multi_array<T, 2>>(
                    extents[Nfixedshape][NRL]
            );
            std::fill(ans.resolved4_fixed_PU->data(), 
                      ans.resolved4_fixed_PU->data() 
                        + ans.resolved4_fixed_PU->num_elements(),
                      0);

            ans.transfer_res4_fixed = std::make_shared<multi_array<T, 4>>(
                    extents[Nfixedshape][NRL][Nfixedshape][NRL]
            );
            std::fill(ans.transfer_res4_fixed->data(), 
                      ans.transfer_res4_fixed->data() 
                        + ans.transfer_res4_fixed->num_elements(),
                      0);
        }
    }
}

#endif
