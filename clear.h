#ifndef EECS_CLEAR_H
#define EECS_CLEAR_H

namespace fastEEC{
    template <typename T, bool doPU, bool doTransfer, unsigned maxOrder, bool doRes3, bool doRes4, bool doRes4Fixed>
    void clear(result<T>& ans, 
               const unsigned NDR,
               const res3axes_t& res3ax,
               const res4shapesAxes_t& res4ax,
               const res4fixedAxes_t& res4fixedax){
               
               

        if constexpr (doRes4Fixed){
            unsigned Nfixedshape = 3;

            ans.resolved4_fixed=make_shared<multi_array<T, 2>>(
                    extents[Nfixedshape][NDR]
            );
            ans.resolved4_fixed_PU=make_shared<multi_array<T, 2>>(
                    extents[Nfixedshape][NDR]
            );
            ans.transfer_res4_fixed=make_shared<multi_array<T, 4>>(
                    extents[Nfixedshape][NDR][Nfixedshape][NDR]
            );
        
            for(unsigned i=0; i<Nfixedshape; ++i){
                for(unsigned j=0; j<NDR; ++j){
                    (*ans.resolved4_fixed)[i][j] = 0;
                    (*ans.resolved4_fixed_PU)[i][j] = 0;
                    for(unsigned a=0; a<Nfixedshape; ++a){
                        for(unsigned b=0; b<NDR; ++b){
                            (*ans.transfer_res4_fixed)[i][j][a][b] = 0;
                        }
                    }
                }
            }
        }

        if constexpr (doRes4){
            unsigned Nshape = 4;
            unsigned NDR_coarse = histogram::axis::traits::extent(*rin.coarseRL);

            unsigned Nr1 = histogram::axis::traits::extent(*rin.r_dipole);
            unsigned Nr2 = histogram::axis::traits::extent(*rin.r_tee);
            unsigned Nr3 = histogram::axis::traits::extent(*rin.r_triangle);
            assert(Nr1 == Nr2 && Nr2 == Nr3);

            unsigned Nct1 = histogram::axis::traits::extent(*rin.ct_dipole);
            unsigned Nct2 = histogram::axis::traits::extent(*rin.ct_tee);
            unsigned Nct3 = histogram::axis::traits::extent(*rin.ct_triangle);
            assert(Nct1 == Nct2 && Nct2==Nct3);

            ans.resolved4_shapes = make_shared<multi_array<T, 4>>(
                    extents[Nshape][NDR_coarse][Nr1][Nct1]
            );
            ans.resolved4_shapes_PU = make_shared<multi_array<T, 4>>(
                    extents[Nshape][NDR_coarse][Nr1][Nct1]
            );
            ans.transfer_res4_shapes = make_shared<multi_array<T, 8>>(
                    extents[Nshape][NDR_coarse][Nr1][Nct1][Nshape][NDR_coarse][Nr1][Nct1]
            );

            for(unsigned i=0; i<Nshape; ++i){
                for(unsigned j=0; j<NDR_coarse; ++j){
                    for(unsigned k=0; k<Nr1; ++k){
                        for(unsigned l=0; l<Nct1; ++l){
                            (*ans.resolved4_shapes)[i][j][k][l] = 0;
                            (*ans.resolved4_shapes_PU)[i][j][k][l] = 0;
                            for(unsigned a=0; a<Nshape; ++a){
                                for(unsigned b=0; b<NDR_coarse; ++b){
                                    for(unsigned c=0; c<Nr1; ++c){
                                        for(unsigned d=0; d<Nct1; ++d){
                                            (*ans.transfer_res4_shapes)[i][j][k][l][a][b][c][d] = 0;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        if constexpr (doRes3){
            unsigned NDR_coarse = histogram::axis::traits::extent(*rin.coarseRL);
            unsigned Nxi = histogram::axis::traits::extent(*rin.xi);
            unsigned Nphi = histogram::axis::traits::extent(*rin.phi);

            ans.resolved3 = make_shared<multi_array<T, 3>>(
                    extents[NDR_coarse][Nxi][Nphi]
            );
            ans.resolved3_PU = make_shared<multi_array<T, 3>>(
                    extents[NDR_coarse][Nxi][Nphi]
            );
            ans.transfer_res3 = make_shared<multi_array<T, 6>>(
                    extents[NDR_coarse][Nxi][Nphi][NDR_coarse][Nxi][Nphi]
            );

            for(unsigned i=0; i<NDR_coarse; ++i){
                for(unsigned j=0; j<Nxi; ++j){
                    for(unsigned k=0; k<Nxi; ++k){
                        (*ans.resolved3)[i][j][k] = 0;
                        (*ans.resolved3_PU)[i][j][k] = 0;
                        for(unsigned a=0; a<NDR_coarse; ++a){
                            for(unsigned b=0; b<Nxi; ++b){
                                for(unsigned c=0; c<Nxi; ++c){
                                    (*ans.transfer_res3)[i][j][k][a][b][c] = 0;
                                }
                            }
                        }
                    }
                }
            }
        }

        ans.wts2 = make_shared<vector<T>>(NDR);
        fill(ans.wts2->begin(), ans.wts2->end(), 0);
        ans.wts3 = make_shared<vector<T>>(NDR);
        fill(ans.wts3->begin(), ans.wts3->end(), 0);
        ans.wts4 = make_shared<vector<T>>(NDR);
        fill(ans.wts4->begin(), ans.wts4->end(), 0);
        ans.wts5 = make_shared<vector<T>>(NDR);
        fill(ans.wts5->begin(), ans.wts5->end(), 0);
        ans.wts6 = make_shared<vector<T>>(NDR);
        fill(ans.wts6->begin(), ans.wts6->end(), 0);

        if constexpr (doPU){
            ans.wts2_PU = make_shared<vector<T>>(NDR);
            fill(ans.wts2_PU->begin(), ans.wts2_PU->end(), 0);
            ans.wts3_PU = make_shared<vector<T>>(NDR);
            fill(ans.wts3_PU->begin(), ans.wts3_PU->end(), 0);
            ans.wts4_PU = make_shared<vector<T>>(NDR);
            fill(ans.wts4_PU->begin(), ans.wts4_PU->end(), 0);
            ans.wts5_PU = make_shared<vector<T>>(NDR);
            fill(ans.wts5_PU->begin(), ans.wts5_PU->end(), 0);
            ans.wts6_PU = make_shared<vector<T>>(NDR);
            fill(ans.wts6_PU->begin(), ans.wts6_PU->end(), 0);
        }

        if constexpr (doTransfer){
            ans.transfer2 = make_shared<multi_array<T, 2>>(extents[NDR][NDR]);
            ans.transfer3 = make_shared<multi_array<T, 2>>(extents[NDR][NDR]);
            ans.transfer4 = make_shared<multi_array<T, 2>>(extents[NDR][NDR]);
            ans.transfer5 = make_shared<multi_array<T, 2>>(extents[NDR][NDR]);
            ans.transfer6 = make_shared<multi_array<T, 2>>(extents[NDR][NDR]);

            for(unsigned i=0; i<NDR; ++i){
                for(unsigned j=0; j<NDR; ++j){
                    (*ans.transfer2)[i][j] = 0;
                    (*ans.transfer3)[i][j] = 0;
                    (*ans.transfer4)[i][j] = 0;
                    (*ans.transfer5)[i][j] = 0;
                    (*ans.transfer6)[i][j] = 0;
                }
            }
        }
    }
};

#endif
