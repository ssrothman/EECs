#ifndef EECS_FASTSTRUCTS_H
#define EECS_FASTSTRUCTS_H

#include "adj.h"

#include <boost/histogram.hpp>
#include <boost/multi_array.hpp>

#include "usings.h"
#include "util.h"

namespace fastEEC{
    template <typename T>
    struct res4shapes{
        std::shared_ptr<multi_array<T, 3>> dipole;
        std::shared_ptr<multi_array<T, 3>> tee;

        void fill(const T val, const unsigned RL,
                  const unsigned shape,
                  const unsigned r, const unsigned theta){
            switch(shape){
                case 0:
                    break;
                case 1:
                    fill_dipole(val, RL, r, theta);
                    break;
                case 2:
                    fill_tee(val, RL, r, theta);
                    break;
                default:
                    printf("SHAPE %u\n", shape);
                    fflush(stdout);
                    assert(false);
            }
        }

        void fill_dipole(const T val, const unsigned RL,
                         const unsigned r, const unsigned theta){
            (*dipole)[RL][r][theta] += val;
        }

        void fill_tee(const T val, const unsigned RL,
                      const unsigned r, const unsigned theta){
            (*tee)[RL][r][theta] += val;
        }

        void setup(unsigned NRL, 
                   unsigned Nr_dipole, unsigned Ntheta_dipole,
                   unsigned Nr_tee, unsigned Ntheta_tee){
            dipole = std::make_shared<multi_array<T, 3>>(
                    extents[NRL][Nr_dipole][Ntheta_dipole]
            );
            tee = std::make_shared<multi_array<T, 3>>(
                    extents[NRL][Nr_tee][Ntheta_tee]
            );
            
            std::fill(dipole->data(), dipole->data() + dipole->num_elements(), 0);
            std::fill(tee->data(), tee->data() + tee->num_elements(), 0);
        }
    };

    template <typename T>
    struct res4shapes_transfer{
        std::shared_ptr<multi_array<T, 6>> dipole;
        std::shared_ptr<multi_array<T, 6>> tee;

        void setup(unsigned NRL, 
                   unsigned Nr_dipole, unsigned Ntheta_dipole,
                   unsigned Nr_tee, unsigned Ntheta_tee){
            dipole = std::make_shared<multi_array<T, 6>>(
                    extents[NRL][Nr_dipole][Ntheta_dipole][NRL][Nr_dipole][Ntheta_dipole]
            );
            tee = std::make_shared<multi_array<T, 6>>(
                    extents[NRL][Nr_tee][Ntheta_tee][NRL][Nr_tee][Ntheta_tee]
            );

            std::fill(dipole->data(), dipole->data() + dipole->num_elements(), 0);
            std::fill(tee->data(), tee->data() + tee->num_elements(), 0);
        }
    };
    
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
         *xi = shortest side/medium side
         *       phi = arcsin(1 - square(RL-RM)/square(RS))
         *
         * cf 2201.07800
         *    2205.02857
         */
        std::shared_ptr<multi_array<T, 3>> resolved3; 
        std::shared_ptr<multi_array<T, 3>> resolved3_PU;
        std::shared_ptr<multi_array<T, 6>> transfer_res3;

        res4shapes<T> resolved4_shapes;
        res4shapes<T> resolved4_shapes_PU;
        res4shapes_transfer<T> transfer_res4_shapes;

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

        result_t& operator+=(const result_t& other){
            for(unsigned i=0; i<5; ++i){
                if(other.wts[i]){
                    addInPlace(*wts[i], *other.wts[i]);
                }
                if(other.wts_PU[i]){
                    addInPlace(*wts_PU[i], *other.wts_PU[i]);
                }
                if(other.transfer_wts[i]){
                    addInPlace(*transfer_wts[i], *other.transfer_wts[i]);
                }
            }

            if (other.resolved3){
                addInPlace(*resolved3, *other.resolved3);
            }
            if (other.resolved3_PU){
                addInPlace(*resolved3_PU, *other.resolved3_PU);
            }
            if (other.transfer_res3){
                addInPlace(*transfer_res3, *other.transfer_res3);
            }

            if(other.resolved4_shapes.dipole){
                addInPlace(*(resolved4_shapes.dipole), 
                           *(other.resolved4_shapes.dipole));
                addInPlace(*(resolved4_shapes.tee),
                           *(other.resolved4_shapes.tee));
            }

            if(other.resolved4_shapes_PU.dipole){
                addInPlace(*(resolved4_shapes_PU.dipole),
                           *(other.resolved4_shapes_PU.dipole));
                addInPlace(*(resolved4_shapes_PU.tee),
                           *(other.resolved4_shapes_PU.tee));
            }
            //if(other.transfer_res4_shapes){
            //    addInPlace(*transfer_res4_shapes, *other.transfer_res4_shapes);
            //}

            if(other.resolved4_fixed){
                addInPlace(*resolved4_fixed, *other.resolved4_fixed);
            }
            if(other.resolved4_fixed_PU){
                addInPlace(*resolved4_fixed_PU, *other.resolved4_fixed_PU);
            }
            if(other.transfer_res4_fixed){
                addInPlace(*transfer_res4_fixed, *other.transfer_res4_fixed);
            }
            return *this;
        }
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
        axisptr RL=nullptr;

        axisptr r_dipole=nullptr;
        axisptr ct_dipole=nullptr;

        axisptr r_tee=nullptr;
        axisptr ct_tee=nullptr;

        axisptr r_triangle=nullptr;
        axisptr ct_triangle=nullptr;

        float shapetol=0;
    };

    struct res4fixedAxes_t{
        axisptr RL=nullptr;

        float shapetol=0;
    };

    template <typename T>
    struct transferInputs{
        jetDetails_t<T> recoJet;

        adjacency adj;
        arma::mat ptrans;
    };
}

#endif
