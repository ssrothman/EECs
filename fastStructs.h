#ifndef EECS_FASTSTRUCTS_H
#define EECS_FASTSTRUCTS_H

#include "adj.h"

#include <boost/histogram.hpp>
#include <boost/multi_array.hpp>

#include "usings.h"
#include "util.h"

#include "SRothman/SimonTools/src/recursive_reduce.h"

namespace fastEEC{
    template <typename T>
    struct res4shapes{
        std::shared_ptr<multi_array<T, 3>> dipole;
        std::shared_ptr<multi_array<T, 3>> tee;
        std::shared_ptr<multi_array<T, 3>> triangle;

        std::shared_ptr<multi_array<T, 3>> minR_phi1;
        std::shared_ptr<multi_array<T, 3>> minR_phi2;
        std::shared_ptr<multi_array<T, 3>> minR_phidiff;
        std::shared_ptr<multi_array<T, 3>> minR_theta;

        void fill(const T val, const unsigned RL,
                  const unsigned shape,
                  const unsigned r, const unsigned theta) noexcept {
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
                    assert(false);
            }
        }

        void fillMinR(const T val, 
                      const unsigned RL, 
                      const unsigned rmax,
                      const unsigned phi1,
                      const unsigned phi2,
                      const unsigned phidiff,
                      const unsigned theta) noexcept{
            (*minR_phi1)[RL][rmax][phi1] += val;
            (*minR_phi2)[RL][rmax][phi2] += val;
            (*minR_phidiff)[RL][rmax][phidiff] += val;
            (*minR_theta)[RL][rmax][theta] += val;
        }

        void fillTri(const T val, const bool isTri,
                const unsigned RL, const unsigned r,
                const unsigned theta) noexcept {
            if(isTri){
                (*triangle)[RL][r][theta] += val;
            } 
        }

        void fill_dipole(const T val, const unsigned RL,
                         const unsigned r, const unsigned theta) noexcept {
            (*dipole)[RL][r][theta] += val;
        }

        void fill_tee(const T val, const unsigned RL,
                      const unsigned r, const unsigned theta) noexcept {
            (*tee)[RL][r][theta] += val;
        }

        void setup(unsigned NRL, 
                   unsigned Nr_dipole, unsigned Ntheta_dipole,
                   unsigned Nr_tee, unsigned Ntheta_tee,
                   unsigned Nr_triangle, unsigned Ntheta_triangle,
                   unsigned Nr_minR, unsigned Nphi_minR,
                   unsigned Nphidiff_minR, unsigned Ntheta_minR) noexcept {
            dipole = std::make_shared<multi_array<T, 3>>(
                    extents[NRL][Nr_dipole][Ntheta_dipole]
            );
            tee = std::make_shared<multi_array<T, 3>>(
                    extents[NRL][Nr_tee][Ntheta_tee]
            );
            triangle = std::make_shared<multi_array<T, 3>>(
                    extents[NRL][Nr_triangle][Ntheta_triangle]
            );
            minR_phi1 = std::make_shared<multi_array<T, 3>>(
                    extents[NRL][Nr_minR][Nphi_minR]
            );
            minR_phi2 = std::make_shared<multi_array<T, 3>>(
                    extents[NRL][Nr_minR][Nphi_minR]
            );
            minR_phidiff = std::make_shared<multi_array<T, 3>>(
                    extents[NRL][Nr_minR][Nphidiff_minR]
            );
            minR_theta = std::make_shared<multi_array<T, 3>>(
                    extents[NRL][Nr_minR][Ntheta_minR]
            );
            
            std::fill(dipole->data(), dipole->data() + dipole->num_elements(), 0);
            std::fill(tee->data(), tee->data() + tee->num_elements(), 0);
            std::fill(triangle->data(), triangle->data() + triangle->num_elements(), 0);
            std::fill(minR_phi1->data(), minR_phi1->data() + minR_phi1->num_elements(), 0);
            std::fill(minR_phi2->data(), minR_phi2->data() + minR_phi2->num_elements(), 0);
            std::fill(minR_phidiff->data(), minR_phidiff->data() + minR_phidiff->num_elements(), 0);
            std::fill(minR_theta->data(), minR_theta->data() + minR_theta->num_elements(), 0);
        }
    };

    template <typename T>
    struct res4shapes_transfer{
        std::shared_ptr<multi_array<T, 6>> dipole;
        std::shared_ptr<multi_array<T, 6>> tee;
        std::shared_ptr<multi_array<T, 6>> triangle;

        void setup(unsigned NRL, 
                   unsigned Nr_dipole, unsigned Ntheta_dipole,
                   unsigned Nr_tee, unsigned Ntheta_tee,
                   unsigned Nr_triangle, unsigned Ntheta_triangle) noexcept {
            dipole = std::make_shared<multi_array<T, 6>>(
                    extents[NRL][Nr_dipole][Ntheta_dipole][NRL][Nr_dipole][Ntheta_dipole]
            );
            tee = std::make_shared<multi_array<T, 6>>(
                    extents[NRL][Nr_tee][Ntheta_tee][NRL][Nr_tee][Ntheta_tee]
            );
            triangle = std::make_shared<multi_array<T, 6>>(
                    extents[NRL][Nr_triangle][Ntheta_triangle][NRL][Nr_triangle][Ntheta_triangle]
            );

            std::fill(dipole->data(), dipole->data() + dipole->num_elements(), 0);
            std::fill(tee->data(), tee->data() + tee->num_elements(), 0);
            std::fill(triangle->data(), triangle->data() + triangle->num_elements(), 0);
        }

        void fillTri(const T val, 
                     const bool isTri_gen,
                     const unsigned RL_gen, 
                     const unsigned r_gen,
                     const unsigned theta_gen,
                     const unsigned RL_reco,
                     [[maybe_unused]] const bool isTri_reco,
                     const unsigned r_reco,
                     const unsigned theta_reco) noexcept {
            if(isTri_gen){
                (*triangle)[RL_gen][r_gen][theta_gen][RL_reco][r_reco][theta_reco] += val;
            }
        }

        void fill(const T val,
                  const unsigned RL_gen,
                  const unsigned shape_gen,
                  const unsigned r_gen,
                  const unsigned ct_gen,
                  const unsigned RL_reco,
                  [[maybe_unused]] const unsigned shape_reco,
                  const unsigned r_reco,
                  const unsigned ct_reco) noexcept {
            /*if (shape_gen != shape_reco){
                return;
            }*/

            switch(shape_gen){
                case 0:
                    break;
                case 1:
                    fill_dipole(val, 
                            RL_gen, r_gen, ct_gen, 
                            RL_reco, r_reco, ct_reco);
                    break;
                case 2:
                    fill_tee(val, 
                            RL_gen, r_gen, ct_gen, 
                            RL_reco, r_reco, ct_reco);
                    break;
                default:
                    assert(false);
            }
        }

        void fill_dipole(const T val,
                         const unsigned RL_gen, 
                         const unsigned r_gen,
                         const unsigned ct_gen,
                         const unsigned RL_reco,
                         const unsigned r_reco,
                         const unsigned ct_reco) noexcept {
            (*dipole)[RL_gen][r_gen][ct_gen][RL_reco][r_reco][ct_reco] += val;
        }

        void fill_tee(const T val,
                      const unsigned RL_gen,
                      const unsigned r_gen,
                      const unsigned ct_gen,
                      const unsigned RL_reco,
                      const unsigned r_reco,
                      const unsigned ct_reco) noexcept {
            (*tee)[RL_gen][r_gen][ct_gen][RL_reco][r_reco][ct_reco] += val;
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

        std::shared_ptr<res4shapes<T>> resolved4_shapes;
        std::shared_ptr<res4shapes<T>> resolved4_shapes_PU;
        std::shared_ptr<res4shapes_transfer<T>> transfer_res4_shapes;

        result_t(const result_t&) = delete;
        result_t() = default;

        result_t& operator+=(const result_t& other) noexcept {
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

            if(other.resolved4_shapes){
                addInPlace(*(resolved4_shapes->dipole), 
                           *(other.resolved4_shapes->dipole));
                addInPlace(*(resolved4_shapes->tee),
                           *(other.resolved4_shapes->tee));
                addInPlace(*(resolved4_shapes->triangle),
                           *(other.resolved4_shapes->triangle));
                addInPlace(*(resolved4_shapes->minR_phi1),
                           *(other.resolved4_shapes->minR_phi1));
                addInPlace(*(resolved4_shapes->minR_phi2),
                           *(other.resolved4_shapes->minR_phi2));
                addInPlace(*(resolved4_shapes->minR_phidiff),
                           *(other.resolved4_shapes->minR_phidiff));
                addInPlace(*(resolved4_shapes->minR_theta),
                           *(other.resolved4_shapes->minR_theta));
            }

            if(other.resolved4_shapes_PU){
                addInPlace(*(resolved4_shapes_PU->dipole),
                           *(other.resolved4_shapes_PU->dipole));
                addInPlace(*(resolved4_shapes_PU->tee),
                           *(other.resolved4_shapes_PU->tee));
                addInPlace(*(resolved4_shapes_PU->triangle),
                           *(other.resolved4_shapes_PU->triangle));
                addInPlace(*(resolved4_shapes_PU->minR_phi1),
                           *(other.resolved4_shapes_PU->minR_phi1));
                addInPlace(*(resolved4_shapes_PU->minR_phi2),
                           *(other.resolved4_shapes_PU->minR_phi2));
                addInPlace(*(resolved4_shapes_PU->minR_phidiff),
                           *(other.resolved4_shapes_PU->minR_phidiff));
                addInPlace(*(resolved4_shapes_PU->minR_theta),
                           *(other.resolved4_shapes_PU->minR_theta));
            }

            if(other.transfer_res4_shapes){
                addInPlace(*(transfer_res4_shapes->dipole), 
                           *(other.transfer_res4_shapes->dipole));
                addInPlace(*(transfer_res4_shapes->tee),
                           *(other.transfer_res4_shapes->tee));
                addInPlace(*(transfer_res4_shapes->triangle),
                           *(other.transfer_res4_shapes->triangle));
            }
            return *this;
        }

        void summarize() const noexcept{
            for(unsigned o=0; o<5; ++o){
                printf("\torder %u: %g\n", o+2, recursive_reduce(*(wts[o]), 0.0));
            }
            if(resolved3){
                printf("\tres3: %g\n", recursive_reduce(*resolved3, 0.0));
            }
            if(resolved4_shapes){
                printf("\tdipole: %g\n", recursive_reduce(*(resolved4_shapes->dipole), 0.0));
                printf("\ttee: %g\n", recursive_reduce(*(resolved4_shapes->tee), 0.0));
                printf("\ttriangle: %g\n", recursive_reduce(*(resolved4_shapes->triangle), 0.0));
            }

            if(wts_PU[0]){
                for(unsigned o=0; o<5; ++o){
                    printf("\torder %u PU: %g\n", o+2, recursive_reduce(*(wts_PU[o]), 0.0));
                }
            }

            if(resolved3_PU){
                printf("\tres3 PU: %g\n", recursive_reduce(*resolved3_PU, 0.0));
            }

            if(resolved4_shapes_PU){
                printf("\tdipole PU: %g\n", recursive_reduce(*(resolved4_shapes_PU->dipole), 0.0));
                printf("\ttee PU: %g\n", recursive_reduce(*(resolved4_shapes_PU->tee), 0.0));
                printf("\ttriangle PU: %g\n", recursive_reduce(*(resolved4_shapes_PU->triangle), 0.0));
            }

            if(transfer_wts[0]){
                for(unsigned o=0; o<5; ++o){
                    double total = recursive_reduce(*(wts[o]), 0.0);
                    printf("\torder %u transfer: %g\n", o+2, recursive_reduce(*(transfer_wts[o]), 0.0));
                }
            }

            if(transfer_res3){
                double total = recursive_reduce(*resolved3, 0.0);
                printf("\tres3 transfer: %g\n", recursive_reduce(*transfer_res3, 0.0));
            }

            if(transfer_res4_shapes){
                double total_dipole = recursive_reduce(*(resolved4_shapes->dipole), 0.0);
                double total_tee = recursive_reduce(*(resolved4_shapes->tee), 0.0);
                double total_triangle = recursive_reduce(*(resolved4_shapes->triangle), 0.0);
                printf("\tdipole transfer: %g\n", recursive_reduce(*(transfer_res4_shapes->dipole), 0.0));
                printf("\ttee transfer: %g\n", recursive_reduce(*(transfer_res4_shapes->tee), 0.0));
                printf("\ttriangle transfer: %g\n", recursive_reduce(*(transfer_res4_shapes->triangle), 0.0));
            }
        }
    };

    template <typename T>
    struct jetDetails_t{
        multi_array<T, 2> floatDRs;
        multi_array<unsigned, 2> dRbins;

        std::vector<T> etas;
        std::vector<T> phis;
        std::vector<T> Es;

        jetDetails_t() noexcept :
            floatDRs(extents[1][1]),
            dRbins(extents[1][1]),
            etas(0),
            phis(0),
            Es(0)
        {}

        jetDetails_t(const jet& J, const axisptr& ax, const normType nt) noexcept:
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

        axisptr r_minR=nullptr;
        axisptr phi_minR=nullptr;
        axisptr phidiff_minR=nullptr;
        axisptr theta_minR=nullptr;

        float shapetol=0;
    };

    template <typename T>
    struct transferInputs{
        std::shared_ptr<jetDetails_t<T>> recoJet;

        std::shared_ptr<adjacency> adj;
        std::shared_ptr<Eigen::MatrixXd> ptrans;

        void setup(const jet * recoJet,
                   const Eigen::MatrixXd * ptrans,
                   const axisptr& ax,
                   const normType nt) noexcept {
            this->recoJet = std::make_shared<jetDetails_t<T>>(*recoJet, ax, nt);
            this->ptrans = std::make_shared<Eigen::MatrixXd>((*ptrans).transpose());
            this->adj = std::make_shared<adjacency>(*ptrans);
        }
    };
}


#endif
