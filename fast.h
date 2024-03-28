#ifndef EECS_FAST_H

#include <boost/histogram.hpp>
#include <boost/multi_array.hpp>

#include <cassert>

#include "SRothman/SimonTools/src/jets.h"

#include "fastStructs.h"
#include "faststart.h"

namespace fastEEC{
    double getNormFact(const jet& J, const normType& nt){
        switch (nt){
            case RAWPT:
                return J.rawpt;
            case SUMPT:
                return J.sumpt;
            case CORRPT:
                return J.pt;
            default:
                throw std::invalid_argument("Invalid normType");
        }
    }

    template <typename T>
    void getEs(vector<T>& ans, const jet& J, const normType nt){
        ans.clear();
        ans.reserve(J.nPart);

        double normFact = getNormFact(J, nt);
        for(unsigned i=0; i<J.nPart; ++i){
            ans.push_back(J.particles.at(i).pt/normFact);
        }
    }

    template <typename T>
    void getEtasPhis(vector<T>& etas, vector<T>& phis, const jet& J){
        etas.clear();
        phis.clear();

        etas.reserve(J.nPart);
        phis.reserve(J.nPart);

        for(unsigned i=0; i<J.nPart; ++i){
            etas.push_back(J.particles.at(i).eta);
            phis.push_back(J.particles.at(i).phi);
        }
    }

    void getDRs(umat& ans, const jet& J, const axisptr& ax){
        ans.resize(extents[J.nPart][J.nPart]);
        for (unsigned i0=0; i0<J.nPart; ++i0){
            for (unsigned i1=0; i1<J.nPart; ++i1){
                double deltaR = dR(J.particles[i0], J.particles[i1]);
                if (deltaR < 1e-10){
                    deltaR = 0;
                }
                unsigned idx = static_cast<unsigned>(ax->index(deltaR) + 1);
                ans[i0][i1] = idx;
            }
        }
    }

    template <typename T>
    void getFloatDRs(multi_array<T, 2>& ans, const jet& J){
        ans.resize(extents[J.nPart][J.nPart]);
        for (unsigned i0=0; i0<J.nPart; ++i0){
            for (unsigned i1=0; i1<J.nPart; ++i1){
                T deltaR = dR(J.particles[i0], J.particles[i1]);
                if(deltaR < 1e-10){
                    deltaR = 0;
                }
                ans[i0][i1] = deltaR;
            }
        }
    }

    template <typename T>
    void getPtrans(multi_array<T, 2>& ans, const arma::mat& ptrans){//NB we transpose for faster iteration
        ans.resize(extents[ptrans.n_cols][ptrans.n_rows]);
        for(unsigned i=0; i<ptrans.n_rows; ++i){
            for(unsigned j=0; j<ptrans.n_cols; ++j){
                ans[j][i] = ptrans(i,j);
            }
        }
    }

    template <typename T, bool doPU, bool doTransfer, bool doRes3, bool doRes4, bool doRes4Fixed>
    void fastEEC(result<T>& ans,

                      const jet& J, const axisptr& ax, 
                      const int order, const normType nt,

                      axisptr& coarseRLax,
                      axisptr& xiax,
                      axisptr& phiax,

                      axisptr& rax_dipole,
                      axisptr& ctax_dipole,

                      axisptr& rax_tee,
                      axisptr& ctax_tee,

                      axisptr& rax_triangle,
                      axisptr& ctax_triangle,
                      T shapetol,

                      const std::vector<bool>* const PU = nullptr,
                      const jet * const J_Reco = nullptr,
                      const arma::mat* ptrans = nullptr){
        assert(order >= 2 && order <= 6);

        umat dRs; 
        std::vector<T> Es;

        getDRs(dRs, J, ax);
        getEs<T>(Es, J, nt);

        unsigned NDR = histogram::axis::traits::extent(*ax);

        struct transferInputs<T> tin;
        if constexpr (doTransfer){
            getDRs(tin.dRs, *J_Reco, ax);
            tin.adj = adjacency(*ptrans);
            getPtrans<T>(tin.ptrans, *ptrans);

            getFloatDRs(tin.rin.floatDRs, *J_Reco);
            getEtasPhis(tin.rin.etas, tin.rin.phis, *J_Reco);

            tin.rin.coarseRL = coarseRLax;

            tin.rin.xi = xiax;
            tin.rin.phi = phiax;

            tin.rin.r_dipole = rax_dipole;
            tin.rin.ct_dipole = ctax_dipole;

            tin.rin.r_tee = rax_tee;
            tin.rin.ct_tee = ctax_tee;

            tin.rin.r_triangle = rax_triangle;
            tin.rin.ct_triangle = ctax_triangle;

            tin.rin.shapetol = shapetol;
        }

        struct resolvedInputs<T> rin;
        getFloatDRs(rin.floatDRs, J);
        getEtasPhis(rin.etas, rin.phis, J);
        rin.coarseRL = coarseRLax;

        rin.xi = xiax;
        rin.phi = phiax;

        rin.r_dipole = rax_dipole;
        rin.ct_dipole = ctax_dipole;

        rin.r_tee = rax_tee;
        rin.ct_tee = ctax_tee;

        rin.r_triangle = rax_triangle;
        rin.ct_triangle = ctax_triangle;

        rin.shapetol = shapetol;

        if (order == 2){
            start<T, doPU, doTransfer, 2, doRes3, doRes4, doRes4Fixed>(
                    dRs, Es, NDR, rin, ans, PU, &tin
            );
        } else if(order == 3){
            start<T, doPU, doTransfer, 3, doRes3, doRes4, doRes4Fixed>(
                    dRs, Es, NDR, rin, ans, PU, &tin
            );
        } else if(order == 4){
            start<T, doPU, doTransfer, 4, doRes3, doRes4, doRes4Fixed>(
                    dRs, Es, NDR, rin, ans, PU, &tin
            );
        } else if(order == 5){
            start<T, doPU, doTransfer, 5, doRes3, doRes4, doRes4Fixed>(
                    dRs, Es, NDR, rin, ans, PU, &tin
            );
        } else if(order == 6){
            start<T, doPU, doTransfer, 6, doRes3, doRes4, doRes4Fixed>(
                    dRs, Es, NDR, rin, ans, PU, &tin
            );
        }
    }
};

#endif
