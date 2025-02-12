#ifndef SROTHMAN_EECS_STANDALONE_EECJET_H
#define SROTHMAN_EECS_STANDALONE_EECJET_H

#include <boost/multi_array.hpp>
#include <vector>

#include "usings.h"
#include "SRothman/SimonTools/src/jet.h"
#include "SRothman/SimonTools/src/deltaR.h"

namespace EEC{
    enum normType{
        NONE,
        SUMPT,
        CORRPT,
        RAWPT
    };

    inline normType normType_from_string(const std::string& s){
        if(s == "NONE") return normType::NONE;
        if(s == "SUMPT") return normType::SUMPT;
        if(s == "CORRPT") return normType::CORRPT;
        if(s == "RAWPT") return normType::RAWPT;
        throw std::invalid_argument("Invalid normType string: " + s);
    }

    struct oneentry{
        double E;
        double eta;
        double phi;
    };

    struct pairentry {
        double deta;
        double dphi;
        double floatDR;
    };

    class Singles {
    public:
        Singles(const simon::jet& J, const normType& nt);

        inline const oneentry& get(const unsigned i) const noexcept {
            //printf("Singles::get(%u)\n",i);
            //fflush(stdout);
            return singles[i];
        }

        size_t size() const noexcept {
            return singles.size();
        }

        double get_pt_denom() const noexcept {
            return pt_denom;
        }

    private:
        std::vector<oneentry> singles;
        double pt_denom;
    };

    class PrecomputedPairs {
    public:
        static constexpr bool distances_squared = false;

        PrecomputedPairs(const Singles& singles) noexcept;

        inline const pairentry& get(const unsigned i, const unsigned j) const noexcept{
            //printf("PrecomputedPairs::get(%u, %u)\n",i,j);
            //fflush(stdout);
            return pairs[i][j];
        }
    private:
        multi_array<pairentry, 2> pairs;
    };

    class JITPairs {
    public:
        static constexpr bool distances_squared = true;

        JITPairs(const Singles& singles) noexcept;

        inline const pairentry get(const unsigned i, const unsigned j) const noexcept {
            //printf("JITPairs::get(%u, %u)\n",i,j);
            //fflush(stdout);
            const oneentry& singleI = singles.get(i);
            const oneentry& singleJ = singles.get(j);

            double dphi = simon::deltaPhi(singleI.phi, 
                                          singleJ.phi);
            double deta = simon::deltaEta(singleI.eta, 
                                          singleJ.eta);
            double dR2 = dphi*dphi + deta*deta;
            return {deta, dphi, dR2};
        }
        Singles singles;
    };

    template <class Pairs>
    class EECjet {
    public:
        const size_t N;
        const Singles singles;
        const Pairs pairs;

        EECjet(const simon::jet& J, const normType& nt) :
            N(J.particles.size()),
            singles(J, nt),
            pairs(singles) {}
    };

    typedef EECjet<PrecomputedPairs> EECjet_Precomputed;
    typedef EECjet<JITPairs> EECjet_JIT;
};

#endif
