#ifndef SROTHMAN_EECS_STANDALONE_EECJET_H
#define SROTHMAN_EECS_STANDALONE_EECJET_H

#include <boost/multi_array.hpp>
#include <vector>

#include "usings.h"
#include "SRothman/SimonTools/src/jet.h"
#include "SRothman/SimonTools/src/deltaR.h"
#include "SRothman/SimonTools/src/particlePair.h"

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
        double pt;
        double eta;
        double phi;
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

    template <bool dsq>
    class PrecomputedPairs {
    public:
        static constexpr bool distances_squared = dsq;

        PrecomputedPairs(const Singles& singles) noexcept{
            pairs.resize(extents[singles.size()][singles.size()]);

            for (unsigned i=0; i<singles.size(); ++i){
                for (unsigned j=0; j<singles.size(); ++j){

                    const oneentry& singleI = singles.get(i);
                    const oneentry& singleJ = singles.get(j);

                    pairs[i][j] = simon::pairentry_resolved<dsq>(
                            singleI, singleJ);
                }
            }
        }

        inline const simon::pairentry_resolved<dsq>& get(const unsigned i, const unsigned j) const noexcept{
            //printf("PrecomputedPairs::get(%u, %u)\n",i,j);
            //fflush(stdout);
            return pairs[i][j];
        }
    private:
        multi_array<simon::pairentry_resolved<dsq>, 2> pairs;
    };

    template <bool dsq>
    class JITPairs {
    public:
        static constexpr bool distances_squared = dsq;

        JITPairs(const Singles& singles) noexcept : 
            singles(singles) {}

        inline const simon::pairentry_resolved<dsq> get(const unsigned i, const unsigned j) const noexcept {
            //printf("JITPairs::get(%u, %u)\n",i,j);
            //fflush(stdout);
            return simon::pairentry_resolved<dsq>(
                singles.get(i), singles.get(j)
            );
        }
        Singles singles;
    };

    template <class Pairs>
    class EECjet {
    public:
        using pairType = Pairs;

        const size_t N;
        const Singles singles;
        const Pairs pairs;

        EECjet(const simon::jet& J, const normType& nt) :
            N(J.particles.size()),
            singles(J, nt),
            pairs(singles) {}
    };

    typedef EECjet<PrecomputedPairs<false>> EECjet_Precomputed;
    typedef EECjet<JITPairs<true>> EECjet_JIT;
};

#endif
