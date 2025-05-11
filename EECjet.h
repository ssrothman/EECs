#ifndef SROTHMAN_EECS_EECJET_H
#define SROTHMAN_EECS_EECJET_H

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
            return singles[i];
        }

        inline const double& getE(const unsigned i) const noexcept {
            return singles[i].pt;
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

    template <typename T>
    class ProjPairs {
    public:
        static constexpr bool BINNED = std::is_same<T, unsigned>::value;
        static constexpr bool distances_squared = false;

        template <typename AX>
        ProjPairs(const Singles& singles,
                  [[maybe_unused]] const AX& axes) noexcept{
            pairs.resize(extents[singles.size()][singles.size()]);

            for (unsigned i=0; i<singles.size(); ++i){
                for (unsigned j=0; j<singles.size(); ++j){
                    if (i==j) {
                        pairs[i][j] = 0;
                        continue;
                    }
                    const oneentry& singleI = singles.get(i);
                    const oneentry& singleJ = singles.get(j);

                    double dR = simon::deltaR(
                        singleI.eta, singleI.phi,
                        singleJ.eta, singleJ.phi
                    );

                    if constexpr(BINNED){
                        pairs[i][j] = simon::getIndex(dR, axes.R);
                    } else {
                        pairs[i][j] = dR;
                    }
                }
            }
        }

        inline const T& getDR(const unsigned i, const unsigned j) const noexcept{
            return pairs[i][j];
        }

    private:
        multi_array<T, 2> pairs;
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

        template <class AX>
        EECjet(const simon::jet& J, const normType& nt, const AX& axes) :
            N(J.particles.size()),
            singles(J, nt),
            pairs(singles, axes) {}
    };

    typedef EECjet<PrecomputedPairs<false>> EECjet_Precomputed;
    typedef EECjet<JITPairs<true>> EECjet_JIT;
    typedef EECjet<ProjPairs<unsigned>> EECjet_ProjBinned;
    typedef EECjet<ProjPairs<double>> EECjet_ProjUnbinned;
};

#endif
