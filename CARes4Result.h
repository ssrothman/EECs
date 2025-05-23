#ifndef SROTHMAN_EECS_CA_RES4_RESULT_H
#define SROTHMAN_EECS_CA_RES4_RESULT_H

#include "usings.h"
#include "CARes4Axes.h"
#include "ResResultContainers.h"

#include "SRothman/SimonTools/src/histutil.h"

#include <vector>
#include <numeric>
#include <boost/multi_array.hpp>

namespace EEC{
    template <typename T>
    struct CAres4_entry{
        bool is_symmetric=false;
        bool is_chain=false;

        T R_idx=0;

        T symmetric_wrtr_r_idx=0, symmetric_wrtr_c_idx=0;

        T symmetric_wrtR1_r_idx=0, symmetric_wrtR1_c_idx=0;
        T symmetric_wrtR2_r_idx=0, symmetric_wrtR2_c_idx=0;

        T chain_r_idx=0, chain_c_idx=0;

        template <bool SHOULD_BIN>
        void fill_chain(T RLidx,
                        double RL, double RS1, double c,
                        const EEC::CARes4Axes& axes) noexcept {
            is_symmetric = false;
            is_chain = true;

            R_idx = RLidx;

            double rs1 = RL==0 ? 0 : RS1/RL;
            if constexpr (SHOULD_BIN){
                chain_r_idx = simon::getIndex(
                        rs1, axes.r_chain
                );
                chain_c_idx = simon::getIndex(
                        c, axes.c_chain
                );
            } else {
                chain_r_idx = rs1;
                chain_c_idx = c;
            }
        }

        template <bool SHOULD_BIN>
        void fill_symmetric(T RLidx, 
                            double RL, double RS1, double RS2,
                            double cRr1, double cRr2, double cr1r2,
                            const EEC::CARes4Axes& axes) noexcept {
            is_symmetric = true;
            is_chain = false;

            R_idx = RLidx;

            double rs1 = RL==0 ? 0 : RS1/RL;
            double rs2 = RL==0 ? 0 : RS2/RL;

            if constexpr (SHOULD_BIN){
                symmetric_wrtr_r_idx = simon::getIndex(
                        rs1, axes.r_symmetric
                );
                symmetric_wrtr_c_idx = simon::getIndex(
                        cr1r2, axes.c_symmetric
                );

                symmetric_wrtR1_r_idx = simon::getIndex(
                        rs1, axes.r_symmetric
                );
                symmetric_wrtR1_c_idx = simon::getIndex(
                        cRr1, axes.c_symmetric
                );
                symmetric_wrtR2_r_idx = simon::getIndex(
                        rs2, axes.r_symmetric
                );
                symmetric_wrtR2_c_idx = simon::getIndex(
                        cRr2, axes.c_symmetric
                );
            } else {
                symmetric_wrtr_r_idx = rs1;
                symmetric_wrtr_c_idx = cr1r2;
                symmetric_wrtR1_r_idx = rs1;
                symmetric_wrtR1_c_idx = cRr1;
                symmetric_wrtR2_r_idx = rs2;
                symmetric_wrtR2_c_idx = cRr2;
            }
        }
    };

    template <class Container>
    class CARes4Result{
    public:
        CARes4Result() : 
            chain(),
            symmetric_wrtR(),
            symmetric_wrtr(),
            pt_denom_set(false),
            pt_denom(-1) {}

        constexpr static bool SHOULD_BIN = Container::SHOULD_BIN;
        using T = typename Container::TYPE;

        CARes4Result(
                const size_t nR, 
                const size_t nr_chain,
                const size_t nc_chain,
                const size_t nr_symmetric, 
                const size_t nc_symmetric) noexcept :
            chain(nR, nr_chain, nc_chain),
            symmetric_wrtR(nR, nr_symmetric, nc_symmetric),
            symmetric_wrtr(nR, nr_symmetric, nc_symmetric),
            pt_denom_set(false),
            pt_denom(-1) {}

        CARes4Result(const CARes4Axes& axes) noexcept :
                CARes4Result(
                        simon::AXextent(axes.R),
                        simon::AXextent(axes.r_chain),
                        simon::AXextent(axes.c_chain),
                        simon::AXextent(axes.r_symmetric),
                        simon::AXextent(axes.c_symmetric)) {}
        
        template <typename CALCULATOR>
        CARes4Result(const CALCULATOR& calculator) noexcept :
            CARes4Result(calculator.get_axes()) {}

        CARes4Result(const Container& chain, 
                   const Container& symmetric_wrtR, 
                   const Container& symmetric_wrtr) noexcept :
            chain(chain),
            symmetric_wrtR(symmetric_wrtR),
            symmetric_wrtr(symmetric_wrtr),
            pt_denom_set(false),
            pt_denom(-1) {}

        CARes4Result(Container&& chain, 
                   Container&& symmetric_wrtR, 
                   Container&& symmetric_wrtr) noexcept :
            chain(std::move(chain)),
            symmetric_wrtR(std::move(symmetric_wrtR)),
            symmetric_wrtr(std::move(symmetric_wrtr)),
            pt_denom_set(false),
            pt_denom(-1) {}

        template <typename T>
        void fill_chain(const CAres4_entry<T>& entry, double wt) noexcept {
            chain.fill(
                entry.R_idx, 
                entry.chain_r_idx, 
                entry.chain_c_idx, 
                wt
            );
        }

        template <typename T>
        void fill_symmetric(const CAres4_entry<T>& entry, double wt) noexcept {
            symmetric_wrtR.fill(
                entry.R_idx, 
                entry.symmetric_wrtR1_r_idx, 
                entry.symmetric_wrtR1_c_idx,
                wt
            );
            symmetric_wrtR.fill(
                entry.R_idx, 
                entry.symmetric_wrtR2_r_idx, 
                entry.symmetric_wrtR2_c_idx,
                wt
            );
            symmetric_wrtr.fill(
                entry.R_idx, 
                entry.symmetric_wrtr_r_idx, 
                entry.symmetric_wrtr_c_idx,
                wt
            );
        }

        template <typename T>
        void fill(const CAres4_entry<T>& entry, double wt) noexcept {
            if (entry.is_chain){
                fill_chain(entry, wt);
            } else {
                fill_symmetric(entry, wt);
            }
        }

        void set_pt_denom(double pt_denom_) {
            if (pt_denom_set){
                throw std::runtime_error("pt_denom already set");
            } else {
                pt_denom = pt_denom_;
                pt_denom_set = true;
            }
        }

        double get_pt_denom() const {
            if (!pt_denom_set){
                throw std::runtime_error("pt_denom not set");
            } else {
                return pt_denom;
            }
        }

        const Container& get_chain() const noexcept{
            return chain;
        }

        const Container& get_symmetric_wrtR() const noexcept{
            return symmetric_wrtR;
        }

        const Container& get_symmetric_wrtr() const noexcept{
            return symmetric_wrtr;
        }

        double total_chain_weight() const noexcept{
            return chain.total_weight();
        }

        double total_symmetric_wrtR_weight() const noexcept{
            return symmetric_wrtR.total_weight();
        }

        double total_symmetric_wrtr_weight() const noexcept{
            return symmetric_wrtr.total_weight();
        }

        template <class OtherContainer>
        bool operator==(const CARes4Result<OtherContainer>& other) const noexcept{
            return chain == other.get_chain() &&
                   symmetric_wrtR == other.get_symmetric_wrtR() &&
                   symmetric_wrtr == other.get_symmetric_wrtr();
        }

        template <class OtherContainer>
        bool operator!=(const CARes4Result<OtherContainer>& other) const noexcept{
            return !(*this == other);
        }

        template <class OtherContainer>
        CARes4Result<Container>& operator+=(const CARes4Result<OtherContainer>& other) noexcept{
            chain += other.get_chain();
            symmetric_wrtR += other.get_symmetric_wrtR();
            symmetric_wrtr += other.get_symmetric_wrtr();
            return *this;
        }

        template <class OtherContainer>
        CARes4Result<Container> operator+(const CARes4Result<OtherContainer>& other) const noexcept{
            CARes4Result<Container> result(*this);
            result += other;
            return result;
        }

        template <class OtherContainer>
        CARes4Result<Container>& operator-=(const CARes4Result<OtherContainer>& other) noexcept{
            if constexpr(std::is_same_v<Container, ResMultiArrayContainer>){
                chain -= other.get_chain();
                symmetric_wrtR -= other.get_symmetric_wrtR();
                symmetric_wrtr -= other.get_symmetric_wrtr();
            } else {
                static_assert(std::is_same_v<Container, ResMultiArrayContainer>);
            }
            return *this;
        }

    private:
        Container chain;
        Container symmetric_wrtR;
        Container symmetric_wrtr;
        
        bool pt_denom_set;
        double pt_denom;
    };

    using CARes4Result_MultiArray = CARes4Result<ResMultiArrayContainer>;
    using CARes4Result_Vector = CARes4Result<ResVectorContainer<unsigned>>;
    using CARes4Result_Unbinned = CARes4Result<ResVectorContainer<double>>;

    template <class OtherContainer>
    inline CARes4Result_MultiArray operator-(CARes4Result_MultiArray a, const CARes4Result<OtherContainer>& b) noexcept{
        CARes4Result_MultiArray result(a);
        result -= b;
        return result;
    }

    template <typename Container>
    inline CARes4Result_MultiArray operator-(const CARes4Result_Vector& a, const CARes4Result<Container> b) noexcept{
        return CARes4Result_MultiArray(a) - b;
    }
};

#endif
