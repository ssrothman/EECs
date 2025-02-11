#ifndef SROTHMAN_EECS_Res4TransferResult_H
#define SROTHMAN_EECS_Res4TransferResult_H

#include "usings.h"
#include "Res4Axes.h"
#include "Res4Calculator.h"

#include "SRothman/SimonTools/src/histutil.h"


namespace EEC{

    class Res4TransferMultiArrayContainer{
    public:
        using data_t = multi_array<double, 6>;

        Res4TransferMultiArrayContainer(
                const size_t nR_reco, 
                const size_t nr_reco, 
                const size_t nc_reco,
                const size_t nR_gen, 
                const size_t nr_gen, 
                const size_t nc_gen) noexcept :
            nR_reco(nR_reco), nr_reco(nr_reco), nc_reco(nc_reco),
            nR_gen(nR_gen), nr_gen(nr_gen), nc_gen(nc_gen) {

            data.resize(boost::extents[nR_reco][nr_reco][nc_reco][nR_gen][nr_gen][nc_gen]);
            std::fill(data.data(), data.data() + data.num_elements(), 0);
        }

        Res4TransferMultiArrayContainer() noexcept :
            Res4TransferMultiArrayContainer(0, 0, 0, 0, 0, 0) {}

        void fill(unsigned iR_reco, unsigned ir_reco, unsigned ic_reco,
                  unsigned iR_gen, unsigned ir_gen, unsigned ic_gen,
                  double wt) noexcept {
            data[iR_reco][ir_reco][ic_reco][iR_gen][ir_gen][ic_gen] += wt;
        }

        const data_t& get_data() const noexcept {
            return data;
        }

        Res4TransferMultiArrayContainer& operator+=(const Res4TransferMultiArrayContainer& other) noexcept {
            assert(data.num_elements() == other.data.num_elements());
            std::transform(data.data(), data.data() + data.num_elements(),
                           other.data.data(), data.data(),
                           std::plus<double>());
            return *this;
        }

        double total_weight() const noexcept {
            return std::accumulate(data.data(), data.data() + data.num_elements(), 0.0);
        }
        const size_t nR_reco, nr_reco, nc_reco;
        const size_t nR_gen, nr_gen, nc_gen;
    private:
        data_t data;
    };

    class Res4TransferVectorContainer{
    public:
        struct entry{
            unsigned iR_reco;
            unsigned ir_reco;
            unsigned ic_reco;
            unsigned iR_gen;
            unsigned ir_gen;
            unsigned ic_gen;
            double wt;
            entry(unsigned iR_reco, unsigned ir_reco, unsigned ic_reco,
                  unsigned iR_gen, unsigned ir_gen, unsigned ic_gen,
                  double wt) noexcept  :
                iR_reco(iR_reco), ir_reco(ir_reco), ic_reco(ic_reco),
                iR_gen(iR_gen), ir_gen(ir_gen), ic_gen(ic_gen), wt(wt) {}
        };
        using data_t = std::vector<entry>;

        Res4TransferVectorContainer(
                const size_t nR_reco, 
                const size_t nr_reco, 
                const size_t nc_reco,
                const size_t nR_gen, 
                const size_t nr_gen, 
                const size_t nc_gen) noexcept :
            nR_reco(nR_reco), nr_reco(nr_reco), nc_reco(nc_reco),
            nR_gen(nR_gen), nr_gen(nr_gen), nc_gen(nc_gen) {

            data.reserve(nR_reco);
        }

        Res4TransferVectorContainer() noexcept :
            Res4TransferVectorContainer(0, 0, 0, 0, 0, 0) {}

        void fill(unsigned iR_reco, unsigned ir_reco, unsigned ic_reco,
                  unsigned iR_gen, unsigned ir_gen, unsigned ic_gen,
                  double wt) noexcept {
            data.emplace_back(iR_reco, ir_reco, ic_reco, iR_gen, ir_gen, ic_gen, wt);
        }

        const data_t& get_data() const noexcept {
            return data;
        }

        Res4TransferVectorContainer& operator+=(const Res4TransferVectorContainer& other) noexcept {
            data.insert(data.end(), other.data.begin(), other.data.end());
            return *this;
        }

        double total_weight() const noexcept {
            return std::accumulate(data.begin(), data.end(), 0.0,
                                   [](double sum, const entry& e) { return sum + e.wt; });
        }
        const size_t nR_reco, nr_reco, nc_reco;
        const size_t nR_gen, nr_gen, nc_gen;
    private:
        data_t data;
    };

    template <class TransferContainer, class BasicContainer>
    class Res4TransferResult{
    public:
        Res4Result<BasicContainer> unmatched_reco;
        Res4Result<BasicContainer> unmatched_gen;

        Res4TransferResult() noexcept  :
            unmatched_reco(), unmatched_gen(),
            dipole(), tee(), triangle(),
            pt_denom_set(false),
            pt_denom_reco(-1), pt_denom_gen(-1) {}

        Res4TransferResult(
                const size_t nR_reco,
                const size_t nr_dipole_reco,
                const size_t nc_dipole_reco,
                const size_t nr_tee_reco,
                const size_t nc_tee_reco,
                const size_t nr_triangle_reco,
                const size_t nc_triangle_reco,
                const size_t nR_gen,
                const size_t nr_dipole_gen,
                const size_t nc_dipole_gen,
                const size_t nr_tee_gen,
                const size_t nc_tee_gen,
                const size_t nr_triangle_gen,
                const size_t nc_triangle_gen) noexcept  :
            unmatched_reco(nR_reco, nr_dipole_reco, nc_dipole_reco,
                       nr_tee_reco, nc_tee_reco,
                       nr_triangle_reco, nc_triangle_reco),
            unmatched_gen(nR_gen, nr_dipole_gen, nc_dipole_gen,
                       nr_tee_gen, nc_tee_gen,
                       nr_triangle_gen, nc_triangle_gen),
            dipole(nR_reco, nr_dipole_reco, nc_dipole_reco, 
                    nR_gen, nr_dipole_gen, nc_dipole_gen),
            tee(nR_reco, nr_tee_reco, nc_tee_reco,
                nR_gen, nr_tee_gen, nc_tee_gen),
            triangle(nR_reco, nr_triangle_reco, nc_triangle_reco,
                     nR_gen, nr_triangle_gen, nc_triangle_gen),
            pt_denom_set(false),
            pt_denom_reco(-1), pt_denom_gen(-1) {}

        Res4TransferResult(
                const Res4Axes& axes_reco,
                const Res4Axes& axes_gen) noexcept  :
            Res4TransferResult(
                    simon::AXextent(axes_reco.R),
                    simon::AXextent(axes_reco.r_dipole),
                    simon::AXextent(axes_reco.c_dipole),
                    simon::AXextent(axes_reco.r_tee),
                    simon::AXextent(axes_reco.c_tee),
                    simon::AXextent(axes_reco.r_triangle),
                    simon::AXextent(axes_reco.c_triangle),
                    simon::AXextent(axes_gen.R),
                    simon::AXextent(axes_gen.r_dipole),
                    simon::AXextent(axes_gen.c_dipole),
                    simon::AXextent(axes_gen.r_tee),
                    simon::AXextent(axes_gen.c_tee),
                    simon::AXextent(axes_gen.r_triangle),
                    simon::AXextent(axes_gen.c_triangle)) {}

        Res4TransferResult(const Res4TransferCalculator& calculator) noexcept:
            Res4TransferResult(calculator.get_axes_reco(), calculator.get_axes_gen()) {}
        
        void fill_dipole(unsigned iR_reco, unsigned ir_reco, unsigned ic_reco,
                         unsigned iR_gen, unsigned ir_gen, unsigned ic_gen,
                         double wt) noexcept {
            dipole.fill(iR_reco, ir_reco, ic_reco, iR_gen, ir_gen, ic_gen, wt);
        }

        void fill_tee(unsigned iR_reco, unsigned ir_reco, unsigned ic_reco,
                      unsigned iR_gen, unsigned ir_gen, unsigned ic_gen,
                      double wt) noexcept {
            tee.fill(iR_reco, ir_reco, ic_reco, iR_gen, ir_gen, ic_gen, wt);
        }

        void fill_triangle(unsigned iR_reco, unsigned ir_reco, unsigned ic_reco,
                           unsigned iR_gen, unsigned ir_gen, unsigned ic_gen,
                           double wt) noexcept {
            triangle.fill(iR_reco, ir_reco, ic_reco, iR_gen, ir_gen, ic_gen, wt);
        }

        void set_pt_denom(double pt_denom_reco_, double pt_denom_gen_) {
            if (pt_denom_set){
                throw std::runtime_error("pt_denom already set");
            } else {
                pt_denom_reco = pt_denom_reco_;
                pt_denom_gen = pt_denom_gen_;
                pt_denom_set = true;
            }
        }

        double get_pt_denom_reco() const {
            if (!pt_denom_set){
                throw std::runtime_error("pt_denom not set");
            } else {
                return pt_denom_reco;
            }
        }

        double get_pt_denom_gen() const {
            if (!pt_denom_set){
                throw std::runtime_error("pt_denom not set");
            } else {
                return pt_denom_gen;
            }
        }

        const TransferContainer& get_dipole() const  noexcept {
            return dipole;
        }

        const TransferContainer& get_tee() const  noexcept {
            return tee;
        }

        const TransferContainer& get_triangle() const  noexcept {
            return triangle;
        }

        Res4TransferResult<TransferContainer, BasicContainer>& operator+=(const Res4TransferResult<TransferContainer, BasicContainer>& other) noexcept {
            unmatched_reco += other.unmatched_reco;
            unmatched_gen += other.unmatched_gen;
            dipole += other.dipole;
            tee += other.tee;
            triangle += other.triangle;
            return *this;
        }

        double total_dipole_weight() const noexcept {
            return dipole.total_weight();
        }

        double total_tee_weight() const noexcept {
            return tee.total_weight();
        }

        double total_triangle_weight() const noexcept {
            return triangle.total_weight();
        }
    private:
        TransferContainer dipole;
        TransferContainer tee;
        TransferContainer triangle;

        bool pt_denom_set;
        double pt_denom_reco, pt_denom_gen;
    };

    using Res4TransferResult_MultiArray_MultiArray = Res4TransferResult<Res4TransferMultiArrayContainer, Res4MultiArrayContainer>;
    using Res4TransferResult_MultiArray_Vector = Res4TransferResult<Res4TransferMultiArrayContainer, Res4VectorContainer>;

    using Res4TransferResult_Vector_MultiArray = Res4TransferResult<Res4TransferVectorContainer, Res4MultiArrayContainer>;
    using Res4TransferResult_Vector_Vector = Res4TransferResult<Res4TransferVectorContainer, Res4VectorContainer>;
};

#endif
