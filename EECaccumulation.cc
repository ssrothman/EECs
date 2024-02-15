#include "EECaccumulation.h"
#include "SRothman/SimonTools/src/util.h"
#include "SRothman/SimonTools/src/deltaR.h"

/*
 * Generic delta-R handler
 * Parent class for all actual cases
 * Includes generic logic for finding max dR of a set 
 * And also for finding the full sorted list of dRs 
 */
/*
 * Parent class for delta-R handlers 
 * that need to track floating-point delta-Rs
 * ie this is needed for some resolved handlers 
 * which need to compute angles at evaluation time
 *
 * This is done because precomputing all the possible angles 
 * would be too memory-intensive
 */

static uvec fullOrd(const uvec& ord, const uvec& comp) {
    uvec ord_full;
    for(unsigned i=0; i<ord.size(); ++i){
        ord_full.insert(ord_full.end(), comp[i], ord[i]);
    }
    return ord_full;
}

template <typename val_t>
class DRType : public EECindexer{
    public:
        typedef std::vector<val_t> valvec;

        DRType(const jet& J) :
            dRs_({J.nPart, J.nPart}, 0) {}
        virtual ~DRType() {};

    protected:
        ArbitraryMatrix<val_t> dRs_;

        val_t maxDR(const uvec& ord) const{
            if(uniform(ord)){
                return 0;
            }

            unsigned M = ord.size();
            uvec ord2 = ord0_nodiag(2);

            val_t maxDR=0;
            do {
                uvec ord_test(2);
                ord_test[0] = ord[ord2[0]];
                ord_test[1] = ord[ord2[1]];

                if(uniform(ord_test)){
                    continue;
                }

                val_t next = dRs_.at(ord_test);
                if(next > maxDR){
                    maxDR = next;
                }
            } while(iterate_nodiag(2u, ord2, M));

            return maxDR;
        }

        valvec sortedDRs(const uvec& ord, bool sort=true) const {
            return sortedDRs(ord, uvec(ord.size(), 1u), sort);
        }

        valvec sortedDRs(const uvec& ord,
                         const uvec& comp, bool sort=true) const {
            uvec ord_full = fullOrd(ord, comp);
            
            unsigned M = ord_full.size();
            
            if(uniform(ord)){
                return valvec(M, 0);
            }

            valvec dRs;
            uvec ord2 = ord0_nodiag(2);
            do {
                uvec ord_test(2);
                ord_test[0] = ord_full[ord2[0]];
                ord_test[1] = ord_full[ord2[1]];

                val_t next = dRs_.at(ord_test);
                dRs.push_back(next);
            } while(iterate_nodiag(2u, ord2, M));

            if(sort){
                std::sort(dRs.begin(), dRs.end());
            }

            return dRs;
        }
};



class DoubleDR : public DRType<double>{
    public:
        DoubleDR(const jet& J) : 
            DRType(J) {
                fillDRs(J);
            }
        virtual ~DoubleDR(){};

    protected:
        void fillDRs(const jet& j){
            auto ord = dRs_.ord0();
            do {
                if(uniform(ord)){
                    dRs_.at(ord) = 0;
                } else {
                    const auto& p1 = j.particles.at(ord[0]);
                    const auto& p2 = j.particles.at(ord[1]);
                    dRs_.at(ord) = dR(p1, p2);
                }
            } while(dRs_.iterate(ord));
        }
};

/*
 * Parent class for delta-R handlers
 * that only care about dR index
 */
class IndexedDR : public DRType<unsigned>{
    public:
        IndexedDR(const jet& J,
                  const axis_t& axis) : 
            DRType(J) {
                fillDRs(J, axis);
            }
        ~IndexedDR(){};

    protected:
        void fillDRs(const jet& j, const axis_t& axis){
            auto ord = dRs_.ord0();
            do {
                if(uniform(ord)){
                    dRs_.at(ord) = 0;
                } else {
                    const auto& p1 = j.particles.at(ord[0]);
                    const auto& p2 = j.particles.at(ord[1]);
                    double deltaR = dR(p1, p2);
                    dRs_.at(ord) = unsigned(axis.index(deltaR) + 1);
                }
            } while(dRs_.iterate(ord));
        }
};

class IndexedMaxDR : public IndexedDR{
    public:
        IndexedMaxDR(const jet& J,
                  const axis_t& axis) : 
            IndexedDR(J, axis) {}
        ~IndexedMaxDR(){};

        uvec getIndex(const uvec& ord,
                      const uvec& comp) const {
            return {maxDR(ord)};
        }
};

class IndexedSortedDRs : public IndexedDR{
    public:
        IndexedSortedDRs(const jet& J,
                         const axis_t& axis,
                         const unsigned Nax) : 
            IndexedDR(J, axis),
            Nax_(Nax) {}
        ~IndexedSortedDRs(){};

        uvec getIndex(const uvec& ord,
                      const uvec& comp) const {
            auto sorted = sortedDRs(ord, comp);
            if(sorted.size() < Nax_){
                throw std::invalid_argument("IndexedSortedDRs::getIndex "
                        "requires Nax_ <= ord.size()");
            }
            return uvec(sorted.begin(), sorted.begin() + Nax_);
        }
    protected:
        unsigned Nax_;
};

//this is the coordinate system described by Ian
//and used in https://arxiv.org/pdf/2205.02857.pdf [non-gaussianities]
//and used in https://arxiv.org/pdf/2201.07800.pdf [Jesse open data]
class ThirdOrderCoords : public DoubleDR{
    public:
        ThirdOrderCoords(const jet& J,
                         const axisvec& axes) :
            DoubleDR(J),
            axes_(axes) {
                if(axes_.size() != 3){
                    throw std::invalid_argument(
                            "ThirdOrderCoords must be "
                            "constructed with 3 axes "
                            "(RL, xi, phi)");
                }
        }
        ~ThirdOrderCoords() {};

        uvec getIndex(const uvec& ord,
                      const uvec& comp) const {
            valvec sorted = sortedDRs(ord, comp);

            if(sorted.size() < 3){
                throw std::invalid_argument(
                        "ThirdOrderCoords::getIndex "
                        "requires ord.size() >= 3");
            }

            double RL = sorted[0];
            double RM = sorted[1];
            double RS = sorted[2];

            double xi = RS/RM;
            double phi = std::asin(std::sqrt(1 - square((RL-RM)/RS)));

            return {unsigned(axes_[0].index(RL) + 1),
                    unsigned(axes_[1].index(xi) + 1),
                    unsigned(axes_[2].index(phi) + 1)};
        }
    protected:
        axisvec axes_;
};

static double angle_from_distances(double a, double b, double c){
    if(a==0 || b==0){
        return 0; //answer is undefined
    } else if(c==0){ 
        return 0; //answer is well-defined but weird
    }
    return std::acos((square(a) + square(b) - square(c)) / (2*a*b));
}

using edge = std::pair<unsigned,unsigned>;
using dist = std::pair<double,edge>;
static bool handle_overlap(dist d1, dist d2,
                           int& i, int& j, int& k,
                           double& RL, double& RM){
    if(d1.second.first == d2.second.first){
        i = d1.second.first;
        j = d1.second.second;
        k = d2.second.second;
        RL = d1.first;
        RM = d2.first;
        return true;
    }else if(d1.second.first == d2.second.second){
        i = d1.second.first;
        j = d1.second.second;
        k = d2.second.first;
        RL = d1.first;
        RM = d2.first;
        return true;
    }else if(d1.second.second == d2.second.first){
        i = d1.second.second;
        j = d1.second.first;
        k = d2.second.second;
        RL = d1.first;
        RM = d2.first;
        return true;
    }else if(d1.second.second == d2.second.second){
        i = d1.second.second;
        j = d1.second.first;
        k = d2.second.first;
        RL = d1.first;
        RM = d2.first;
        return true;
    }else{
        return false;
    }
}

static void try_l(const dist& d, 
                  const int i, const int j, const int k, int& l, 
                  const double RL, const double RM, 
                  double& RS, double& R2, double& R2L){
    //we know we've worked out (i,j) and (i,k)
    //remaining are (i, l), (j,k), (j,l), or (k,l)
    //just test to see which one we've got

    if(int(d.second.first) == i){ //(i, l)
        R2 = d.first;
        l = int(d.second.second);
    }else if(int(d.second.second)==i){//(l, i)
        R2 = d.first;
        l = int(d.second.first);
    }else if(int(d.second.first)==j && int(d.second.second)==k){//(j,k)
        RS = d.first;
    }else if(int(d.second.second)==j && int(d.second.first)==k){//(k,j)
        RS = d.first;
    }else if(int(d.second.first)==j){//(j, l)
        l = int(d.second.second);
        R2L = d.first;
    }else if(int(d.second.second)==j){//(l, j)
        l = int(d.second.first);
        R2L = d.first;
    }else if(int(d.second.first)==k){//(k, l)
        l = int(d.second.second);
    }else if(int(d.second.second)==k){//(l,k)
        l = int(d.second.first);
    }else{
        throw std::logic_error("try_l should always succeed in one of the cases");
    }
}

static void find_ijkl(const std::vector<dist>& distances,
                      int& i, int& j, int& k, int& l,
                      double& RL, double& RM, double& RS, double& R2,
                      double& R2L){
    if(handle_overlap(distances[0], distances[1], i, j, k, RL, RM)){
        try_l(distances[2], i, j, k, l, RL, RM, RS, R2, R2L);
    }else if(handle_overlap(distances[0], distances[2], i,j,k,RL,RM)){
        try_l(distances[1], i, j, k, l, RL, RM, RS, R2, R2L);
    } else{
        throw std::logic_error("There should be overlap between either the first and second or first and third distances");
    }
    for(int idist=3; idist<6; ++idist){
        try_l(distances[idist], i, j, k, l, RL, RM, RS, R2, R2L);
    }
}

/*
 * Proposal from Ian was to fix a triangle
 * and then move the point inside of it
 *
 * we form a triangle by 
 * 1. insisting that the longest side is on the triangle
 * 2. selecting the other point in the triangle to maximize
 *      the length of the second-longest side
 *      on the triangle
 * 
 * Note that the second-longest side of the triangle
 * might not be the second-longest overall
 *
 * To identify triangle shapes we consider write the
 * side lengths of the triangle as RL, RM, RS (sorted)
 * and place restrictions on RM/RL and RS/RL, with some tolerance
 *
 * We identify the distance R2 as the distance from the 
 * vertex between the two longest sides of the triangle and the
 * fourth point (not on the triangle)
 *
 * The binning variables are then RL, R2/RL, 
 *      and phi=the angle between RL and R2
 */


class FourthOrderFixedTriangle : public DoubleDR{
    public:

        FourthOrderFixedTriangle(const jet& J,
                      const axisvec& axes,
                      const struct trianglespec& spec) :
            DoubleDR(J),
            axes_(axes),
            spec_(spec) {}

        ~FourthOrderFixedTriangle() {};

        uvec getIndex(const uvec& ord,
                      const uvec& comp) const {
            uvec ord_full = fullOrd(ord, comp);

            if(ord_full.size() < 4){
                throw std::invalid_argument(
                        "FourthOrderFixedTriangle::getIndex "
                        "requires ord.size() >= 4");
            }

            /* 
             * Variable definitions:
             * RL = longest side of triangle
             * RM = second-longest side of triangle
             * RS = shortest side of triangle
             * R2 = distance from vertex between RL and RM to 4th point
             *
             * i, j, k, l are the four points, ordered such that
             * RL = D(i, j)
             * RM = D(i, k)
             * RS = D(j, k)
             * R2 = D(i, l)
             * R2L = D(j, l)
             */
            std::vector<dist> distances;
            uvec ord_pair = ord0_nodiag(2);
            do{
                uvec ord_test(2);
                ord_test[0] = ord_full[ord_pair[0]];
                ord_test[1] = ord_full[ord_pair[1]];
                edge e(ord_test[0], ord_test[1]);
                distances.emplace_back(dRs_.at(ord_test), e);
            } while(iterate_nodiag(2u, ord_pair, 4u));

            std::sort(distances.begin(), distances.end(),
                    [](const dist& a,
                       const dist& b){
                        return a.first > b.first;
                    }
            );


            int i=-1, j=-1, k=-1, l=-1;
            double RL=-1, RM=-1, RS=-1, R2=-1, R2L = -1;
            find_ijkl(distances, i, j, k, l, RL, RM, RS, R2, R2L);

            printf("Found triangle:\n");
            printf("\ti=%d, j=%d, k=%d, l=%d\n", i, j, k, l);
            printf("\tRL=%f, RM=%f, RS=%f, R2=%f, R2L=%f\n", RL, RM, RS, R2, R2L);

            if(std::abs(RM/RL - spec_.RMoRL) > spec_.tol ||
                    std::abs(RS/RL - spec_.RSoRL) > spec_.tol){
                return {};
            } else {
                double phi = angle_from_distances(RL, R2, R2L);
                return {unsigned(axes_[0].index(RL) + 1),
                        unsigned(axes_[1].index(R2/RL) + 1),
                        unsigned(axes_[2].index(phi) + 1)};
            }
            return {};
        }

    protected:
        axisvec axes_;
        struct trianglespec spec_;
};


EECweightAccumulator::EECweightAccumulator(unsigned Naccu, 
                                           const axisvec& axes):
    indexer_(nullptr),
    axes_(axes){
        fillSizes(axes);
        for (unsigned i=0; i<Naccu; ++i){
            accus_.emplace_back(sizes_);
        }
}

EECweightAccumulator::EECweightAccumulator(unsigned Naccu,
                                           const axis_t& axis): 
    EECweightAccumulator(Naccu, axisvec({axis})) {} 

void EECweightAccumulator::setupMaxDRIndexer(const jet& j){
    if(axes_.size() != 1){
        throw std::logic_error(
                "maxDR indexer only works with 1D binning");
    }

    indexer_ = std::make_unique<IndexedMaxDR>(j, axes_[0]);
}

void EECweightAccumulator::setupSortedDRIndexer(const jet& j){
    setupSortedDRIndexer(j, axes_.size());
}

void EECweightAccumulator::setupSortedDRIndexer(const jet& j, 
                                                const unsigned Nax){
    if(axes_.size() != 1){
        throw std::logic_error(
                "Sorted DR indexer clones a 1D binning "
                "and is not compatible with different "
                "binnings in different axes");
    }
    indexer_ = std::make_unique<IndexedSortedDRs>(
            j, axes_[0], Nax);
}

void EECweightAccumulator::setupThirdOrderIndexer(const jet& j){
    if(axes_.size() != 3){
        throw std::logic_error(
            "Third order indexer only works with 3D binning");
    }
    indexer_ = std::make_unique<ThirdOrderCoords>(j, axes_);
}

void EECweightAccumulator::setupFourthOrderIndexer(const jet& j,
        const struct trianglespec& spec){
    if(axes_.size() != 3){
        throw std::logic_error(
            "Fourth order indexer only works with 3D binning");
    }
    indexer_ = std::make_unique<FourthOrderFixedTriangle>(
            j, axes_, spec);
}
        
void EECweightAccumulator::setupIndexer(const edm::ParameterSet& conf,
                                        const jet& j){
    auto type = conf.getParameter<std::string>("type");
    if(type == "maxDR"){
        setupMaxDRIndexer(j);
    } else if(type == "sortedDR"){
        auto Nax = conf.getParameter<unsigned>("Nax");
        setupSortedDRIndexer(j, Nax);
    } else if(type == "thirdOrder"){
        setupThirdOrderIndexer(j);
    } else if(type == "fourthOrder"){
        trianglespec spec(conf);
        setupFourthOrderIndexer(j, spec);
    } else {
        throw std::logic_error(
                "Invalid indexer type: " + type);
    }
}

void EECweightAccumulator::accumulate(const uvec& ord, const uvec& comp, 
        const std::vector<double>& weights,
        unsigned offset,
        const uvec* const input_idx,
        uvec* return_idx){

    if (input_idx){
        accumulate(*input_idx, weights, offset);
        if(return_idx){
            *return_idx = *input_idx;
        }
    } else {
        uvec idx = indexer_->getIndex(ord, comp);
        accumulate(idx, weights, offset);
        if(return_idx){
            *return_idx = idx;
        }
    }
}

void EECweightAccumulator::accumulate(const uvec& idx,
        const std::vector<double>& wts, unsigned offset){
    for(unsigned i=offset; i<wts.size(); ++i){
        accus_.at(i).at(idx)(wts[i]);
    }
}

void EECweightAccumulator::accumulate(const uvec& ord,
        const uvec& comp, double wt, unsigned iAcc){

    uvec idx = indexer_->getIndex(ord, comp);
    accus_.at(iAcc).at(idx)(wt);
}

double EECweightAccumulator::get(unsigned iAcc, const uvec& ord) const{
    return boost::accumulators::extract_result<stats_tag>(
            accus_.at(iAcc).at(ord));
}

ArbitraryMatrix<double> EECweightAccumulator::data(unsigned iAcc) const{
    ArbitraryMatrix<double> result(sizes_);
    uvec ord0 = result.ord0();
    do {
        result.at(ord0) = get(iAcc, ord0);
    } while(result.iterate(ord0));
    return result;
}

void EECweightAccumulator::fillSizes(const axisvec& axes){
    sizes_.resize(axes.size());
    for(size_t i = 0; i < axes.size(); ++i){
        sizes_[i] = boost::histogram::axis::traits::extent(axes[i]);
    }
}

EECtransferAccumulator::EECtransferAccumulator(unsigned Naccu, 
                                               const axisvec& axes):
    indexerReco_(nullptr),
    indexerGen_(nullptr),
    axes_(axes){
        fillSizes(axes);
        for(const auto& s : sizes_){
            sizes_.push_back(s);
        }
        for (unsigned i=0; i<Naccu; ++i){
            accus_.emplace_back(sizes_);
        }
}

EECtransferAccumulator::EECtransferAccumulator(unsigned Naccu,
                                               const axis_t& axis): 
    EECtransferAccumulator(Naccu, axisvec({axis})) {}

void EECtransferAccumulator::setupMaxDRIndexers(const jet& jReco, 
                                                const jet& jGen){
    if(axes_.size() != 1){
        throw std::logic_error(
                "maxDR indexer only works with 1D binning");
    }

    indexerReco_ = std::make_unique<IndexedMaxDR>(jReco, axes_[0]);
    indexerGen_ = std::make_unique<IndexedMaxDR>(jGen, axes_[0]);
}

void EECtransferAccumulator::setupSortedDRIndexers(const jet& jReco,
                                                const jet& jGen){
    setupSortedDRIndexers(jReco, jGen, axes_.size());
}

void EECtransferAccumulator::setupSortedDRIndexers(const jet& jReco,
                                                const jet& jGen,
                                                const unsigned Nax){
    if(axes_.size() != 1){
        throw std::logic_error(
                "Sorted DR indexer clones a 1D binning "
                "and is not compatible with different "
                "binnings in different axes");
    }
    indexerReco_ = std::make_unique<IndexedSortedDRs>(
            jReco, axes_[0], Nax);
    indexerGen_ = std::make_unique<IndexedSortedDRs>(
            jGen, axes_[0], Nax);
}

void EECtransferAccumulator::setupThirdOrderIndexers(const jet& jReco,
                                                  const jet& jGen){
    if(axes_.size() != 3){
        throw std::logic_error(
            "Third order indexer only works with 3D binning");
    }
    indexerReco_ = std::make_unique<ThirdOrderCoords>(jReco, axes_);
    indexerGen_ = std::make_unique<ThirdOrderCoords>(jGen, axes_);
}

void EECtransferAccumulator::setupFourthOrderIndexers(const jet& jReco,
        const jet& jGen, const struct trianglespec& spec){
    if(axes_.size() != 3){
        throw std::logic_error(
            "Fourth order indexer only works with 3D binning");
    }
    indexerReco_ = std::make_unique<FourthOrderFixedTriangle>(
            jReco, axes_, spec);
    indexerGen_ = std::make_unique<FourthOrderFixedTriangle>(
            jGen, axes_, spec);
}
        
void EECtransferAccumulator::setupIndexers(const edm::ParameterSet& conf,
                                        const jet& jReco, 
                                        const jet& jGen){
    auto type = conf.getParameter<std::string>("type");
    if(type == "maxDR"){
        setupMaxDRIndexers(jReco, jGen);
    } else if(type == "sortedDR"){
        auto Nax = conf.getParameter<unsigned>("Nax");
        setupSortedDRIndexers(jReco, jGen, Nax);
    } else if(type == "thirdOrder"){
        setupThirdOrderIndexers(jReco, jGen);
    } else if(type == "fourthOrder"){
        trianglespec spec(conf);
        setupFourthOrderIndexers(jReco, jGen, spec);
    } else {
        throw std::logic_error(
                "Invalid indexer type: " + type);
    }
}

void EECtransferAccumulator::accumulate(const uvec& ordReco, 
                                        const uvec& ordGen, 
                                        const uvec& compReco,
                                        const uvec& compGen,
                                        const std::vector<double>& weights){
    uvec indexReco = indexerReco_->getIndex(ordReco, compReco);
    uvec indexGen = indexerGen_->getIndex(ordGen, compGen);
    uvec index = indexReco;
    index.insert(index.end(), indexGen.begin(), indexGen.end());
    for(unsigned i=0; i<weights.size(); ++i){
        accus_.at(i).at(index)(weights[i]);
    }
}

void EECtransferAccumulator::accumulate(const uvec& ordGen,
                                        const uvec& compGen,
                                        const adjacency& adj,
                                        const arma::mat& ptrans,
                                        double nextwt,
                                        unsigned iAcc){
    printf("accumulating transfer\n");
    uvec indexGen = indexerGen_->getIndex(ordGen, compGen);

    uvec ord_J1 = fullOrd(ordGen, compGen);
    unsigned order = ord_J1.size();
    uvec nadj = adj.nadj(ord_J1); 
    if (nadj.empty()){
        printf("breaking early\n");
        return; //break early if any have no neighbors
    } 

    uvec ord_iter = ord0_full(order);
    do{
        double tfact = 1.0;
        uvec ord_J2(order);
        for(unsigned i=0; i<order; ++i){
            unsigned p1 = ord_J1[i];
            unsigned p2 = adj.at(p1)[ord_iter[i]];
            ord_J2[i] = p2;

            printf("\tadj: ");
            printOrd(adj.at(p1));
            printf("\n");
            printf("\tord_iter[i] = %u\n", ord_iter[i]);
            printf("\tnadj[i] = %u\n", nadj[i]);
            printf("\tp1: %u, p2: %u\n", p1, p2);
            tfact *= ptrans(p2, p1);
        }
        uvec index = indexerReco_->getIndex(ord_J2);
        index.insert(index.end(), indexGen.begin(), indexGen.end());
        accus_.at(iAcc).at(index)(nextwt * tfact);
    } while(iterate_awkward(nadj, ord_iter));
    printf("done\n");
}

double EECtransferAccumulator::get(unsigned iAcc, const uvec& ord) const{
    return boost::accumulators::extract_result<stats_tag>(
            accus_.at(iAcc).at(ord));
}

ArbitraryMatrix<double> EECtransferAccumulator::data(unsigned iAcc) const{
    ArbitraryMatrix<double> result(sizes_);
    uvec ord0 = result.ord0();
    do {
        result.at(ord0) = get(iAcc, ord0);
    } while(result.iterate(ord0));
    return result;
}

void EECtransferAccumulator::fillSizes(const axisvec& axes){
    sizes_.resize(axes.size());
    for(size_t i = 0; i < axes.size(); ++i){
        sizes_[i] = boost::histogram::axis::traits::extent(axes[i]);
    }
}
