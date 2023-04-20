#ifndef EECs_EEC_H
#define EECs_EEC_H

#include <vector>
#include <memory>
#include <armadillo>

#include "simon_util_cpp/vecND.h"
#include "simon_util_cpp/combinatorics.h"

#include "toyjets/common.h"

//output struct
struct projectedEEC{
    unsigned order;
    unsigned nPart;
    std::unique_ptr<vecND::nodiagvec> dR2s, wts;
    float ptAtZero;
    std::unique_ptr<arma::fmat> cov;
    std::unique_ptr<vecND::nodiagvec> dR2s_o;
    std::unique_ptr<arma::fmat> transfer;

    projectedEEC(unsigned order, unsigned nPart,
                 std::unique_ptr<vecND::nodiagvec>&& dR2s,
                 std::unique_ptr<vecND::nodiagvec>&& wts,
                 float ptAtZero,
                 std::unique_ptr<arma::fmat>&& cov,
                 std::unique_ptr<vecND::nodiagvec>&& dR2s_o,
                 std::unique_ptr<arma::fmat>&& transfer) :
        order(order),
        nPart(nPart),
        dR2s(std::move(dR2s)),
        wts(std::move(wts)),
        ptAtZero(ptAtZero),
        cov(std::move(cov)),
        dR2s_o(std::move(dR2s_o)),
        transfer(std::move(transfer)){}
};

//data to be carried around the EEC calculation
//needs only to be computed once per jet
//and is constant from then on
struct constdata{
    unsigned order;        //correlator order
    unsigned nPart;        //number of particles in jet
    std::unique_ptr<std::vector<float>> Es; //jet constituent momenta 
                           //normalized to jet momentum
    std::unique_ptr<vecND::nodiagvec> dR2s; //deltaR^2 between pairs of particles in jet
    std::unique_ptr<comp_t> comps;          //integer compositions of order
    std::shared_ptr<const arma::fmat> ptrans; //particle-level transfer matrix
    unsigned nPart_o; //nPart in other jet
    std::unique_ptr<vecND::nodiagvec> dR2s_o; //other dR&2 array

    constdata(unsigned order, unsigned nPart, 
              std::unique_ptr<std::vector<float>>&& Es,
              std::unique_ptr<vecND::nodiagvec>&& dR2s,
              std::unique_ptr<comp_t>&& comps):
        order(order),
        nPart(nPart),
        Es(std::move(Es)),
        dR2s(std::move(dR2s)),
        comps(std::move(comps)),
        ptrans(nullptr),
        dR2s_o(nullptr),
        nPart_o(0){}

    constdata(unsigned order, unsigned nPart, 
              std::unique_ptr<std::vector<float>>&& Es,
              std::unique_ptr<vecND::nodiagvec>&& dR2s,
              std::unique_ptr<comp_t>&& comps,
              std::shared_ptr<const arma::fmat> ptrans,
              unsigned nPart_o,
              std::unique_ptr<vecND::nodiagvec>&& dR2s_o):
        order(order),
        nPart(nPart),
        Es(std::move(Es)),
        dR2s(std::move(dR2s)),
        comps(std::move(comps)),
        ptrans(ptrans),
        dR2s_o(std::move(dR2s_o)),
        nPart_o(nPart_o) {}
};

//data that is accumulated over the course of the EEC calculation
struct computedata{
    std::unique_ptr<vecND::nodiagvec> wts;
    std::unique_ptr<arma::fmat> cov;
    std::unique_ptr<arma::fmat> transfer;
    float ptAtZero;

    computedata(std::unique_ptr<vecND::nodiagvec>&& wts,
                std::unique_ptr<arma::fmat>&& cov,
                std::unique_ptr<arma::fmat>&& transfer):
        wts(std::move(wts)),
        cov(std::move(cov)),
        transfer(std::move(transfer)),
        ptAtZero(0) {} 

    computedata(std::unique_ptr<vecND::nodiagvec>&& wts,
                std::unique_ptr<arma::fmat>&& cov):
        wts(std::move(wts)),
        cov(std::move(cov)),
        transfer(nullptr),
        ptAtZero(0) {}
};

//get a computedata with everything appropriately sized
//and filled with zeros
computedata initializeComputedata(const constdata& cd);

//run the full computation
projectedEEC doProjected(const std::shared_ptr<const jet> j,
                         const unsigned order,
                         const std::shared_ptr<const arma::fmat> ptrans=nullptr,
                         const std::shared_ptr<const jet> j_o=nullptr);

//setup the constdata struct from the given jet(s)
constdata getConstdata(const std::shared_ptr<const jet> j,
                       const unsigned order,
                       const std::shared_ptr<const jet> j_o = nullptr,
                       const std::shared_ptr<const arma::fmat> ptrans = nullptr);

//compute the contribution from the 0-DR contact terms
void pointAtZero(const constdata& cd, computedata& ans);

//compute the contribution from M-tuples of particles
void MwayContribution(const constdata& cd,
                      const unsigned M,
                      computedata& ans);

//accumulate weight for a given configuration
void addWt(const constdata& cd,
           const unsigned M,
           const std::vector<unsigned>& ord,
           computedata& ans);

//get largest dR of a configuration
unsigned getMaxDR(const vecND::nodiagvec& dR2s,
                  const std::vector<unsigned> ord);

//populate dR^2 matrix
void getDR2(vecND::nodiagvec& result,
            const std::vector<float>& eta,
            const std::vector<float>& phi);

//normalize momentum vector pts to have unit sum
void normalizePt(const std::vector<float>& pts,
                 std::vector<float>& Eout);

#endif
