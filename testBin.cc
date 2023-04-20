#include <stdio.h>
#include "eec.h"
#include "toyjets/gen.h"
#include "toyjets/gaus.h"
#include "bin.h"
#include <armadillo>
#include <boost/histogram.hpp>
#include <boost/histogram/axis.hpp>
#include <boost/histogram/ostream.hpp>
#include <boost/histogram/indexed.hpp>

int main(){
    unsigned N=50u;
    unsigned order=2u;
    unsigned Njet=100u;

    using namespace boost::histogram;
    auto dRaxis = axis::variable<double>{0.01, 0.10, 0.30, 0.50, 1.0};
    auto EEChist = make_histogram_with(std::vector<double>(), dRaxis);
    auto COVhist = make_histogram_with(std::vector<double>(), dRaxis, dRaxis);
    auto TRANShist = make_histogram_with(std::vector<double>(), dRaxis, dRaxis);

    auto EEChist_o = make_histogram_with(std::vector<double>(), dRaxis);
    auto COVhist_o = make_histogram_with(std::vector<double>(), dRaxis, dRaxis);

    for(unsigned i=0; i<Njet; ++i){
        printf("%u\n", i);
        auto j = std::make_shared<jet>();
        auto j_o = std::make_shared<jet>();

        gausJet(N, *j_o);
        auto ptrans = std::make_shared<arma::fmat>(genJet(*j_o, *j, 
                        0.15, 0.05, 0.05,
                        0.20, 0.80, 0.00, 0.00, 0.2));

        projectedEEC EEC = doProjected(j, order, ptrans, j_o);
        projectedEEC EEC_o = doProjected(j_o, order);

        arma::fmat cov = (*EEC.cov) * arma::trans(*EEC.cov);

        fill1d(EEChist, *EEC.dRs, *EEC.wts);
        fill2d(COVhist, *EEC.dRs, *EEC.dRs, cov);
        fill2d(TRANShist, *EEC.dRs_o, *EEC.dRs, *EEC.transfer);

        arma::fmat cov_o = (*EEC_o.cov) * arma::trans(*EEC_o.cov);

        fill1d(EEChist_o, *EEC_o.dRs, *EEC_o.wts);
        fill2d(COVhist_o, *EEC_o.dRs, *EEC.dRs_o, cov_o);
    }
    printf("gen EEC\n");
    std::cout << EEChist;
    printf("reco EEC\n");
    std::cout << EEChist_o;
    /* printf("\nbinned COV\n"); */
    /* std::cout << COVhist; */
    /* printf("\nbinned transfer\n"); */
    /* std::cout << TRANShist; */

    arma::fvec EECvec(6, arma::fill::zeros);
    arma::fmat COVmat(6, 6, arma::fill::zeros);
    arma::fmat TRANSmat(6, 6, arma::fill::zeros);

    arma::fvec EECvec_o(6, arma::fill::zeros);
    arma::fmat COVmat_o(6, 6, arma::fill::zeros);

    for(auto &&x : indexed(EEChist, coverage::all)){
        EECvec(x.index(0)+1) = *x;
    }
    for(auto &&x : indexed(COVhist, coverage::all)){
        COVmat(x.index(0)+1, x.index(1)+1) = *x;
    }
    for(auto &&x : indexed(TRANShist, coverage::all)){
        TRANSmat(x.index(0)+1, x.index(1)+1) = *x;
    }

    for(auto &&x : indexed(EEChist_o, coverage::all)){
        EECvec_o(x.index(0)+1) = *x;
    }
    for(auto &&x : indexed(COVhist_o, coverage::all)){
        COVmat_o(x.index(0)+1, x.index(1)+1) = *x;
    }

    printf("GEN\n");
    std::cout << arma::trans(EECvec);
    printf("RECO\n");
    std::cout << arma::trans(EECvec_o);
    printf("T column sum\n");
    std::cout << arma::trans(arma::sum(TRANSmat, 1));
    printf("T row sum\n");
    std::cout << arma::sum(TRANSmat, 0);
    arma::frowvec norm = arma::trans(EECvec);
    norm.replace(0, 1);
    TRANSmat.each_row() /= arma::trans(EECvec);
    printf("normalized T\n");
    std::cout << TRANSmat;
    printf("T * GEN\n");
    std::cout << arma::trans(TRANSmat * EECvec);
}
