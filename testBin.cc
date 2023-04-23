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
#include <chrono>

using namespace std::chrono;

int main(){
    unsigned N=20u;
    unsigned order=3u;
    unsigned Njet=10u;

    using namespace boost::histogram;
    auto dRaxis = axis::variable<double>{0.01, 0.10, 0.30, 0.50, 1.0};
    auto EEChist = make_histogram_with(std::vector<double>(), dRaxis);
    auto COVhist = make_histogram_with(std::vector<double>(), dRaxis, dRaxis);
    auto TRANShist = make_histogram_with(std::vector<double>(), dRaxis, dRaxis);

    auto EEChist_o = make_histogram_with(std::vector<double>(), dRaxis);
    auto COVhist_o = make_histogram_with(std::vector<double>(), dRaxis, dRaxis);

    size_t jettime = 0;
    size_t EECtime = 0;
    size_t covtime = 0;
    size_t hittime = 0;

    for(unsigned i=0; i<Njet; ++i){
        auto t0 = high_resolution_clock::now();
        printf("%u\n", i);
        auto j = std::make_shared<jet>();
        auto j_o = std::make_shared<jet>();

        gausJet(N, *j_o);
        auto ptrans = std::make_shared<arma::mat>(genJet(*j_o, *j, 
                        0.15, 0.01, 0.01,
                        0.00, 0.00, 0.00, 0.00, 0.2));
        auto t1 = high_resolution_clock::now();

        printf("before EEC\n");
        projectedEEC EEC = doProjected(j, order, ptrans, j_o);
        printf("between EECs\n");
        projectedEEC EEC_o = doProjected(j_o, order);
        printf("after EEC\n");

        auto t2 = high_resolution_clock::now();

        arma::mat cov = (*EEC.cov) * arma::trans(*EEC.cov);
        arma::mat cov_o = (*EEC_o.cov) * arma::trans(*EEC_o.cov);

        auto t3 = high_resolution_clock::now();

        fill1d(EEChist, *EEC.dRs, *EEC.wts);
        fill2d(COVhist, *EEC.dRs, *EEC.dRs, cov);
        fill2d(TRANShist, *EEC.dRs_o, *EEC.dRs, *EEC.transfer);


        fill1d(EEChist_o, *EEC_o.dRs, *EEC_o.wts);
        fill2d(COVhist_o, *EEC_o.dRs, *EEC.dRs_o, cov_o);

        auto t4 = high_resolution_clock::now();

        jettime += duration_cast<microseconds>(t1-t0).count();
        EECtime += duration_cast<microseconds>(t2-t1).count();
        covtime += duration_cast<microseconds>(t3-t2).count();
        hittime += duration_cast<microseconds>(t4-t3).count();
    }
    printf("JET TIME %lu\nEEC TIME %lu\nCOV TIME %lu\nHST TIME %lu\n\n", jettime, EECtime, covtime, hittime);

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
