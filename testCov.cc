#include <stdio.h>
#include "eec.h"
#include "toyjets/gen.h"
#include "toyjets/gaus.h"

int main(){
    unsigned N=5u;
    unsigned order=3u;
    auto j = std::make_shared<jet>();
    auto j_o = std::make_shared<jet>();

    gausJet(N, *j);
    j_o->nPart = j->nPart-1;
    j_o->pt = std::vector<float>(j->pt.begin(), j->pt.end()-1);
    j_o->eta = std::vector<float>(j->eta.begin(), j->eta.end()-1);
    j_o->phi = std::vector<float>(j->phi.begin(), j->phi.end()-1);
    j_o->sumpt = j->sumpt - j->pt.at(j->pt.size()-1);

    projectedEEC EEC_o = doProjected(j_o, order);
    printf("did j_o\n\n");
    projectedEEC EEC = doProjected(j, order);
    printf("did j\n\n");

    printf("j (%0.3f)\n", j->sumpt);
    std::cout << arma::frowvec(j->pt);
    printf("j_o (%0.3f)\n", j_o->sumpt);
    std::cout << arma::frowvec(j_o->pt);

    printf("DR J\n");
    std::cout << arma::frowvec(*EEC.dRs);
    printf("DR dropped\n");
    std::cout << arma::frowvec(*EEC_o.dRs);
    printf("\n\n");
    printf("WT\n");
    std::cout << arma::frowvec(*EEC.wts);
    printf("WT dropped\n");
    std::cout << arma::frowvec(*EEC_o.wts);
    printf("\n\n");

    printf("COV\n");
    std::cout << arma::trans(EEC.cov->col(j_o->nPart));
}
