#include <stdio.h>
#include "eec.h"
#include "toyjets/gen.h"
#include "toyjets/gaus.h"

int main(){
    unsigned N=5u;
    unsigned order=2u;
    auto j = std::make_shared<jet>();
    auto j_o = std::make_shared<jet>();

    gausJet(N, *j_o);
    auto ptrans = std::make_shared<arma::fmat>(genJet(*j_o, *j, 
                    0.15, 0.05, 0.05,
                    0.20, 0.80, 0.10, 0.20, 0.2));

    projectedEEC EEC = doProjected(j, order);
    printf("did gen\n\n");
    projectedEEC EEC_o = doProjected(j_o, order);
    printf("did reco\n\n");
    projectedEEC trans = doProjected(j, order, ptrans, j_o);
    printf("did both\n\n");

    printf("j\n");
    std::cout << arma::frowvec(j->pt)/j->sumpt;
    printf("j_o\n");
    std::cout << arma::frowvec(j_o->pt)/j_o->sumpt;
    printf("mat * j\n");
    std::cout << arma::trans(*ptrans *arma::fvec(j->pt)/j->sumpt);
    printf("\n\n");

    std::cout << *ptrans;
    printf("\n\n");

    printf("EEC_Gen DR\n");
    std::cout << arma::frowvec(EEC.dR2s->data());
    printf("EEC_Gen DR 2\n");
    std::cout << arma::frowvec(trans.dR2s->data());
    printf("\n");
    printf("EEC_Reco DR\n");
    std::cout << arma::frowvec(EEC_o.dR2s->data());
    printf("EEC_Reco DR 2\n");
    std::cout << arma::frowvec(trans.dR2s_o->data());
    printf("\n\n");
    printf("EEC_Gen WT\n");
    std::cout << arma::frowvec(EEC.wts->data());
    printf("EEC_Gen WT 2\n");
    std::cout << arma::frowvec(trans.wts->data());
    printf("EEC_Reco WT\n");
    std::cout << arma::frowvec(EEC_o.wts->data());
    printf("\n\n");

    std::cout << *trans.transfer;
}
