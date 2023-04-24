#include <stdio.h>
#include "toyjets/gen.h"
#include "toyjets/gaus.h"
#include "eec_oo.h"

int main(){
    unsigned N=5u;
    unsigned order=2u;
    auto j = std::make_shared<jet>();
    auto j_o = std::make_shared<jet>();

    gausJet(N, *j_o);
    auto ptrans = std::make_shared<arma::mat>(genJet(*j_o, *j, 
                    0.15, 0.05, 0.05,
                    0.15, 0.80, 0.00, 0.10, 0.3));

    printf("j\n");
    std::cout << arma::trans(j->ptvec())/j->sumpt;
    printf("j_o\n");
    std::cout << arma::trans(j_o->ptvec())/j_o->sumpt;
    printf("mat * j\n");
    std::cout << arma::trans(*ptrans *j->ptvec()/j->sumpt);
    printf("\n\n");
    std::cout << *ptrans;
    printf("\n\n");

    auto customComps = getCustomComps(2, 1);
    NonIRCEECCalculator trans(j, order, ptrans, j_o, customComps);
    NonIRCEECCalculator reco(j_o, order, customComps);

    printf("made\n");
    trans.run();
    printf("did trans\n\n");
    reco.run();
    printf("did reco\n\n");

    printf("GEN DR\n");
    std::cout << arma::rowvec(trans.getdRs());
    printf("\n");
    printf("RECO DR\n");
    std::cout << arma::rowvec(trans.getdRs_J2());
    printf("RECO DR 2\n");
    std::cout << arma::rowvec(reco.getdRs());
    printf("\n\n");
    printf("GEN WT\n");
    std::cout << arma::rowvec(trans.getwts(2));
    printf("RECO WT\n");
    std::cout << arma::rowvec(reco.getwts(2));
    printf("\n\n");
    std::cout << trans.getTransfer(2);
    printf("\n\n");
    std::cout << arma::trans(arma::sum(trans.getTransfer(2), 1));
}
