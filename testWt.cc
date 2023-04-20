#include <stdio.h>
#include "eec.h"
#include "simon_util_cpp/combinatorics.h"
#include "toyjets/gaus.h"

int main(){
    unsigned N = 5u;
    unsigned order = 4u;
    auto j = std::make_shared<jet>();

    gausJet(N, *j);

    projectedEEC result = doProjected(j, order);

    float sumWt = result.ptAtZero;
    for(unsigned i=0; i<result.wts->size(); ++i){
        sumWt += result.wts->at(i);
    }

    printf("Total wt = %0.3g\n", sumWt);
}
