#include <stdio.h>
#include "eec.h"
#include "simon_util_cpp/combinatorics.h"

int main(){
    unsigned N = 5;
    comp_t comps;
    fillCompositions(N, comps);
    for(unsigned i=0; i<comps.size(); ++i){
        printf("ORDER %d\n", i);
        for(unsigned j=0; j<comps[i].size(); ++j){
            printOrd(comps[i][j].composition);
            printf(": %u\n", comps[i][j].factor);
        }
        printf("\n");
    }
}
