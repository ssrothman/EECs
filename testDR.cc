#include <stdio.h>
#include "eec.h"

int main(){
    std::vector<float> Es = {1.0, 2.0, 3.0, 2.0, 5.0};
    std::vector<float> etas = {0.0, 0.2, -0.4, 1.0, 0.1};
    std::vector<float> phis = {3.0, -3.0, 3.1, 2.7, 2.8};

    vecND::nodiagvec dR2s(5, 2);
    getDR2(dR2s, etas, phis);

    auto ord = dR2s.ord0();
    bool loop;
    for(size_t i=0, loop=true; loop; loop=dR2s.iterate(ord), ++i){
        printOrd(ord);
        printf(": %0.3f\n", dR2s.at(i));
    }
}
