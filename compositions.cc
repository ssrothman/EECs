#include "compositions.h"


std::shared_ptr<std::vector<comp_t>> getCustomComps(unsigned p1, unsigned p2){

    auto result = std::make_shared<std::vector<comp_t>>();
    result->resize(1);
    result->at(0).resize(2);
    std::vector<unsigned> c0 = {p1+p2};
    if(p2!=p1){
        result->at(0)[0].emplace_back(c0, 2u);
        std::vector<unsigned> c12 = {p1, p2};
        std::vector<unsigned> c21 = {p2, p1};
        result->at(0)[1].emplace_back(c12, 1u);
        result->at(0)[1].emplace_back(c21, 1u);
    } else {
        result->at(0)[0].emplace_back(c0, 1u);
        std::vector<unsigned> c11 = {p1, p2};
        result->at(0)[1].emplace_back(c11, 2u);
    }
    return result;
}

