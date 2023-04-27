#ifndef EECs_ADJ_H
#define EECs_ADJ_H

#include <vector>

class adjacency{
public:
    adjacency(const arma::mat ptrans){
        data.clear();
        data.resize(ptrans.n_cols);
        for(unsigned i=0; i<ptrans.n_cols; ++i){
            for(unsigned j=0; j<ptrans.n_rows; ++j){
                if(ptrans(j, i)>0){
                    data[i].emplace_back(j);
                }
            }
        }
    }

    adjacency(){
        data.clear();
    }

    std::vector<std::vector<unsigned>> data;
};

#endif
