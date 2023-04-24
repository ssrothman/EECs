#ifndef EECs_ADJ_H
#define EECs_ADJ_H

#include <vector>

class adjacency{
public:
    adjacency(std::shared_ptr<const arma::mat> ptrans):
        data()
    {
        printf("constructing adjacency\n");
        if(ptrans){
            printf("about to clear data\n");
            data.clear();
            printf("about to resize data\n");
            data.resize(ptrans->n_cols);
            printf("about to enter loop\n");
            for(unsigned i=0; i<ptrans->n_cols; ++i){
                for(unsigned j=0; j<ptrans->n_rows; ++j){
                    if((*ptrans)(j, i)>0){
                        printf("about to emplace back\n");
                        data[i].emplace_back(j);
                    }
                }
            }
        }
        printf("constructed adjacency\n");
    }

    std::vector<std::vector<unsigned>> data;
};

#endif
