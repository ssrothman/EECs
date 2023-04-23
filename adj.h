#ifndef EECs_ADJ_H
#define EECs_ADJ_H

class adjacency{
public:
    adjacency(std::shared_ptr<const arma::mat> ptrans){
        if(ptrans){
            data.clear();
            data.resize(ptrans->n_cols);
            for(unsigned i=0; i<ptrans->n_cols; ++i){
                for(unsigned j=0; j<ptrans->n_rows; ++j){
                    if((*ptrans)(j, i)>0){
                        data[i].emplace_back(j);
                    }
                }
            }
        }
    }

    std::vector<std::vector<unsigned>> data;
};

#endif
