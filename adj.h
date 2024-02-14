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

    const std::vector<unsigned>& at(unsigned i) const{
        return data.at(i);
    }

    std::vector<unsigned> nadj(const std::vector<unsigned>& ord) const {
        std::vector<unsigned> ans(ord.size(), 0);
        for(unsigned i=0; i<ord.size(); ++i){
            ans[i] = at(ord[i]).size();
            if(ans[i] == 0){
                return {};
            }
        }
        return ans;
    }

private:
    std::vector<std::vector<unsigned>> data;
};

#endif
