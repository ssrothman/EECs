#ifndef EECs_ADJ_H
#define EECs_ADJ_H

#include <vector>
#include <Eigen/Dense>

class adjacency{
public:
    adjacency(const Eigen::MatrixXd& ptrans) noexcept {
        data.clear();
        data.resize(ptrans.cols());
        for(unsigned iGen=0; iGen<ptrans.cols(); ++iGen){
            for(unsigned iReco=0; iReco<ptrans.rows(); ++iReco){
                if(ptrans(iReco, iGen)>0){
                    data[iGen].emplace_back(iReco);
                }
            }
        }
    }

    adjacency() noexcept {
        data.clear();
    }

    const std::vector<unsigned>& at(unsigned i) const noexcept{
        return data.at(i);
    }

    std::vector<unsigned> nadj(const std::vector<unsigned>& ord) const noexcept {
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
