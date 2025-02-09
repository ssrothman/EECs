#ifndef SROTHMAN_EECS_STANDALONE_ADJACENCY_H
#define SROTHMAN_EECS_STANDALONE_ADJACENCY_H

#include <vector>
#include <Eigen/Dense>

namespace standaloneEEC{
    struct neighbor{
        unsigned idx;
        double wt;
    };
    using neighborhood = std::vector<neighbor>; // vector of neighbors

    class adjacency {
    public:

        std::vector<neighborhood> adj; //the adjacency list
        std::vector<bool> hasMatch; // = adj[i].size() > 0
        
        adjacency(const Eigen::MatrixXd& tmat) {
            adj.resize(tmat.rows());
            hasMatch.resize(tmat.rows());
            for (unsigned i = 0; i < tmat.rows(); i++) {
                hasMatch[i] = false;
                for (unsigned j = 0; j < tmat.cols(); j++) {
                    if (tmat(i, j) != 0) {
                        adj[i].push_back({j, tmat(i, j)});
                        hasMatch[i] = true;
                    };
                };
            };
        }; 
    };
};

#endif
