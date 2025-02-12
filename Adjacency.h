#ifndef SROTHMAN_EECS_STANDALONE_ADJACENCY_H
#define SROTHMAN_EECS_STANDALONE_ADJACENCY_H

#include <vector>
#include <Eigen/Dense>
#include "usings.h"

namespace EEC{
    struct neighbor{
        unsigned idx;
        double wt;
    };
    using neighborhood = std::vector<neighbor>; // vector of neighbors

    class Adjacency {
    public:
        
        Adjacency(const Eigen::MatrixXd& tmat)  noexcept {
            unsigned ngen = tmat.cols();
            unsigned nreco = tmat.rows();

            adj.resize(ngen);
            hasMatch.resize(ngen);
            for(unsigned i=0; i<ngen; ++i){
                adj[i].clear();
                hasMatch[i] = false;
            }

            for (unsigned igen = 0; igen < ngen; igen++) {
                hasMatch[igen] = false;
                for (unsigned ireco = 0; ireco < nreco; ireco++) {
                    if (tmat(ireco, igen) != 0) {
                        adj[igen].push_back({ireco, tmat(ireco, igen)});
                        hasMatch[igen] = true;
                    };
                };
            };
        }; 

        void print() const noexcept {
            printf("adj.size() = %lu\n", adj.size());
            printf("hasMatch.size() = %lu\n", hasMatch.size());
            for(unsigned i=0; i<adj.size(); ++i){
                printf("adj[%u].size() = %lu\n", i, adj[i].size());
                printf("hasMatch[%u] = %d\n", i, hasMatch[i]);
                for(unsigned j=0; j<adj[i].size(); ++j){
                    printf("adj[%u][%u] = {%u, %g}\n", i, j, adj[i][j].idx, adj[i][j].wt);
                }
            }
        }

        const neighborhood& get_neighborhood(const unsigned i) const noexcept {
            /*printf("getting neighborhood of %u\n", i);
            fflush(stdout);
            printf("\tsize: %lu\n", adj[i].size());
            fflush(stdout);
            if(adj[i].size() > 0){
                printf("\tadj[%u][0] = {%u, %g}\n", i, adj[i][0].idx, adj[i][0].wt);
                fflush(stdout);
            }*/
            return adj[i];
        }

        const bool has_match(const unsigned i) const noexcept {
            return hasMatch[i];
        }
    private:
        std::vector<neighborhood> adj; //the adjacency list
        std::vector<bool> hasMatch; // = adj[i].size() > 0
    };
};

#endif
