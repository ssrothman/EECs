#ifndef EECS_FAST_UTIL_H
#define EECS_FAST_UTIL_H

#include <boost/histogram.hpp>
#include <boost/multi_array.hpp>

#include "usings.h"

#include "adj.h"

namespace fastEEC{
    template <typename AX>
    unsigned extent(const AX& ax){
        return histogram::axis::traits::extent(ax);
    }

    template <typename T>
    void dumpToFile(const std::vector<T>& vec, const std::string& filename){
        std::ofstream file;
        file.open(filename);
        file << std::setprecision(10);
        file << std::scientific;
        file << std::showpos;
        file << std::uppercase;
        file << "[";
        file << vec.size();
        file << "]\n";
        for(unsigned i=0; i<vec.size(); ++i){
            file << vec[i];
            if (i < vec.size()-1){
                file << ", ";
            }
        }
        file.close();
    }

    template <typename T, long unsigned D>
    void dumpToFile(const boost::multi_array<T, D>& arr, const std::string& filename){
        std::ofstream file;
        file.open(filename);
        file << std::setprecision(10);
        file << std::scientific;
        file << std::showpos;
        file << std::uppercase;
        file << "[";
        for(unsigned d=0; d<D; ++d){
            file << arr.shape()[d];
            if (d < D-1){
                file << ", ";
            }
        }
        file << "]\n";
        for(unsigned i=0; i<arr.num_elements(); ++i){
            file << arr.data()[i];
            if (i < arr.num_elements()-1){
                file << ", ";
            }
        }
        file.close();
    }

    template <typename T>
    void addInPlace(T& lhs, const T& rhs){
        lhs += rhs;
    }

    template <typename T>
    void addInPlace(std::vector<T>& lhs,
                    const std::vector<T>& rhs){
        if(lhs.size() != rhs.size()){
            throw std::invalid_argument("Vectors must be the same size");
        }
        for(unsigned i=0; i<lhs.size(); ++i){
            lhs[i] += rhs[i];
        }
    }

    template <typename T, long unsigned D>
    void addInPlace(boost::multi_array<T, D>& lhs,
                    const boost::multi_array<T, D>& rhs){
        if(lhs.num_elements() != rhs.num_elements()){
            throw std::invalid_argument("Arrays must be the same shape");
        }
        for(unsigned i=0; i<lhs.num_elements(); ++i){
            lhs.data()[i] += rhs.data()[i];
        }
    }

    template <typename T>
    unsigned getIndex(const T& val, const axisptr& ax){
        return static_cast<unsigned>(ax->index(val) + 1);
    }

    template <typename VEC>
    void printVec(const VEC& vec){
        for (const auto& v : vec){
            std::cout << v << ", ";
        }
    }

    enum normType {
        RAWPT, 
        CORRPT,
        SUMPT, 
    };

    template <typename T>
    void getFloatDRs(multi_array<T, 2>& ans, const jet& J){
        ans.resize(extents[J.nPart][J.nPart]);
        for (unsigned i0=0; i0<J.nPart; ++i0){
            for (unsigned i1=0; i1<J.nPart; ++i1){
                T deltaR = dR(J.particles[i0], J.particles[i1]);
                if(deltaR < 1e-10){
                    deltaR = 0;
                }
                ans[i0][i1] = deltaR;
            }
        }
    }
    
    inline void getDRbins(umat& ans, const jet& J, const axisptr& ax){
        ans.resize(extents[J.nPart][J.nPart]);
        for (unsigned i0=0; i0<J.nPart; ++i0){
            for (unsigned i1=0; i1<J.nPart; ++i1){
                double deltaR = dR(J.particles[i0], J.particles[i1]);
                if (deltaR < 1e-10){
                    deltaR = 0;
                }
                unsigned idx = static_cast<unsigned>(ax->index(deltaR) + 1);
                ans[i0][i1] = idx;
            }
        }
    }

    inline double getNormFact(const jet& J, const normType& nt){
        switch (nt){
            case RAWPT:
                return J.rawpt;
            case SUMPT:
                return J.sumpt;
            case CORRPT:
                return J.pt;
            default:
                throw std::invalid_argument("Invalid normType");
        }
    }

    template <typename T>
    void getEs(vector<T>& ans, const jet& J, const normType nt){
        ans.clear();
        ans.reserve(J.nPart);

        double normFact = getNormFact(J, nt);
        for(unsigned i=0; i<J.nPart; ++i){
            ans.push_back(J.particles.at(i).pt/normFact);
        }
    }

    template <typename T>
    void getEtasPhis(vector<T>& etas, vector<T>& phis, const jet& J){
        etas.clear();
        phis.clear();

        etas.reserve(J.nPart);
        phis.reserve(J.nPart);

        for(unsigned i=0; i<J.nPart; ++i){
            etas.push_back(J.particles.at(i).eta);
            phis.push_back(J.particles.at(i).phi);
        }
    }

    template <typename T>
    void getPtrans(multi_array<T, 2>& ans, const arma::mat& ptrans){//NB we transpose for faster iteration
        ans.resize(extents[ptrans.n_cols][ptrans.n_rows]);
        for(unsigned i=0; i<ptrans.n_rows; ++i){
            for(unsigned j=0; j<ptrans.n_cols; ++j){
                ans[j][i] = ptrans(i,j);
            }
        }
    }
}

#endif
