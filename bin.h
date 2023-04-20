#ifndef EECs_BIN_H
#define EECs_BIN_H

#include <boost/histogram.hpp>

template <typename H, typename X=float, typename W=float>
void fill1d(H& hist, 
            const std::vector<X>& x,
            const std::vector<W>& w){
    for(unsigned i=0; i<x.size(); ++i){
        hist(x[i], boost::histogram::weight(w[i]));
    }
}

template <typename H, typename X1=float, typename X2=float, typename W=float>
void fill2d(H& hist, 
            const std::vector<X1>& x1, 
            const std::vector<X2>& x2, 
            const arma::Mat<W>& w){
    for(unsigned i1=0; i1<x1.size(); ++i1){
        for(unsigned i2=0; i2<x2.size(); ++i2){
            hist(x1[i1], x2[i2], boost::histogram::weight(w(i1, i2)));
        }
    }
}

#endif
