/*
 * The plan here is to write some c++ code
 * Atm just hoping and praying that wrapping it won't be too much of a pain in the ass
 * :)
 */
#include <vector>
#include <tuple>
#include <iostream>
#include <array>
#include <math.h>
#include <stdint.h>
#include <algorithm>

#include "eec_back.h"

//#define VERBOSE

//TODO: think carefully about datatype size
//      Things probably don't need to be full width
//      idx_t needs to be larger than max nParts**2

size_t fact(idx_t n){
  size_t result=1;
  for(idx_t i=2;i<n+1;++i){
    result*=i;
  }
  return result;
}

void printOrd(const std::vector<idx_t> ord){
#ifdef VERBOSE
  std::cout << "(";
  idx_t i =0;
  for(i=0; i<ord.size()-1; ++i){
    std::cout << ord[i] << ", ";
  }
  std::cout << ord[ord.size()-1] << ")";
#endif
}

void fillCompositions(const idx_t n, comp_t& out){
  out.clear();
  for(idx_t i=0; i<n; ++i){
    std::vector<std::vector<idx_t>> next;
    out.push_back(next);
  }

  idx_t a[n+1];
  idx_t k, y, x, l;
  for(idx_t q=0; q<n+1; ++q){
    a[q] = 0;
  }
  k = 1;
  y = n - 1;
  while(k != 0){
    x = a[k - 1] + 1;
    k -= 1;
    while(2 * x <= y){
      a[k] = x;
      y -= x;
      k += 1;
    }
    l = k + 1;
    while(x <= y){
      a[k] = x;
      a[l] = y;
      //yield a[:k + 2];
      do{
        std::vector<idx_t> next;
        next.reserve(k+1);
        for(idx_t q=0; q<k+2; ++q){
          next.push_back(a[q]);
        } 
        out[k+1].push_back(next);
      } while (std::next_permutation(a, a+k+2));
      x += 1;
      y -= 1;
    }
    a[k] = x + y;
    y = x + y - 1;
    //yield a[:k + 1];
    do{
      std::vector<idx_t> next;
      next.reserve(k);
      for(idx_t q=0; q<k+1; ++q){
        next.push_back(a[q]);
      } 
      out[k].push_back(next);
    } while (std::next_permutation(a, a+k+1));
  }
}

void fillSymFactors(const idx_t N, const comp_t& compositions, factor_t& out){
  out.clear();
  out.reserve(N);
  size_t factN = fact(N);
  size_t nextFactor;
  for(idx_t i=0; i<N; ++i){//for each possible length of composition
    std::vector<idx_t> next;
    next.reserve(compositions[i].size());
    for(idx_t j=0; j<compositions[i].size(); ++j){//for each composition of that length
      nextFactor = factN;
      for(idx_t k=0; k<compositions[i][j].size(); ++k){//for each int in that composition
        nextFactor/=fact(compositions[i][j][k]);
      }//end loop over individual composition
      next.push_back(nextFactor);
    }//end loop over compositions of a given length
    out.push_back(next);
  }//end loop over composition length
}

size_t intPow(idx_t a, idx_t b){
  //mine, very stupid
  //should be upgraded to to square multiply
  size_t result=1;
  for(idx_t i=0; i<b; ++i){
    result*=a;
  }
  return result; 
}

coord_t intPow(coord_t a, idx_t b){
  //mine, very stupid
  //should be upgraded to to square multiply
  coord_t result=1;
  for(idx_t i=0; i<b; ++i){
    result*=a;
  }
  return result; 
}

size_t simp(idx_t n, idx_t k){
  size_t result = 1;
  for(idx_t i=0; i<k; ++i){
    result*=n+i;
  }
  result/=fact(k);
  return result;
}

size_t choose(idx_t n, idx_t k){
  if (k > n) return 0;
  if (k * 2 > n) k = n-k;
  if (k == 0) return 1;

  size_t result = n;
  for(size_t i = 2; i <= k; ++i ) {
    result *= (n-i+1);
    result /= i;
  }
  return result;
}

inline void iterate(const idx_t dims, std::vector<idx_t>& ordinates,
                      const idx_t nParts){
  // iterate over dimensions in reverse...
  for (int dim = dims - 1; dim >= 0; --dim){
    if (ordinates[dim] < nParts-dims+dim){
      ++ordinates[dim];
      for (idx_t d2=dim+1; d2<dims; ++d2){
        ordinates[d2] = ordinates[d2-1]+1;
      }
      return;
    }
    if(dim>0){
      ordinates[dim] = ordinates[dim-1]+1;
    } else{
      ordinates[dim] = 0;
    }
  }
}

size_t getIndex(const std::vector<idx_t>& ordinates, const idx_t nPart){
  size_t result=0;
  size_t partPow=1;
  for(idx_t dim=0; dim<ordinates.size(); ++dim){
    result += ordinates[dim]*partPow;
    partPow*=nPart;
  }
  return result;
}

void getPairs(const idx_t nPart, std::vector<pair>* result){
  idx_t i, j;

  result->clear();
  for(i=0; i<nPart-1; ++i){
    for(j=i+1; j<nPart; ++j){
      result->emplace_back(i,j);
    }
  }
}

coord_t dR2(coord_t eta1, coord_t phi1, coord_t eta2, coord_t phi2){
  /*
   * Compute delta R^2 between (eta1, phi1) and (eta2, phi2)
   */
  coord_t deta = eta2-eta1;
  coord_t dphi = phi2-phi1;
  if(dphi > M_PI){
    dphi = 2*M_PI - dphi;
  }else if(dphi < -M_PI){
    dphi = 2*M_PI + dphi; //correct up to sign
  }
  return deta*deta + dphi*dphi;
}

void fillDR2(const coord_t* const jet, const idx_t nPart, 
               coord_t* dRs){
  /*
   *  Compute the pairwise delta r^2 between all the particles in the jet
   *  
   *  jet: (nPart x 3) array of (pT, eta, phi), where pT has been morned to jet pT
   *  nPart: number of jet constituents
   *  dRs: array to fill with dR^2 values
   */
  idx_t i, j, n=0;
  for(i=0;i<nPart-1;++i){
    for(j=i+1;j<nPart;++j){
      dRs[n++] = dR2(jet[3*i+1], jet[3*i+2], jet[3*j+1], jet[3*j+2]);
#ifdef VERBOSE
      std::cout << n-1 << ": dR2 (" << jet[3*i+1] << ", " <<  jet[3*i+2] << ") X (" 
        << jet[3*j+1] << ", " << jet[3*j+2] << "): " << dRs[n-1] << std::endl;
#endif
    }
  }
}

coord_t getWt(const coord_t* const jet, const idx_t nPart, const idx_t N,
                  const std::vector<idx_t>& ord, const idx_t M,
                  const comp_t& compositions, const factor_t& symFactors){
#ifdef VERBOSE
  std::cout << std::endl;
  std::cout << "Getting weight for ";
  printOrd(ord);
  std::cout << std::endl;
#endif
  coord_t result = 0;
  for(idx_t i=0; i<compositions[M-1].size(); ++i){ //for each composition
#ifdef VERBOSE
    std::cout << "\tComposition ";
    printOrd(compositions[M-1][i]);
    std::cout << std::endl;
#endif 
    coord_t nextWt = symFactors[M-1][i];
    for(idx_t j=0; j<compositions[M-1][i].size(); ++j){ //for each element
      nextWt*=intPow(jet[3*ord[j]+0], compositions[M-1][i][j]);
    } //end for each element
    result += nextWt;
  }//end for each composition
  
  return result;
}

coord_t getWt(const coord_t* const jet, const idx_t nPart, const idx_t N,
                  const idx_t i, const idx_t j, const idx_t M,
                  const comp_t& compositions, const factor_t& symFactors){
  std::vector<idx_t> ord;
  ord.push_back(i);
  ord.push_back(j);
  return getWt(jet, nPart, N, ord, 2, compositions, symFactors);
}


void doM(const coord_t* const jet, const idx_t nPart, const idx_t N, 
          const idx_t M,
          coord_t *dRs, 
          const comp_t& compositions, const factor_t& symFactors,
          const std::vector<idx_t>* const cache, const idx_t L,
          std::vector<idx_t>* newCache,
          coord_t *wts){
  /*
   * Compute the weights for correlators with M distinct particles
   * 
   * jet: (nPart x 3) array of (pT, eta, phi), where pT has been normed to jet pT
   * nPart: number of jet constituents
   * N: correlator order
   * M: number of distinct particles in the correlator
   * dRs: array of pairwise dRs
   * compositions: compositions[M] = the M-fold compositions of N. Needed by getWt()
   * cache: indices from previous computation L<M
   *    Pointer instead of ref to allow passing nullptr when N<=2
   * L: number of distinct particles in cache
   * newCache: vector in which to store the new cache indices
   *    Pointer instead of ref to allow passing nullptr if you don't want a newCache
   * wts: vector to add weights to
   */

#ifdef VERBOSE
  std::cout << "doing " << M << "..." << std::endl;
#endif

  //setup new cache
  if(newCache)
    newCache->resize(intPow(nPart, M));
  
  idx_t i, j, idx;
  if(M==2){ 
    /* 
     * Special case for M=2. Ignores cache
     * Ignore cache, L.
     * Pairwise indices are already correct for indexing dRs and wts vectors
     */
    idx = 0;
    for(i=0;i<nPart-1;++i){
      for(j=i+1;j<nPart;++j){
        wts[idx] += getWt(jet, nPart, N, i, j, M, compositions, symFactors);
        if(newCache)
          (*newCache)[i + nPart*j] = idx;
#ifdef VERBOSE
        std::cout << "(" << i << ", " << j << ") " << std::endl
          << "\tdR: " << sqrt(dRs[idx]) << std::endl
          << "\twt: " << wts[idx] << std::endl;
#endif
        ++idx;
      } //end for j
    } // end for i
  } // end if M==2
  else {
    //looping over arbitrary dimensions is complicated...
    if(!cache){
      std::cerr << "Error: cache can only be nullptr when M==2" << std::endl;
      return;
    }
    size_t maxIter = choose(nPart,M);
    std::vector<idx_t> ord(M); //which M-fold combination?
    for(i=0; i<M; ++i){
      ord[i]=i;
    }
    for(size_t iter=0; iter<maxIter; ++iter){//iterate over M-fold combinations
      std::vector<idx_t> ordL(L); //which subcombinations are in the cache? 
      for(i=0; i<L; ++i){
        ordL[i]=i;
      }
      size_t maxIterL = choose(M, L);
      idx_t bestIdx=0;
      coord_t bestDR=-1;
      for(size_t iterL=0; iterL<maxIterL; ++iterL){ //iterate over L-fold combinations of items of the M-fold combination
        std::vector<idx_t> sub(L); //which sub-combination?
        for(i=0; i<L; ++i){
          sub[i] = ord[ordL[i]];
        }
        idx = (*cache)[getIndex(sub, nPart)];//get cached index
        if(dRs[idx]>bestDR){
          bestIdx=idx;
          bestDR=dRs[idx];
        }
        iterate(L, ordL, M);
      }//end iterate over L-fold combnations. We should now have identified the maxDR

      float newWt = getWt(jet, nPart, N, ord, M, compositions, symFactors);
#ifdef VERBOSE
      printOrd(ord);
      std::cout << ":" << std::endl 
        << "\tdR: " << sqrt(dRs[bestIdx]) << std::endl 
        << "\twt: " << newWt << std::endl;
#endif
      wts[bestIdx] += newWt; //placeholder wts call
      if(newCache)
        (*newCache)[getIndex(ord, nPart)] = bestIdx;
      iterate(M, ord, nPart);
    }
  }
}

void eec_onejet(float* jet, int nPart, int nFeat, int N,
                float* dRs, int nDR,
                float* wts, int nWT){
  /*
   * Compute EEC for one jet
   * 
   * double *jet: (nPart x 3) array of (pT, eta, phi) for particles in the jet
   * int nPart: number of particles in the jet
   * int N: correlator order
   *
   * assumes pT has already been normalized by the jet pT
   */
  if(nFeat!=3){
    std::cerr << "Error: nFeat must be 3" << std::endl;
    return;
  }

#ifdef VERBOSE
  std::cout << "doing a jet with " << nPart << " particles, the first of which has pT " << jet[0] << std::endl;
#endif
  
  if(nDR!=nWT || nDR!=choose(nPart, 2)){
    std::cerr << "Error: must have len(dRs) == len(wts) == nPart choose 2" << std::endl;
    return;
  }
  //std::vector<coord_t> dRs(nDR, -1.0); 
  //std::vector<coord_t> wts(nDR, 0.0);

  //fill dR vector
  fillDR2(jet, nPart, dRs);

  for(size_t i=0; i<nDR; ++i){
    wts[i] = 0;
  }

  comp_t compositions;
  fillCompositions(N, compositions);

  factor_t symFactors;
  fillSymFactors(N, compositions, symFactors);

#ifdef VERBOSE
  std::cout<<std::endl;
  for(int i=0; i<N; ++i){
    std::cout << "Compositions of length " << i+1 << std::endl;
    for(int j=0;j<compositions[i].size();++j){
      std::cout << "\t";
      printOrd(compositions[i][j]);
      std::cout << " " << symFactors[i][j] << std::endl;
    }
  }
  std::cout<<std::endl;
#endif

  std::vector<idx_t> cache(intPow((idx_t)nPart,2), 0);

  doM(jet, nPart, N, 2, dRs, compositions, symFactors, nullptr, 0, &cache, wts);
  for(idx_t M=3; M<=N; ++M){
    doM(jet, nPart, N, M, dRs, compositions, symFactors, &cache, 2, nullptr, wts);
  }

#ifdef VERBOSE
  std::cout<<std::endl<<"dR\twt"<<std::endl;
  for(idx_t i=0; i<nDR; ++i){
    std::cout << sqrt(dRs[i]) << ",\t" << wts[i] << std::endl;
  }
#endif
  for(idx_t i=0; i<nDR; ++i){
    if(abs(wts[i])>1e10){
      std::cout <<"found a bad one"<<std::endl<<std::endl;
    }
  }

  /*
   * Plan:
   *  For M-1:
   *    store where in the dR array each combination points
   *  For M:
   *    for each combination:
   *      identify the 3 component (M-1) combinations
   *      identify which has max dR
   *      record index in dR vector
   *      compute appropriate weights by iterating over partitions
   *      add to wTs vector in appropriate index
   *
   *  If M>some limit (depending on maximum allowable memory consumption, 5 could be a reasonable choice)
   *    Stop caching dR location of each combination
   *    Instead, fallback on largest cached dR combination list, and do something clever
   *    
   */

  /*
   * Food for thought: 
   *  might as well compute all lower order correlators while we're at it?
   *  Histograms are very memory-efficient
   */

  /*
   * Probably most efficient to pass in pT, eta, phi individually 
   * rather than stacking them
   */
}

void eec(float* jets, int nPartTot, int nFeat,
          int* jetIdxs, int nJets,
          int N,
          float* dRs, int nDRTot,
          float* wts, int nWTTot,
          int* dRIdxs, int nDRIdxs){
  //TODO: need some size checks for everything else too
  if(nFeat!=3){
    std::cerr << "Error: nFeat must be 3" << std::endl;
    return;
  }
  for(size_t i=0; i<nDRTot; ++i){
    wts[i]=0;
    dRs[i]=0;
  }


  size_t prevJetIdx = 0, prevDRIdx = 0;
  idx_t nPart, nDR;
  for(size_t i=0; i<nJets; ++i){
    nPart = jetIdxs[i] - prevJetIdx;
    nDR = dRIdxs[i] - prevDRIdx;
    eec_onejet(&jets[3*prevJetIdx], nPart, nFeat, N, 
                &dRs[prevDRIdx], nDR, 
                &wts[prevDRIdx], nDR);
    prevJetIdx = jetIdxs[i];
    prevDRIdx = dRIdxs[i];
  }
}

/*
int main(){
  coord_t jet[15];
  jet[3*0+0] = 1.0;
  jet[3*1+0] = 2.0;
  jet[3*2+0] = 0.5;
  jet[3*3+0] = 2.0;
  jet[3*4+0] = 3.0;

  jet[3*0+1] = 0.0;
  jet[3*1+1] = 0.1;
  jet[3*2+1] = 0.4;
  jet[3*3+1] = 1.0;
  jet[3*4+1] = 0.4;

  jet[3*0+2] = 0.0;
  jet[3*1+2] = 0.2;
  jet[3*2+2] = 0.4;
  jet[3*3+2] = 0.0;
  jet[3*4+2] = -0.5;

  eec_onejet(jet, 5, 3, 5);

  return 0;
}
*/
