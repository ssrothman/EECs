%module eec_back

%{
  #define SWIG_FILE_WITH_INIT
  #include "eec_back.h"
%}

%include "numpy.i"

%init %{
import_array();
%}

%apply (float* IN_ARRAY2, int DIM1, int DIM2) {(float* jet, int nPart, int nFeat), (float *jets, int nPartTot, int nFeat)}
%apply (int* IN_ARRAY1, int DIM1) {(int* jetIdxs, int nJets), (int* dRIdxs, int nDRIdxs)}

%apply (float* ARGOUT_ARRAY1, int DIM1) {(float* dRs, int nDR), (float* wts, int nWT), (float* dRs, int nDRTot), (float* wts, int nWTTot)}

%include "eec_back.h"
