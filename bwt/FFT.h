//
//  FFT.h
//  MOST
//
//  Created by zhang yizhe on 12-7-3.
//  Copyright (c) 2012年 SJTU. All rights reserved.
//

#ifndef bwt_FFT_h
#define bwt_FFT_h
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex>
#include "common.h"
using std::vector;
using std::complex;

class FFT {
public:
    int inputSize;
    vector<complex<float> > transformed;
    vector<complex<float> > origin;
    vector<complex<float> > invTrans;
    FFT(const vector<float> &input);
    bool checkInput();
    //mode =0 for fft mode=1 for ifft
    void transform(const vector<complex<float> > &from, vector<complex<float> > &to, int mode);
    float denoise(int rmCount);
    void bitReversal (vector<complex<float> > &vecC);
};

#endif
