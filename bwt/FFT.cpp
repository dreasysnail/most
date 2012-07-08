//
//  FFT.cpp
//  bwt
//
//  Created by zhang yizhe on 12-7-3.
//  Copyright (c) 2012年 SJTU. All rights reserved.
//

#include <iostream>
#include <algorithm>
#include "FFT.h"
using std::iter_swap;

//Bit-reversal Permutation
void FFT::bitReversal (vector<complex<float> > &vecC)
{
    int a, b;
    int p = 0;
    for (int i=1; i<inputSize; i*=2)
    {
        p++;
    }
    for (int i=0; i<inputSize; i++)
    {   
        a = i;
        b = 0;
        for (int j=0; j<p; j++)
        {
            // b = b * 2 + a % 2;
            b = (b << 1) + (a & 1); 
            // a = a / 2;
            a >>= 1;         
        }
        if (b>i)
        {
            iter_swap(vecC.begin()+b,vecC.begin()+i);
        }
    }
}

//Cooley-Tukey algorithm

FFT::FFT(const vector<float> &input){
    
    
    inputSize = input.size();
    if (!checkInput())
        cerr<<"FFT should have input size as a power of 2"<<std::endl;
    for (int i=0; i<input.size(); i++) {
        origin.push_back(complex<float> (input[i],0));
    }
    
}
float FFT::denoise(int rmCount){
    transform(origin,transformed,0);
    float tagNoise=0;
    for (int i=inputSize/2-rmCount; i<inputSize/2; i++) {
        tagNoise += norm(transformed[i]);
        transformed[i]=complex<float> (0.0,0.0);
        transformed[inputSize-i]=complex<float> (0.0,0.0);
        
    }
    transform(transformed,invTrans,1);
    return tagNoise;
}
bool FFT::checkInput(){
    int temp = inputSize;
    while (temp!=1) {
        if (temp%2)
            return false;
        temp = temp/2;
    }
    return true;
}

void FFT::transform(const vector<complex<float> > &from, vector<complex<float> > &to, int mode){
    complex<float> tc;
    //first inputsize/2 n sqrt 's conjugate complex
    vector<complex<float> > w;
    //no check for speed consideration but keep in mind that inputSize should be a power of 2! 
    to = from;
    bitReversal(to);

    float arg;
    if (mode == 0) {
        arg = -2*PI/inputSize;
    }
    else {
        arg = 2*PI/inputSize;
    }
    
    tc = complex<float> (cos(arg),sin(arg));
    
    w.push_back(complex<float> (1.0,0.0));
    
    for (int j=1; j<inputSize/2; j++)
    {
        w.push_back(complex<float> (w.back()*tc));
    }
    for (int m=2; m<=inputSize; m*=2)
    {
        for (int k=0; k<inputSize; k+=m)
        {
            for (int j=0; j<m/2; j++)
            {
                int index1 = k+j;
                int index2 = index1 + m/2;
                //subscript for rotation factor w 
                int t=inputSize*j/m;   
                assert(t<inputSize/2);
                assert(index1<inputSize);  
                assert(index2<inputSize); 
                tc = complex<float> (w[t]*to[index2]);
                to [index2] = to[index1]-tc;
                to [index1] += tc;                
            }
        }
    }
    if (mode != 0) {
        for (int j=0; j<inputSize; j++)
        {
            to[j] = to[j]/float(inputSize);
        }
    }
}

