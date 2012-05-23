//
//  motif.h
//  bwt
//
//  Created by zhang yizhe on 12-5-20.
//  Copyright (c) 2012年 SJTU. All rights reserved.
//

#ifndef bwt_motif_h
#define bwt_motif_h

#include <iostream>
#include <algorithm>
#include <vector>
#include <math.h>

using namespace std;

#define K 10
#define K_5 9765625
//#define K 6
//#define K_5 15625

#define DELTA 6 



int getIndex();

float fillMotif(const int index);

string translate(const int index);

float pow1(float base,int index);

vector<int> locateMotif(string query,const char T[]);

float testMotifTag(const vector<int> &loci,const vector<int>& tag);

inline bool ascending(int x[],const vector<int> &traverseP,int &travIndex){
    int tempIndex = travIndex;
    for (int i=0; i<traverseP.size(); i++) {
        x[traverseP[traverseP.size()-i-1]]=tempIndex%4+1;
        tempIndex = int(tempIndex/4);
    } 
    travIndex++;
    if (travIndex>pow1(4,traverseP.size())){
        //traversal end 
        travIndex=0;
        return true;
    }
    return false;
}






#endif
