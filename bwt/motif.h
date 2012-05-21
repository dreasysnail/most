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

#define K 6
#define K_5 15625
#define DELTA 2.5



int getIndex();
float fillMotif(const int index);

string translate(const int index);

inline bool ascending(int x[],const vector<int> &traverseP,int &travIndex){
    int tempIndex = travIndex;
    for (int i=0; i<traverseP.size(); i++) {
        x[traverseP[traverseP.size()-i-1]]=tempIndex%4+1;
        tempIndex = int(tempIndex/4);
    } 
    travIndex++;
    if (travIndex>pow(4,traverseP.size())){
        //traversal end 
        travIndex=0;
        return true;
    }
    return false;
}


#endif
