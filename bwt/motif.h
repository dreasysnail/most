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
#include "common.h"

class Motif{
    public:
    int pwm[4][K];
    string query;
    int index;
    vector<int> loci;
    float signif;
    vector<int> expMotifs;
    //conserve score;
    float score;
    //index
    Motif(int i):index(i),signif(0){query = translate(index);}
    void fillMotif();
    vector<int> explainMotif();
    static string translate(int i);
    void locateMotif(const char T[]);
    void testMotifTag(const vector<int>& tag);
    string antisense(const string& tempString);
    bool ascending(int x[],const vector<int> &traverseP,int &travIndex);
    static int motif[K_5];
    void printMotif();
};




inline bool Motif::ascending(int x[],const vector<int> &traverseP,int &travIndex){
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
