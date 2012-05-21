//
//  motif.cpp
//  bwt
//
//  Created by zhang yizhe on 12-5-20.
//  Copyright (c) 2012年 SJTU. All rights reserved.
//

#include "motif.h"
int motif[K_5]={0};
int GenomeSize;


float fillMotif(const int index){
    int tempIndex=index;
    int currentNum;
    int kmer[K]={0};
    int count = 0;
    int query=0;
    int traverseIndex = 0;
    //for 0: 1-4 traverse
    vector<int> traversalPoint;
    int x[K]={0};
    //0:* 1:A 2:C 3:G 4:T 
    float probThresh=pow(0.25,K)*DELTA;
    //    cout<<pow(0.25, K)<<endl;
    for (int i=0; i<K; i++){
        currentNum = tempIndex%5;
        kmer[K-i-1] = currentNum;
        tempIndex = int(tempIndex/5);
        if (currentNum==0){
            probThresh *= 4;
            traversalPoint.push_back(K-i-1);
            x[K-i-1] = 1;
        }
    }
    sort(traversalPoint.begin(), traversalPoint.end());
    while (true) {
        for (int j=0; j<K; j++) {
            if (kmer[j]==0) {
                query = query*5+x[j];
            }
            else {
                query = query*5+kmer[j];
            }
            
        }
        count += motif[query];
        query = 0;
        //if over
        
        if(ascending(x,traversalPoint,traverseIndex))
            break;
    }
    
    /* 
    for (int j=0; j<(K_5); j++) {
        int bitPos=K-1;
        bool skip=false;
        int query=j;
        for (int i=0; i<K; i++){
            currentNum = query%5;
            query = int(query/5);
            if (currentNum==0||(kmer[bitPos]!=0&&kmer[bitPos]!=currentNum)) {
                skip = true;
                break;
            }
            bitPos--;
        }
        if (skip) {
            continue;
        }
        //if match
        count += motif[j];   
    }
    */
    motif[index] = count;
    if (count==0||count>probThresh*GenomeSize) {
        return int(count*1000/probThresh/GenomeSize)/1000.0;
    }
    else {
        return 0;
    }
    
}


string translate(const int index){
    string temp;
    int query = index;
    int currentNum;
    for (int i=0; i<K; i++){
        currentNum = query%5;
        query = int(query/5);
        switch (currentNum) {
            case 0:
                temp.push_back('X');
                break;
            case 1:
                temp.push_back('A');
                break;
            case 2:
                temp.push_back('C');
                break;
            case 3:
                temp.push_back('G');
                break;
            case 4:
                temp.push_back('T');
                break;
            default:
                break;
        }
    }
    return temp;
}