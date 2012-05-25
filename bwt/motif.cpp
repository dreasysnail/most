//
//  motif.cpp
//  bwt
//
//  Created by zhang yizhe on 12-5-20.
//  Copyright (c) 2012年 SJTU. All rights reserved.
//


#include "motif.h"

using namespace std;

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
    float probThresh=pow1(0.25,K)*DELTA;
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
    if (count>probThresh*GenomeSize) {
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
                temp.push_back('_');
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

float pow1(float base,int index){
    float pow = 1;
    while (index>0) {
        if (index%2==1) {
            pow = pow * base;
        }
        base = base * base;
        index = int(index/2);
    }
    return pow;
}

vector<int> locateMotif(string query,const char T[]){
    vector<int> loci;
    bool flag = true;
    //    cout<<GenomeSize<<endl;
    for (int i=0; i<GenomeSize; i++) {
        for (int index=0; index<query.size(); index++) {
            if(query[index]!='_'&&T[i+index]!=query[index]){
                flag = false;
                break;
            }
        }
        if (flag) {
            loci.push_back(i);
        }
        else {
            flag = true;
        }
    }
    return loci;    
}

float testMotifTag(const vector<int> &loci,const vector<int>& tag){
    float sumMotif=0;
    float sumAround=0;
    for (int i=0; i<loci.size(); i++) {
        if (loci[i]<300||loci[i]>tag.size()-300) {
            continue;
        }
        for (int j=-100; j<100; j++) {
            sumMotif += tag[loci[i]+j];
        }
        for (int j=-300; j<-100; j++) {
            sumAround += tag[loci[i]+j];
        }
        for (int j=100; j<300; j++) {
            sumAround += tag[loci[i]+j];
        }
    }
    //    cout<<sumMotif<<"\t"<<sumAround<<endl;
    float result = logf(sumMotif*2/sumAround);
    return result>0?result:(-result);
}

string antisense(const string& tempString){
    string output("");
    for (int i=tempString.size()-1; i>=0; i--) {
        switch (tempString[i]) {
            case 'A':
                output.push_back('T');
                break;
            case 'T':
                output.push_back('A');
                break;
            case 'C':
                output.push_back('G');
                break;
            case 'G':
                output.push_back('C');
                break;
            case '#':
                output.push_back('#');
                break;
            case 'N':
                output.push_back('N');
                break;
            default:
                break;
        }
    }
    return output;
}
