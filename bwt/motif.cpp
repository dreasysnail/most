//
//  motif.cpp
//  bwt
//
//  Created by zhang yizhe on 12-5-20.
//  Copyright (c) 2012年 SJTU. All rights reserved.
//


#include "motif.h"
#include "bwt.h"

using namespace std;

int GenomeSize;
int Motif::motif[K_5] = {0};

void Motif::fillMotif(){
    explainMotif();
    vector<int>::iterator it;
    int count=0;
    float probThresh=pow1(0.25,K)*DELTA;
    for (int i=0; i<K; i++){
        if (query[i]=='_'){
            probThresh *= 4;
        }
    }
    for (it=expMotifs.begin(); it!=expMotifs.end();it++ ) {
        count += motif[(*it)];
    }     
    if (count>probThresh*GenomeSize){
        //        return count;
    score = int(count*1000/probThresh/GenomeSize)/1000.0;
    }
    else {
        score = 0;
    }
    
}

vector<int> Motif::explainMotif(){
    //find all possible explanation
    int tempIndex=index;
    int query=0;
    int currentNum;
    int kmer[K]={0};
    int traverseIndex = 1;
    //for 0: 1-4 traverse
    vector<int> traversalPoint;
    int x[K]={0};
    //0:* 1:A 2:C 3:G 4:T 
    for (int i=0; i<K; i++){
        currentNum = tempIndex%5;
        kmer[K-i-1] = currentNum;
        tempIndex = int(tempIndex/5);
        if (currentNum==0){
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
        expMotifs.push_back(query);
        query = 0;
        //if over
        
        if(ascending(x,traversalPoint,traverseIndex))
            break;
    }
    return expMotifs;
}


string Motif::translate(int i){
    string temp;
    int query = i;
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

void Motif::locateMotif(const char T[]){
    //implement1 low efficiency
    bool flag = true;
    loci.clear();
    //cout<<GenomeSize<<endl;
    for (int i=0; i<GenomeSize; i++) {
        for (int index=0; index<query.size(); index++) {
            if((query[index]!='_'&&T[i+index]!=query[index])||T[i+index]=='#'||T[i+index]=='N'){
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
}






void Motif::testMotifTag(const vector<int>& tag){
    const int dist=20;
	const int offset=2;
    float sumMotif=0;
    float sumAround=0;
    for (int i=0; i<loci.size(); i++) {
        if (loci[i]<3*dist+K+offset||loci[i]>tag.size()-3*dist-2*K-1-offset) {
            continue;
        }
        for (int j=-dist; j<dist+K; j++) {
            sumMotif += tag[loci[i]+j];
        }
        for (int j=-3*dist-K-offset; j<-dist-offset; j++) {
            sumAround += tag[loci[i]+j];
        }
        for (int j=dist+K+offset; j<3*dist+2*K+offset; j++) {
            sumAround += tag[loci[i]+j];
        }
    }
    //   cout<<"summotif:"<<sumMotif<<"\tsumaround:"<<sumAround<<endl;
    float result = logf(sumMotif*2/sumAround);
    signif = result>0?result:(-result);
}

void Motif::printMotif(){
    cout<<query<<"\t"<<score<<"\t"<<signif<<endl;
}
