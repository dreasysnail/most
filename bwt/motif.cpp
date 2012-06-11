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

void Motif::fillMotif(const vector<Motif>& allmotifs){
    explainMotif();
    vector<int>::iterator it;
    int count=0;
    //float motifProb;
    for (it=expMotifs.begin(); it!=expMotifs.end();it++ ) {
        count += motif[(*it)];
        //add generalization's prob
        probThresh += allmotifs[(*it)].probThresh;
    } 
    float Thresh = probThresh*GenomeSize*DELTA;
    if (count>Thresh){
        score = int(count*1000/Thresh)/1000.0;
    }
    else {
        score = 0;
    }
    
}

vector<int> Motif::explainMotif(){
    //find all possible explanation(generalization)
    vector<int> tempV;
    expMotifs.swap(tempV);
    if (noWildcard()) {
        expMotifs.push_back(index);
        return expMotifs;
    }
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


void Motif::locateMotif(const char T[]){
    //implement1 low efficiency
    bool flag = true;
    loci.clear();
    //cout<<GenomeSize<<endl;
    for (int i=0; i<GenomeSize; i++) {
        for (int index=0; index<query.size(); index++) {
            if((query[index]!='*'&&T[i+index]!=query[index])||T[i+index]=='#'||T[i+index]=='N'){
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
    float sumMotif=0;
    float sumAround[2*SAMPLESIZE]={0};
    //float sumShoulder=0;
    int counter = loci.size();
    for (int i=0; i<loci.size(); i++) {
        if (loci[i]<(2*SAMPLESIZE+1)*dist+SAMPLESIZE*K+offset||loci[i]>tag.size()-(2*SAMPLESIZE+1)*dist-(SAMPLESIZE+1)*K-1-offset) {
            counter--;
            continue;
        }
        for (int j=-dist; j<dist+K; j++) {
            sumMotif += tag[loci[i]+j];
        }
        //     for (int j=-3*dist-K-offset1; j<-dist-offset1; j++) sumShoulder += tag[loci[i]+j];
        //     for (int j=dist+K+offset1; j<3*dist+2*K+offset1; j++) sumShoulder += tag[loci[i]+j];
        for (int l=0; l<SAMPLESIZE; l++) {
            for (int j=-(2*l+3)*dist-(l+1)*K-offset; j<-(2*l+1)*dist-l*K-offset; j++) {
                sumAround[l] += tag[loci[i]+j];
            }
            for (int j=(2*l+1)*dist+(l+1)*K+offset; j<(2*l+3)*dist+(l+2)*K+offset; j++) {
                sumAround[l+SAMPLESIZE] += tag[loci[i]+j];
            }
        }
    }
    float totalAround=0;
    float variance = 0;
    for (int i=0; i<SAMPLESIZE*2; i++) {
        totalAround += sumAround[i];
    }
    float average = (totalAround+sumMotif)/(SAMPLESIZE*2+1);
    for (int i=0; i<SAMPLESIZE*2; i++) {
        variance += (sumAround[i]-average)*(sumAround[i]-average);
    }
    variance += (sumMotif-average)*(sumMotif-average);
    variance = sqrt(variance/SAMPLESIZE*2);
    //   vector<float> temp(sumAround,sumAround+2*SAMPLESIZE);
    //cout<<sumMotif<<temp<<average<<"var "<<variance<<" Counter"<<counter<<" height"<<logf(sumMotif/float(counter))<<endl;
    //signif.push_back((sumMotif-average)/variance*logf(sumMotif/float(counter))) ;
    signif.push_back((sumMotif-average)/variance);
    if (signif.back()<0) {
        signif.back()=-signif.back();
    }
}


void Motif::initProb(const genomeRegions& gR,int order){
    if (order==0) {
        float motifProb = 1;
        for (int i=0; i<K; i++){
            motifProb *= gR.prob[0][queryNum[i]];        
        }
        probThresh = motifProb;  
    }
    else if (order==1) {
        float motifProb = gR.prob[0][queryNum[0]];
        for (int i=1; i<K; i++){
            motifProb *= gR.prob[queryNum[i-1]+1][queryNum[i]];       
        }
        probThresh = motifProb;  
    }
    else if (order==-1) {
        probThresh = pow1(0.25, K);
    }
      
}

void Motif::initPWM(){
    for (int i=0; i<K; i++) {            
        switch (query[i]) {
            case 'A':
                pwm[0][i] += motif[index];
                break;
            case 'C':
                pwm[1][i] += motif[index];
                break;
            case 'G':
                pwm[2][i] += motif[index];
                break;
            case 'T':
                pwm[3][i] += motif[index];
                break;
            default:
                break;
        }
    }
}
void Motif::calPWM(const vector<Motif>& allmotifs){
    explainMotif();
    vector<int>::iterator it;
    //float motifProb;
    for (it=expMotifs.begin(); it!=expMotifs.end();it++ ) {
        for (int i=0; i<K; i++) {
            for (int j=0; j<4; j++) {
                pwm[j][i] += allmotifs[(*it)].pwm[j][i];
            }
        }
    } 
}

bool Motif::ascending(int x[],const vector<int> &traverseP,int &travIndex){
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

void Motif::sumScore(){
    overallScore = score*score;
    float tempScore = 0;
    for(int i=0;i<signif.size();i++){
        tempScore += signif[i];
    }
    overallScore *= tempScore/10.0;
}

ostream &operator<<( ostream &s, const Motif &motif ){
    s<<">"<<motif.query<<"_Conscore_"<<motif.score<<"_Signif_"<<motif.signif;
    for (int i=0; i<4; i++) {
        for (int j=0; j<K; j++) {
            s<<motif.pwm[i][j]<<'\t';
        }
        s<<endl;
    }
    return s;
}