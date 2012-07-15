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
#include <iomanip>
//#include <stdio.h>
//#include <stdlib.h>
#include "common.h"
using std::cout;
using std::endl;
using std::string;
using std::vector;

class Motif{
    public:
    string query;
    //char alp2num(query[i-1])[K];
    int index;
    vector<int> pwm[4];
    //bool overThresh;
    vector<int> loci;
    vector<float> tagBiPeak;
    vector<float> tagSymmetry;
    float motifProb;
    //each vector is a bin   sumBin[bin][tag]
    vector<float> sumBin[2*SAMPLESIZE+1];
    vector<float> signalIntensity;
    //vector<int> expMotifs;
    //std::vector<string> Clusterexp;
    //conserve score;
    float score;
    float overallScore;
    vector<float> tagNoise;
    
    //index
    Motif(int i):index(i),overallScore(0),score(0),motifProb(0)
                {query = translate(index);
                pwm[0].assign(K, 0);
                pwm[1].assign(K, 0);
                pwm[2].assign(K, 0);
                pwm[3].assign(K, 0);
                }
    
    void locateMotif(const char T[]);
    inline int mapLoci(int genomePos,const vector<int>& tag);
    void testMotifTag(genomeRegions &gR,bool ifDraw);
    void drawDist(genomeRegions &gR);
    void initBin(genomeRegions &gR);
    //order 0 order 1 
    void initProb(const genomeRegions& gR,int order);
    void initPWM();
    //void calPWM(const vector<Motif>& allmotifs);
    void calConscore(int nSize){
        float Thresh = motifProb*nSize;
        if (loci.size()>Thresh*DELTA){
            score = loci.size()/Thresh;
        }
        else {
            score = 0;
        }
    };
    float pvalue(){
        return 1-calPhi(score-1);
    }
    string antisense(const string& tempString);
    bool ascending(vector<char> &x,const vector<int> &traverseP,int &travIndex); 
    //write loci to file
    //score for each occurence
    vector<int> lociScore;
    bool writeLoci(ostream &s,genomeRegions &gR);
    void initLociScore();
    //inlines
    inline bool noWildcard();
    inline int wildcardNum();
    //if its counterpart's (e.g. A*TTA <=>TAA*T ) string is smaller then it's implicit
    inline bool implicit();
    inline bool isRepeat();
    void sumOverallScore();
    float sumTagScore();
    static int distance[4][15]; 
    //print info
    void inline printMotif();  
    inline static string translate(int i);
    //cluster method

};

class Cluster: public Motif{
public:
    Cluster(const Motif &m):Motif(m){index=-1;};
    inline void addProb(const Motif& m,int prevSize);
    void concatenate(const Motif& m,int index ,int optimShift);
    std::pair<int,int> editDistance(const Motif& m);
    float tagDistrDistance(const Motif& m);
    void calPWM(const Motif& m,int optimShift);
    inline void appendLoci(const Motif& m);
    void reCalSumBin(const Motif& m,const genomeRegions &gR);
    void mergeLoci();
    void trim();
    bool trivial(int pos);
    void writeCLusterLog(ostream &s,const Motif& m);
    //sumpwm for trim
    vector<int> totalPWM;
    int sumPWM();
};

//temp roc ploting
class tempLociScore: public Motif{
public:
    tempLociScore(const Motif &m,int num):Motif(m),counts(num),TP(false){index=-2;};
    void CountOnly(){myScore=counts;}
    // problematic
    void CountAndBipeak(){myScore=counts+LogOnePlusX(tagBiPeak[0]+tagBiPeak[1]+tagBiPeak[2]);}
    void CountAndIntensity(){myScore=counts+LogOnePlusX(signalIntensity[0]+signalIntensity[1]+signalIntensity[2]);}
    int counts;
    int myScore;
    bool TP;
};

class compareLoci {
public:
    bool operator()(tempLociScore lhs, tempLociScore rhs){
        return lhs.myScore>rhs.myScore;
    }
};




inline bool Motif::noWildcard(){
    for (int i=0; i<K; i++) 
        if (query[i]!='A'&&query[i]!='C'&&query[i]!='G'&&query[i]!='T') 
            return false;
    return true;
}

inline bool Motif::implicit(){
    //first it has a different counterpart
    if (::antisense(query)>=query) {
        return false;
    }
    else {
        return true;
    }
}
//remove TTTTTT TTATTT ATTTTT ATATAT
inline bool Motif::isRepeat(){
    int counter = 0;
    for (int i=1; i<K; i++) {
        if (query[i]!=query[0]) {
            counter++;
            if (i>=2&&query[i]!=query[1]&&query[1]!=query[0]) {
                return false;
            }
        }
    }
    //TTTTTT TTATTT
    if (counter==1||counter==0) {
        return true;
    }
    //ATTTTT
    bool flag=true;
    for (int i=2; i<K; i++) {
        if (query[i]!=query[1]) {
            flag=false;
        }
    }
    if (flag) {
        return true;
    }
    //ATATAT
    flag = true;
    for (int i=2; i<K; i+=2) {
        if (query[i]!=query[0]) {
            flag=false;
        }
    }
    for (int i=3; i<K; i+=2) {
        if (query[i]!=query[1]) {
            flag=false;
        }
    }
    if (flag) {
        return true;
    }
    return false;
}

inline int Motif::wildcardNum(){
    int count=0;
    for (int i=0; i<K; i++) 
        if (query[i]=='N') 
            count++;
    return count;
}
 

//m1>m2?
/*
class compareMotifHeap {
public:
    bool operator()(Motif lhs, Motif rhs){
        return lhs.overallScorerhs.overallScore;
    }
};
*/

class compareMotif {
public:
    bool operator()(Motif lhs, Motif rhs){
        return lhs.overallScore>rhs.overallScore;
    }
};



inline string Motif::translate(int query){
    string tempS("");
    char currentNum;
    
    for (int i=0; i<K; i++){
        currentNum = query%4;
        query = int(query/4);
        switch (currentNum) {
            case 0:
                tempS.push_back('A');
                break;
            case 1:
                tempS.push_back('C');
                break;
            case 2:
                tempS.push_back('G');
                break;
            case 3:
                tempS.push_back('T');
                break;
            default:
                break;
        }
    }
    return tempS;
}

void inline Motif::printMotif(){
    cout.precision(4);
    cout<<std::setw(25)<<setiosflags(std::ios::left)<<query<<"\t"<<pvalue()<<"\t"<<score<<"\t";
    for (int i=0; i<tagBiPeak.size(); i++) {
        cout<<tagBiPeak[i]<<"\t";
    }
    for (int i=0; i<tagSymmetry.size(); i++) {
        cout<<tagSymmetry[i]<<"\t";
    }
    if (option["FFT"]=="T") {
        for (int i=0; i<tagNoise.size(); i++) {
            cout<<tagNoise[i]<<"\t";
        }
    }
    for (int i=0; i<signalIntensity.size(); i++) {
        cout<<signalIntensity[i]<<"\t";
    }
    cout<<loci.size()<<endl;
}
//print pwm to stream
ostream &operator<<( ostream &s, Motif &motif );
//clustering
inline void Cluster::appendLoci(const Motif& m){
    //reCal insentity if tag
    if (option["mode"]=="tag") {
        for (int i=0; i<signalIntensity.size(); i++) {
            signalIntensity[i] = (signalIntensity[i]*loci.size()+m.signalIntensity[i]*m.loci.size())/(loci.size()+m.loci.size());
        }
    }
    loci.insert(loci.end(),m.loci.begin(),m.loci.end());
}
inline void Cluster::addProb(const Motif& m,int prevSize){
    motifProb = ((m.motifProb/K)*m.loci.size()+(motifProb/prevSize)*loci.size())/(m.loci.size()+loci.size())*query.size();
}


#endif
