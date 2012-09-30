//
//  motif.h
//  MOST
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
#include "bwt.h"
using std::cout;
using std::endl;
using std::string;
using std::vector;
class Suffix;

class Motif{
    public:
    string query;
    //char alp2num(query[i-1])[K];
    int index;
    vector<int> pfm[4];
    //bool overThresh;
    vector<int> loci;
    vector<float> tagBiPeak;
    vector<float> tagSymmetry;
    float motifProb;
    //each vector is a bin   sumBin[bin][tag]
    vector<float> sumBin[2*SAMPLESIZE+1];
    //vector<float> signalIntensity;
    //vector<int> expMotifs;
    //std::vector<string> Clusterexp;
    //conserve score;
    float score;
    float overallScore;
    vector<float> tagNoise;
    
    //index
    Motif(int i):index(i),overallScore(0),score(0),motifProb(0)
                {query = index2Str(index);
                pfm[0].assign(K, PEUSUDO_COUNT);
                pfm[1].assign(K, PEUSUDO_COUNT);
                pfm[2].assign(K, PEUSUDO_COUNT);
                pfm[3].assign(K, PEUSUDO_COUNT);
                }
    Motif(string i):query(i),overallScore(0),score(0),motifProb(0)
                {index = str2index(query);
                pfm[0].assign(K, PEUSUDO_COUNT);
                pfm[1].assign(K, PEUSUDO_COUNT);
                pfm[2].assign(K, PEUSUDO_COUNT);
                pfm[3].assign(K, PEUSUDO_COUNT);
                }
    
    void locateMotif(const char T[]);
    inline int mapLoci(int genomePos,const vector<int>& tag);
    void CalSumBin(genomeRegions &gR);
    void testMotifTag(genomeRegions &gR,bool ifDraw);
    void drawDist(genomeRegions& gR,ostream &distFile);
    void initBin(genomeRegions &gR);
    //order 0 order 1 
    void initProb(const genomeRegions& gR,int order);
    void initPFM();
    //void calPFM(const vector<Motif>& allmotifs);
    void calConscore(int nSize){
        //protocol: motifProb loci
        if (loci.size()>motifProb*nSize*DELTA)
            score = loci.size()/motifProb/nSize;
        else 
            score = 0;
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
    void printMotif(ostream &s);  
    inline string index2Str(int i);
    inline int str2index(string &i);
    //pfm
    float PearsonCorrPFM(const Motif &m,int offset,int clustersize,bool strand);
    void generateIUPAC();
    //sumpfm for trim
    vector<int> totalPFM;
    int sumPFM();


};

struct Member{
    int shift;
    std::vector<float> tagBiPeak;
    std::vector<float> tagSymmetry;
    float dist;
    string query;
};

class Cluster: public Motif{
public:
    Cluster(const Motif &m):Motif(m){index=-1;
                                     sumPFM();
                                     Member tempM ={0,m.tagBiPeak,m.tagSymmetry,0,m.query};
                                     Members.push_back(tempM);};
    inline void addProb(const Motif& m,int prevSize);
    void concatenate(const Motif& m,int index ,int optimShift);
    std::pair<float,int> editDistance(const Motif& m);
    float tagDistrDistance(const Motif& m);
    void calPFM(const Motif& m,int optimShift);
    inline void appendLoci(const Motif& m);
    void reCalSumBin(const Motif& m,const genomeRegions &gR);
    void mergeLoci();
    void mergeProb(const Motif& m);
    void trim();
    bool trivial(int pos);
    bool oligo(int pos);
    void writeCLusterLog(ostream &s,const Motif& m);
    //get extension from word m
    void getExtended(const Motif &m, genomeRegions &gR, Suffix & active);
    //motif member
    vector<Member> Members;
    void printMember(ostream &s);
};

//temp roc ploting
class tempLociScore: public Motif{
public:
    tempLociScore(const Motif &m,int num):Motif(m),counts(num),TP(false){index=-2;};
    void CountOnly(){myScore=counts;}
    // problematic
    void CountAndBipeak(){
        myScore=counts+LogOnePlusX(tagBiPeak[0]+tagBiPeak[1]+tagBiPeak[2]);

        //myScore=counts+(tagBiPeak[0]+tagBiPeak[1]+tagBiPeak[2])/600;
    }
    //void CountAndIntensity(){myScore=counts+(signalIntensity[0]+signalIntensity[1]+signalIntensity[2])/30000;}

    int counts;
    float myScore;
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
//remove simple tandem repeat TTTTTT TTATTT ATTTTT ATATAT
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
    if (option["rmrepeat"]=="2") {
        char tempChar=query[0];
        for (int i=1; i<K; i++) {
            if (tempChar!=query[0]&&query[i]!=tempChar&&query[i]!=query[0]) {
                return false;
            }
            if (query[i]!=query[0]) {
                tempChar=query[i];
            }
        }
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
 


class compareMotif {
public:
    bool operator()(Motif lhs, Motif rhs){
        return lhs.overallScore>rhs.overallScore;
    }
};



inline string Motif::index2Str(int query){
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

inline int Motif::str2index(string &query){
    int tempIndex = 0;
    for (int i=K-1; i>=0; i--){
        tempIndex = tempIndex*4 + alp2num(query[i]);
        
    }
    return tempIndex;
}


//print pfm to stream
ostream &operator<<( ostream &s, Motif &motif );
//clustering
inline void Cluster::appendLoci(const Motif& m){
    loci.insert(loci.end(),m.loci.begin(),m.loci.end());
}
inline void Cluster::addProb(const Motif& m,int prevSize){
    //motifProb = ((m.motifProb/K)*m.loci.size()+(motifProb/prevSize)*loci.size())/(m.loci.size()+loci.size())*query.size();
    motifProb = m.motifProb+motifProb;
}





#endif


