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
//#include <stdio.h>
//#include <stdlib.h>
#include "common.h"
using std::cout;
using std::endl;
class Motif{
    public:
    string query;
    //char alp2num(query[i-1])[K];
    int index;
    vector<int> pwm[4];
    //bool overThresh;
    vector<int> loci;
    vector<float> signif;
    float motifProb;
    //each vector is a bin   sumBin[bin][tag]
    vector<float> sumBin[2*SAMPLESIZE+1];
    //vector<int> expMotifs;
    //std::vector<string> Clusterexp;
    //conserve score;
    float score;
    float overallScore;
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
    void testMotifTag(genomeRegions &gR,const string& outPutDir,bool ifDraw);
    void initBin(genomeRegions &gR);
    //order 0 order 1 
    void initProb(const genomeRegions& gR,int order);
    void initPWM();
    //void calPWM(const vector<Motif>& allmotifs);
    void calConscore(int nSize){
        float Thresh = motifProb*nSize*DELTA;
        if (loci.size()>Thresh){
            score = int(loci.size()*1000/Thresh)/1000.0;
        }
        else {
            score = 0;
        }
    };
    string antisense(const string& tempString);
    bool ascending(vector<char> &x,const vector<int> &traverseP,int &travIndex);   
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
};

class Cluster: public Motif{
public:
    Cluster(const Motif &m):Motif(m){};
    inline void addProb(const Motif& m,int prevSize);
    void concatenate(const Motif& m,int index ,int optimShift);
    std::pair<int,int> editDistance(const Motif& m);
    float histoneDistrDistance(const Motif& m);
    void calPWM(const Motif& m,int optimShift);
    inline void appendLoci(const Motif& m);
    void trim();
    bool unif(int pos);
    //sumpwm for trim
    vector<int> totalPWM;
    int sumPWM();
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

inline bool Motif::isRepeat(){
    int counter = 0;
    for (int i=1; i<K; i++) {
        if (query[i]!=query[0]) {
            counter++;
        }
        if (counter>=MAXREPEATCNT) {
            return false;
        }
    }
    return true;
}

inline int Motif::wildcardNum(){
    int count=0;
    for (int i=0; i<K; i++) 
        if (query[i]=='N') 
            count++;
    return count;
}
 

//m1>m2?
class compareMotif {
public:
    bool operator()(Motif m1, Motif m2){
        return m1.score>m2.score;
        //return m1.overallScore>m2.overallScore;
    }
};

class compareCluster {
public:
    bool operator()(Motif m1, Motif m2){
        return m1.overallScore>m2.overallScore;
    }
};
/*
inline bool compareMotif(Motif m1, Motif m2)
{   
    if (m1.signif!=m2.signif) {
        return m1.signif>m2.signif>0?true:false;
    }
	return (m1.score>m2.score);
}
*/

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
    cout<<query<<"\t"<<score<<"\t";
    for (int i=0; i<signif.size(); i++) {
        cout<<signif[i]<<"\t";
    }
    cout<<loci.size()<<endl;
}
//print pwm to stream
ostream &operator<<( ostream &s, const Motif &motif );
//clustering
inline void Cluster::appendLoci(const Motif& m){
    loci.insert(loci.end(),m.loci.begin(),m.loci.end());
    //cout<<"+"<<Motif::translate(currentMotif.expMotifs[i])<<allmotifs[currentMotif.expMotifs[i]].loci.size()<<"="<<loci.size()<<endl;
}
inline void Cluster::addProb(const Motif& m,int prevSize){
    motifProb = ((m.motifProb/K)*m.loci.size()+(motifProb/prevSize)*loci.size())/(m.loci.size()+loci.size())*query.size();
}



#endif
