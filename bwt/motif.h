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
    float probThresh;
    vector<int> expMotifs;
    std::vector<string> Clusterexp;
    //conserve score;
    float score;
    float overallScore;
    //index
    Motif(int i):index(i),overallScore(0),score(0),probThresh(0)
                {query = translate(index);
                pwm[0].assign(K, 0);
                pwm[1].assign(K, 0);
                pwm[2].assign(K, 0);
                pwm[3].assign(K, 0);
                }
    void fillMotif(const vector<Motif>& allmotifs);
    vector<int> explainMotif();
    
    void locateMotif(const char T[]);
    void testMotifTag(genomeRegions &gR,const string& outPutDir);
    //order 0 order 1 
    void initProb(const genomeRegions& gR,int order);
    void initPWM();
    //void calPWM(const vector<Motif>& allmotifs);
    void calConscore(int nSize){
        float Thresh = probThresh*nSize*DELTA;
        if (loci.size()>Thresh){
            score = int(loci.size()*1000/Thresh)/1000.0;
        }
        else {
            score = 0;
        }
    };
    string antisense(const string& tempString);
    bool ascending(vector<char> &x,const vector<int> &traverseP,int &travIndex);
    //clustering methods
    std::pair<int,int> editDistance(const Motif& cluster);

    
    void concatenate(const Motif& m,int index ,int optimShift);
    void calPWM(const Motif& m,int optimShift);
    void trim();
    
    
    //inlines
    inline bool noWildcard();
    inline int wildcardNum();
    //if its counterpart's (e.g. A*TTA <=>TAA*T ) string is smaller then it's implicit
    inline bool implicit();
    inline bool isRepeat();
    void sumScore();
    static int motif[K_5];
    static int distance[4][15];

    
    
    //print info
    void inline printMotif();
    
    
    
    
    inline static string translate(int i);
};





inline bool Motif::noWildcard(){
    for (int i=0; i<K; i++) 
        if (query[i]=='N') 
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
        currentNum = query%5;
        query = int(query/5);
        switch (currentNum) {
            case 0:
                tempS.push_back('N');
                //alp2num(query[i-1])[i]=-1;
                break;
            case 1:
                tempS.push_back('A');
                //alp2num(query[i-1])[i]=0;
                break;
            case 2:
                tempS.push_back('C');
                //alp2num(query[i-1])[i]=1;
                break;
            case 3:
                tempS.push_back('G');
                //alp2num(query[i-1])[i]=2;
                break;
            case 4:
                tempS.push_back('T');
                //alp2num(query[i-1])[i]=3;
                break;
            default:
                //alp2num(query[i-1])[i]=-1;
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


#endif
