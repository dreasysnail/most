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
#include <stdio.h>
#include <stdlib.h>
#include "common.h"
using std::cout;
using std::endl;
class Motif{
    public:
    int pwm[4][K];
    string query;
    char queryNum[K];
    int index;
    bool overThresh;
    vector<int> loci;
    vector<float> signif;
    float probThresh;
    vector<int> expMotifs;
    //conserve score;
    float score;
    float overallScore;
    //index
    Motif(int i):index(i),overallScore(0),score(0),probThresh(0),overThresh(false)
                {query = translate(index);
                 memset(pwm, 0, sizeof(int [4][K]));}
    void fillMotif(const vector<Motif>& allmotifs);
    vector<int> explainMotif();
    
    void locateMotif(const char T[]);
    void testMotifTag(const vector<int>& tag);
    //order 0 order 1 
    void initProb(const genomeRegions& gR,int order);
    void initPWM();
    void calPWM(const vector<Motif>& allmotifs);
    string antisense(const string& tempString);
    bool ascending(int x[],const vector<int> &traverseP,int &travIndex);
    inline bool noWildcard();
    inline int wildcardNum();
    //if its counterpart's (e.g. A*TTA <=>TAA*T ) string is smaller then it's implicit
    inline bool implicit();
    inline bool isRepeat();
    void sumScore();
    static int motif[K_5];
    //print info
    void inline printMotif();
    
    
    
    
    inline string translate(int i);
};



inline bool Motif::noWildcard(){
    for (int i=0; i<K; i++) 
        if (query[i]=='*') 
            return false;
    return true;
}

inline bool Motif::implicit(){
    //first it has a different counterpart
    if (query[K-1]=='*'||::antisense(query)>=query) {
        return false;
    }
    else {
        return true;
    }
}

inline bool Motif::isRepeat(){
    for (int i=1; i<K; i++) {
        if (query[i]!='*'&&query[i]!=query[0]) {
            return false;
        }
    }
    return true;
}

inline int Motif::wildcardNum(){
    int count=0;
    for (int i=0; i<K; i++) 
        if (query[i]=='*') 
            count++;
    return count;
}
 

//m1>m2?
class compareMotif {
public:
    bool operator()(Motif m1, Motif m2){
        //return m1.score>m2.score;
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
    string temp;
    char currentNum;
    
    for (int i=0; i<K; i++){
        currentNum = query%5;
        query = int(query/5);
        switch (currentNum) {
            case 0:
                temp.push_back('*');
                queryNum[i]=-1;
                break;
            case 1:
                temp.push_back('A');
                queryNum[i]=0;
                break;
            case 2:
                temp.push_back('C');
                queryNum[i]=1;
                break;
            case 3:
                temp.push_back('G');
                queryNum[i]=2;
                break;
            case 4:
                temp.push_back('T');
                queryNum[i]=3;
                break;
            default:
                queryNum[i]=-1;
                break;
        }

    }
    return temp;
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
