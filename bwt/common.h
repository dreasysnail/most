//
//  common.h
//  bwt
//
//  Created by zhang yizhe on 12-5-20.
//  Copyright (c) 2012年 SJTU. All rights reserved.
//

#ifndef bwt_common_h
#define bwt_common_h

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <map>


//define motif length
#define K 9
#define K_5 1953125
//#define K 5
//#define K_5 3125
#define DELTA 4
//#define PSEUDO 0.001
//#define MAXGENOME 2000000
#define SAMPLESIZE 5
//200M
const long int MAX_LENGTH = 1e6;
const long int HASH_TABLE_SIZE = 1e6;  //A prime roughly 10% larger
//tag counter parameters
const int dist=100;     //range               
const int offset1=100;    //gap
const int offset=20;    //gap
const float SINGIFTHRESH = 2;
const float SCORETHRESH =5;



#undef OUTFASTA
#define DRAW_MOTIF
//display for node string

#undef display
//#define display


using std::string;
using std::vector;
using std::map;
using std::cerr;
using std::endl;
using std::ostream;

struct genomeRegion {
    int startP;
    int endP;
    string chr;
};

class genomeRegions{
public:
    
    genomeRegions(int num):extend(num){};
    vector <genomeRegion> genomes;
    //final seqs
    map<string,string> genomeSeqs;
    //final tags
    map<string,vector <int> > genomeTags;
    //genome prob for trans and single site
    float prob[5][4];
    //total size of each genome(<200M).
    map <string,int> genomeLength;
    map <string,string> rawGenome;
    map <string,map<string,vector<int> > > rawTag;
    vector<string> chromeNames;
    vector<string> tagName;

    int extend;
    void getSeq(const string& outPutDir);
    void inline writeSeq(ostream &outFile,int startP,int endP,string &currentChr);
    int catSeq(char* T);
    void getTagBed(const string& thisHistone,const string& currentChr );    //get 0-1 sequence from bed file
    void appendTag(int a,int b,const string& chr,const string& thisHistone);
    void appendReverse();
    void initProb(int mode);
    bool existChr(const string& chr){return find(chromeNames.begin(),chromeNames.end(),chr)==chromeNames.end()?false:true;}
    void printProb();
    void writeRawTag(genomeRegions &tagBed);
    void mergeOverlap();
    bool readBed(const string &filename);
    //read and write rawTag;
    bool readWig(const string &filename);
    
    bool readFasta(const string &filename);
   
};

/* dnaRegion comparison function */
inline bool compareGenome(genomeRegion gr1, genomeRegion gr2)
{   
    int isFirstBigger = gr1.chr.compare(gr2.chr);
    if (isFirstBigger!=0) {
        return isFirstBigger>0?false:true;
    }
	return (gr1.startP<gr2.startP);
}


inline void printAndExit(string errorInfo)
{
	cerr<<errorInfo<<endl;
	exit(1);
}

inline string trim(std::string& str)
{
	str.erase(0, str.find_first_not_of(' '));       //prefixing spaces
	str.erase(str.find_last_not_of(' ')+1);         //surfixing spaces
	return str;
}
inline string trimN(string &temp){
    string result("");
    for (int i=0; i<temp.size(); i++) {
        if (temp[i]=='N') {
            continue;
        }
        result.push_back(temp[i]);
    }
    return result;
}

void printProgress(const int i,const int total,const string& message);

template <class T>
ostream &operator<<( ostream &s, const vector<T> &v){
    if (v.size()==0){
        s<<"NULL vector"<<endl;
        return s;
    }
    s<<"(";
    for (int i=0; i<v.size()-1; i++)
    {
        s<<v[i]<<",";
    }    
    s<<v.back()<<")"<<endl;
    return s;
};

template <class T>
ostream &operator<<( ostream &s, const vector<vector<T> > &m){
    for (int i=0; i<m[0].size(); i++) {
        for (int j=0; j<m.size(); j++) {
            s<<m[i][j]<<'\t';
        }
        s<<endl;
    }
    return s;
};

inline int alp2num(const char& name){
    switch (name) {
        case 'A':
            return 0;
        case 'C':
            return 1;
        case 'G':
            return 2;
        case 'T':
            return 3;
        default:
            return -1;
    }
}
float pow1(float base,int index);
string antisense(const string& tempString);
#endif
