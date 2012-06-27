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
#include <string>
#include <map>

using std::string;
using std::vector;
using std::map;
using std::cerr;
using std::endl;
using std::ostream;
using std::pair;

void parseCommandLine(int argc,char** argv,map<string, string> &option);

//define motif length
extern int K;
#define DELTA 3
//#define PSEUDO 0.001
//#define MAXGENOME 2000000
#define SAMPLESIZE 75
//200M
const long int MAX_LENGTH = 3e6;
const long int HASH_TABLE_SIZE = 6e6;  //A prime roughly 10% larger
//tag counter parameters
const int dist=2;     //range               
const int offset=20;    //gap


//const float SINGIFTHRESH = 2;
//const float SCORETHRESH =5;

//cluster
const int CLUSTERMAX = 25;
const int MOTIFMAX = 150;
const int MAXDISTANCE = 12;
const int SHIFT = 3;
const int MAXCLUSTERSIZE = 20;
const int MAXREPEATCNT = 3;
const float MAXSIGNDIST = 3;
const int MINOVERALLSCORE = 0;

//cut tag extend bound
const int EXTENDBOUND = (2*SAMPLESIZE+2)*dist+SAMPLESIZE*(MAXCLUSTERSIZE+1)+offset*SAMPLESIZE+10;
#undef OUTFASTA
//display for node string

#undef display
//#define display




struct genomeRegion {
    int startP;
    int endP;
    string chr;
};

class genomeRegions{
public:
    
    genomeRegions(int num):extend(num),segmentStartPos(0){segmentStartPos.push_back(0);};
    vector <genomeRegion> genomes;
    //final seqs
    map<string,string> regionSeqs;
    //final tags
    map<string,vector <short int> > regionTags;
    //genome prob for trans and single site
    float prob[5][4];
    //total size of each genome(<200M).
    map <string,int> genomeLength;
    map <string,string> rawGenome;
    map <string,map<string,vector<short int> > > rawTag;
    vector<string> chromeNames;
    vector<string> tagName;

    int extend;
    void getSeq(const string& outPutDir);
    void inline appendSeq(ostream &outFile,int startP,int endP,string &currentChr);
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
    int segmentCount;
    //store subscript in genometag for starting pos of segment
    vector<int> segmentStartPos;
   
};

bool chrCompare(const string& chr1,const std::string &chr2);
/* dnaRegion comparison function */
inline bool compareGenome(genomeRegion gr1, genomeRegion gr2)
{   
    if (gr2.chr!=gr1.chr) {
        return chrCompare(gr1.chr,gr2.chr);;
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
        return s;
    }
    s<<"(";
    for (int i=0; i<v.size()-1; i++)
    {
        s<<v[i]<<",";
    }    
    s<<v.back()<<")";
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

template <typename T>
ostream &operator<<( ostream &s, const std::pair<T, T> &_p){
    s<<"("<<_p.first<<","<<_p.second<<")"<<endl;
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
        case 'M':
            return 4;
        case 'R':
            return 5;
        case 'W':
            return 6;
        case 'S':
            return 7;
        case 'Y':
            return 8;
        case 'K':
            return 9;
        case 'B':
            return 10;
        case 'D':
            return 11;
        case 'H':
            return 12;
        case 'V':
            return 13;
        case 'N':
            return 14;      
        default:
            return -1;
    }
}


float pow1(float base,int index);
string antisense(const string& tempString);
char degenerate(char a,char b);
//assume sorted , locate subscript return <sub,dist>

pair<int,int> locateSubscript(const vector<int> &listObj, vector<int>::const_iterator begin,vector<int>::const_iterator end, int queryVal);


#endif
