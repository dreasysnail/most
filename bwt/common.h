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
#include <map>

using namespace std;
struct genomeRegion {
    int startP;
    int endP;
    string chr;
};

class genomeRegions{
public:
    
    genomeRegions(int num):extend(num){};
    vector <genomeRegion> genomes;
    vector <string> genomeSeqs;
    vector <int> genomeTags;
    map <string,string> rawGenome;
    map <string,vector<int> > rawTag;
    int extend;
    void getSeq();
    void getTagBed();    //get 0-1 sequence from bed file
    void appendTag(int a,int b,const string& chr);
    void writeRawTag(genomeRegions &tagBed);
    bool readBed(const string &filename);
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

void printProgress(const int i,const string& message);

template <class T>
ostream &operator<<( ostream &s, const vector<T> &v){
    for (int i=0; i<v.size(); i++)
    {
        s<<v[i];
    }    
    s<<endl;
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
#endif
