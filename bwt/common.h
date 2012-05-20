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
    vector <genomeRegion> genomes;
    vector <string> genomeSeqs;
    map <string,string> rawGenome;
    void getSeq();
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

#endif
