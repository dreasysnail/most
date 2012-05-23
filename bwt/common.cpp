//
//  common.cpp
//  bwt
//
//  Created by zhang yizhe on 12-5-20.
//  Copyright (c) 2012年 SJTU. All rights reserved.
//

#include "common.h"
#include "motif.h"


bool genomeRegions::readBed(const string &filename){
    ifstream bedFile(filename.c_str());
    if (!bedFile) {
        string errorInfo = "Error! Fail to open bed file \"" + filename + "\" !";
		printAndExit(errorInfo);
    }
    
    string line;
    genomeRegion segment;
    while (getline(bedFile, line)) {
        trim(line);
		if(line.size() == 0) continue;
		/* Read Bed line */
		istringstream ss(line);
		
        
		string temp;
		ss>>segment.chr;
		ss>>temp;
		segment.startP = atoi(temp.c_str())-extend;
		ss>>temp;
		segment.endP = atoi(temp.c_str())+extend;
        
		/* Store the bed line */
 		genomes.push_back(segment);
    }
    
	
	return true;
}


bool genomeRegions::readFasta(const string &filename){
    ifstream fastaFile(filename.c_str());
    if (!fastaFile) {
        string errorInfo = "Error! Fail to open fasta file \"" + filename + "\" !";
		printAndExit(errorInfo);
    }
    
    string chr("");
    string line;
    while(getline(fastaFile,line))
	{
		trim(line);
		if(line.size() == 0) continue;
		if(line[0] == '>')//Coordinates of the region
		{
            chr = line.substr(1,line.size()-1);
            if (rawGenome.count(chr)) {
                string errorInfo = "Error! chrome "+chr+" information collides\n";
				cout<<errorInfo;
                continue;
            }
            else {
                rawGenome[chr] = "";
                continue;
            }
        }
        //if seq line
        if (chr=="") {
            cout<<"Error! no chrome name?"<<endl;
            continue;
        }
        else {
            rawGenome[chr] += line;
        }
    }
    //sort genomes
    
    sort(genomes.begin(), genomes.end(), compareGenome);
    getSeq();
	return true;
}



void genomeRegions::getSeq(){

    
    
    //assert sorted genomes, deal with overlaps
    int startCut=0;
    int endCut=0;
    string currentChr("");
    vector<genomeRegion>::iterator it;
    try {
        for (it = genomes.begin(); it!=genomes.end(); it++) {
            if (it->chr!=currentChr) {
                genomeSeqs.push_back(rawGenome[currentChr].substr(startCut,endCut-startCut+1));
                startCut=it->startP;
                endCut=it->endP;
                //clear memory
                rawGenome[currentChr].clear();
                currentChr=it->chr;
            }
            else {
                if (it->startP>endCut) {
                    genomeSeqs.push_back(rawGenome[currentChr].substr(startCut,endCut-startCut+1));
                    startCut = it->startP;
                    endCut = it->endP;
                }
                else {
                    endCut=it->endP>endCut?it->endP:endCut;
                }
            }
        }
        //last one
        genomeSeqs.push_back(rawGenome[currentChr].substr(startCut,endCut-startCut+1));
        //cout<<genomeSeqs[1].size()<<endl; 
        //clear memory
        rawGenome[currentChr].clear();
    } catch (exception &e) {
        cerr<<e.what()<<endl;
        cerr<<"chrome:"<<currentChr<<" start:"<<startCut<<" end:"<<endCut<<" genomesize"<<rawGenome[currentChr].size()<<endl;
    }
    
}

void genomeRegions::writeRawTag(genomeRegions &tagBed){
    //sort tagBed
    sort(tagBed.genomes.begin(), tagBed.genomes.end(), compareGenome);
    string currentChr("");
    for (int i=0; i<tagBed.genomes.size(); i++) {
        if (tagBed.genomes[i].chr!=currentChr) {
            currentChr = tagBed.genomes[i].chr;
            vector<int> temp;
            if (rawTag.count(currentChr)) {
                string errorInfo = "Error! chrome "+currentChr+" information collides in Tag file\n";
				cout<<errorInfo;
                continue;
            }
            temp.assign(tagBed.genomes.back().endP+1, 0);
            rawTag[currentChr] = temp;

        }
        for (int j=tagBed.genomes[i].startP-1; j<tagBed.genomes[i].endP; j++) {
            rawTag[currentChr][j]++;
            //            cout<<j;
        }
    }
}


void genomeRegions::getTagBed(){
        
    //assert sorted genomes, deal with overlaps
    int startCut=0;
    int endCut=0;
    bool first = true;
    string currentChr("");
    vector<genomeRegion>::iterator it;
    try {
        for (it = genomes.begin(); it!=genomes.end(); it++) {
            if (it->chr!=currentChr) {
                if (first) {
                    first = false;
                }
                else {
                    appendTag(startCut, endCut, currentChr);
                }
                startCut=it->startP;
                endCut=it->endP;
                //clear memory
                rawGenome[currentChr].clear();
                currentChr=it->chr;
            }
            else {
                if (it->startP>endCut) {
                    appendTag(startCut, endCut, currentChr);            
                    startCut = it->startP;
                    endCut = it->endP;
                }
                else {
                    endCut=it->endP>endCut?it->endP:endCut;
                }
            }
        }
        //last one
        appendTag(startCut, endCut, currentChr);
        //clear memory
        rawTag[currentChr].clear();
    } catch (exception &e) {
        cerr<<e.what()<<endl;
        cerr<<"chrome:"<<currentChr<<" start:"<<startCut<<" end:"<<endCut<<" genomesize"<<rawTag[currentChr].size()<<endl;
    }
    
}

void genomeRegions::appendTag(int a,int b,const string& chr){
    int tempsize = rawTag[chr].size();
    //    cout<<tempsize<<endl;
    for (int i=a-1; i<b; i++) {
        if (i>=tempsize) {
            genomeTags.push_back(0);
            //            cout<<i;
        }
        else {
            genomeTags.push_back(rawTag[chr][i]);
            //    cout<<endl<<i;
        }
    }
    genomeTags.push_back(2);
}

void printProgress(const int i,const string& message){
    if (i==0) {
        cout<<endl<<message<<endl;
    }
    if (i%(int((K_5)/10)+1)==0) {
        cout<<">>>";
    }
}




