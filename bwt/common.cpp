//
//  common.cpp
//  bwt
//
//  Created by zhang yizhe on 12-5-20.
//  Copyright (c) 2012年 SJTU. All rights reserved.
//

#include "common.h"
#include "motif.h"
#include <cassert>
using namespace std;


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
        if (!existChr(segment.chr)) continue;
		ss>>temp;
		segment.startP = atoi(temp.c_str())-extend;
		ss>>temp;
		segment.endP = atoi(temp.c_str())+extend;
        if (segment.startP<0||segment.startP>segment.endP) continue;
		/* Store the bed line */
 		genomes.push_back(segment);
    }
    
	
	return true;
}
//simple function ,less api
bool genomeRegions::readWig(const string &filename){
    ifstream wigFile(filename.c_str());
    //assume sorted wig file
    if (!wigFile) {
        string errorInfo = "Error! Fail to open wig file \"" + filename + "\" !";
		printAndExit(errorInfo);
    }

    string line;
    genomeRegion segment;
    int wigStep;
    string currentChr("");
    string thisHistone("");
        
    while (getline(wigFile, line)) {
        //new histone/tag
        if (line[0]=='v') {
            //find info
            size_t found = line.find("histone=");
            if (found!=string::npos){
                size_t found2 = line.find("span=");
                if (found2 == line.npos) {
                    cerr<<"FATAL: Wigfile's info line contains no \"span=\""<<endl;
                    exit(1);
                }
                tagName.push_back(line.substr(found+8,found2-found-9));
                cerr<<"find Tag "<<tagName.back()<<endl;
            }
            else {
                char temptagName[20];
                sprintf (temptagName, "TagNo.%lu", tagName.size()+1);
                tagName.push_back(temptagName);
                cerr<<"find Tag "<<tagName.back()<<" You'd better contain histone=**"<<endl;
            }
            //check histone contains all chrome info that is needed
            
             
            if (thisHistone!="") {
                for (vector<string>::iterator fastaIt=chromeNames.begin(); fastaIt!=chromeNames.end(); fastaIt++) {
                    if (rawTag[thisHistone].count(*fastaIt)<=0) {
                        cerr<<"WARNING:"<<thisHistone<<" in TagFile doesn't contain chrome "<<(*fastaIt)<<" appears in fastaFile"<<endl;
                        rawTag[thisHistone][*fastaIt].assign(genomeLength[*fastaIt], 0);
                        getTagBed(thisHistone,*fastaIt);
                    }
                }
                getTagBed(thisHistone,currentChr);
            }

            
            thisHistone = tagName.back();
            
            found = line.find("span=");
            if (found!=string::npos){
                wigStep=atoi(line.substr(found+5,line.size()-found-5).c_str());
            }
            else {
                cerr<<"check wig file's info line: should contain span=**"<<endl;
                exit(0);
            }
            currentChr = "";
            continue;
        }
		if(line.size() == 0) continue;
		/* Read wig line */
		istringstream ss(line);
		
        
		string temp;
		ss>>segment.chr;
        if (!existChr(segment.chr)) continue;
        ss>>temp;
		segment.startP = atoi(temp.c_str());
        if (segment.startP<0) continue;
		segment.endP = segment.startP + wigStep;
        ss>>temp;
        int score = atof(temp.c_str())*4;
        if (segment.chr!=currentChr) {
            if (currentChr!="") {
                getTagBed(thisHistone,currentChr);
            }
            currentChr = segment.chr;
            //or find  != .end()
            if (rawTag[thisHistone].count(currentChr)) {
                string errorInfo = "Error! chrome "+currentChr+" information collides in Tag file\n";
				cout<<errorInfo;
                continue;
            }
            //cout<<line<<currentChr<<endl;
            

            rawTag[thisHistone][currentChr].assign(genomeLength[currentChr], 0);            
        }
        if (score<=0) {
            continue;
        }
        for (int j=segment.startP-1; j<segment.endP; j++) {
            rawTag[thisHistone][currentChr][j]+=score;
            //            cout<<j;
        }
        
    }
    //last histone's last chrome
    getTagBed(thisHistone,currentChr);
    //append reverse for antisense,for each histone
    appendReverse();
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
                chromeNames.push_back(chr);
                continue;
            }
        }
        //if seq line
        if (chr=="") {
            cout<<"Error! no chrome name?"<<endl;
            continue;
        }
        else {
            transform(line.begin(), line.end(), line.begin(), ::toupper);
            rawGenome[chr] += line;
        }
    }
    
    
//  initProb(1);
	return true;
}



void genomeRegions::getSeq(const string& outPutDir){
    string fileName = outPutDir+"/allregions.fa";
    ofstream regionFile(fileName.c_str());
    vector<genomeRegion>::iterator it;
    
    for (it = genomes.begin(); it!=genomes.end(); it++) {
        writeSeq(regionFile,it->startP,it->endP,it->chr);
    }
    //clear memory
    for (vector<string>::iterator it=chromeNames.begin(); it!=chromeNames.end(); it++) {
        genomeLength[*it] = rawGenome[*it].size();
        string tempS("");
        rawGenome[*it].swap(tempS);
    }
    
    //output is also sorted
    
}


void inline genomeRegions::writeSeq(ostream &outFile,int startCut,int endCut,string &currentChr){
    try {
        if (endCut<=startCut+1) {
            return;
        }
        genomeSeqs[currentChr] += rawGenome[currentChr].substr(startCut,endCut-startCut+1)+"#";
        //generate fasta
    #ifdef OUTFASTA
        outFile<<">"<<currentChr<<":"<<startCut<<"-"<<endCut<<"\n";
        for (int i = startCut; i<endCut; i+=50) {
            outFile<<rawGenome[currentChr].substr(i,50)<<"\n";
            if (i+50 > endCut) {
                outFile<<rawGenome[currentChr].substr(i,endCut-i)<<"\n";
            }
        }
        outFile<<endl;
    #endif
    } catch (exception &e) {
        cerr<<e.what()<<endl;
        cerr<<"Seq insert error: chrome:"<<currentChr<<" start:"<<startCut<<" end:"<<endCut<<" chromesize"<<rawGenome[currentChr].size()<<endl;
        return;
    }

}

int genomeRegions::catSeq(char* T){
    map<string, string>::iterator chrit;
    long int offset=0;
    string tempString("");
    for (chrit=genomeSeqs.begin(); chrit!=genomeSeqs.end(); chrit++) {
        tempString += chrit->second;
        string tempS("");
        chrit->second.swap(tempS);
    }
    tempString += antisense(tempString)+"#";
    offset = tempString.size();
    cerr<<"regionSize:"<<offset<<endl;
    assert(offset<=MAX_LENGTH);
    strcpy(T, tempString.c_str()); 
    tempString.clear();
    return offset;
}



void genomeRegions::writeRawTag(genomeRegions &tagBed){
    //sort tagBed
    sort(tagBed.genomes.begin(), tagBed.genomes.end(), compareGenome);
    string currentChr("");
    string thisHistone = "bed";
    map<string, vector<int> > tempMap;
    rawTag[thisHistone] = tempMap;
    for (int i=0; i<tagBed.genomes.size(); i++) {
        if (tagBed.genomes[i].chr!=currentChr) {
            currentChr = tagBed.genomes[i].chr;
            vector<int> temp;
            if (rawTag[thisHistone].count(currentChr)) {
                string errorInfo = "Error! chrome "+currentChr+" information collides in Tag file\n";
                cout<<errorInfo;
                continue;
            }
            int lastThisChr=tagBed.genomes.size()-1;
            for (int k=i; k<tagBed.genomes.size(); k++) {
                if (tagBed.genomes[k].chr!=currentChr) {
                    lastThisChr = k-1;
                    break;
                }
            }
            temp.assign(tagBed.genomes[lastThisChr].endP+1, 0);
            rawTag[thisHistone][currentChr] = temp;

        }
        for (int j=tagBed.genomes[i].startP-1; j<tagBed.genomes[i].endP; j++) {
            rawTag[thisHistone][currentChr][j]++;
            //            cout<<j;
        }
    }

}


void genomeRegions::mergeOverlap(){
    sort(genomes.begin(), genomes.end(), compareGenome);
    //assert sorted genomes, deal with overlaps
    vector<genomeRegion> mergedGenome; 
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
                    genomeRegion tempG;
                    tempG.startP = startCut;
                    tempG.endP = endCut;
                    tempG.chr = currentChr;
                    mergedGenome.push_back(tempG);
                }
                startCut=it->startP;
                endCut=it->endP;
                //clear memory
                rawTag[currentChr].clear();
                currentChr=it->chr;
            }
            else {
                if (it->startP>endCut) {
                    genomeRegion tempG;
                    tempG.startP = startCut;
                    tempG.endP = endCut;
                    tempG.chr = currentChr;
                    mergedGenome.push_back(tempG);           
                    startCut = it->startP;
                    endCut = it->endP;
                }
                else {
                    endCut=it->endP>endCut?it->endP:endCut;
                }
            }
        }
        //last one
        genomeRegion tempG;
        tempG.startP = startCut;
        tempG.endP = endCut;
        tempG.chr = currentChr;
        mergedGenome.push_back(tempG);
        //clear memory
        rawTag[currentChr].clear();
        //append reverse for antisense
        appendReverse();
    } catch (exception &e) {
        cerr<<e.what()<<endl;
        cerr<<"chrome:"<<currentChr<<" start:"<<startCut<<" end:"<<endCut<<" genomesize"<<rawTag[currentChr].size()<<endl;
    }
    genomes.swap(mergedGenome);
}

void genomeRegions::getTagBed(const string& thisHistone,const string& currentChr){   
    vector<genomeRegion>::iterator it;
    bool rear = false;
    for (it = genomes.begin(); it!=genomes.end(); it++) {
        //cout<<thisHistone<<it->chr<<it->startP<<it->endP<<endl;
        if (it->chr!=currentChr) {
            if (rear) 
                break;
            else 
                continue;
        }
        else {
            appendTag(it->startP,it->endP,currentChr,thisHistone);
            rear=true;

        }
    }
    //clear memory
    vector<int> tempV;
    rawTag[thisHistone][currentChr].swap(tempV);
}

void genomeRegions::appendTag(int a,int b,const string& chr,const string& thisHistone){
    try {
        for (int i=a-1; i<b; i++) {
            if (i>=genomeLength[chr]-1) {
                break;
            }
            else {
                if (i>=rawTag[thisHistone][chr].size()) {
                    genomeTags[thisHistone].push_back(0);
                }
                else {
                    genomeTags[thisHistone].push_back(rawTag[thisHistone][chr][i]);
                }
            }
        }
        genomeTags[thisHistone].push_back(0);
    } catch (exception &e) {
        cerr<<e.what()<<endl;
        cerr<<"Tag insert error: chrome:"<<chr<<" start:"<<a<<" end:"<<b<<" chromesize"<<rawTag[thisHistone][chr].size()<<endl;
    }
}

void genomeRegions::appendReverse(){
    for (vector<string>::iterator tagIt=tagName.begin(); tagIt!=tagName.end(); tagIt++) {
        string thisHistone = *tagIt;
        vector<int> temp(genomeTags[thisHistone]);
        reverse(temp.begin(), temp.end());
        copy(temp.begin(), temp.end(),back_inserter(genomeTags[thisHistone]));
        genomeTags[thisHistone].push_back(0);
    }
}

void genomeRegions::initProb(int mode){
    long int countNum[5][4]={0};
    //whole genome
    if (mode==1) {
        map <string,string>::iterator it;
        for (it=rawGenome.begin(); it!=rawGenome.end(); it++) {
            for (int i=0; i<it->second.size()-1; i++ ) {
                if (alp2num(it->second[i])==-1||alp2num(it->second[i+1])==-1) {
                    continue;
                }
                countNum[0][alp2num(it->second[i])]++;
                countNum[alp2num(it->second[i])+1][alp2num(it->second[i+1])]++;
            }
        }
    }
    //only selected region
    else if (mode==2) {
        map<string,string>::iterator it;
        for (it=genomeSeqs.begin(); it!=genomeSeqs.end(); it++ ) {
            if ((it->second).size()==0) {
                continue;
            }
            for (int i=0; i<(it->second).size()-1; i++ )  {
                
                if (alp2num((it->second)[i])==-1||alp2num((it->second)[i+1])==-1) {
                    continue;
                }
                countNum[0][alp2num((it->second)[i])]++;
                countNum[alp2num((it->second)[i])+1][alp2num((it->second)[i+1])]++;
            }
        }

    }
   
    for (int i=0; i<5; i++) {
        for (int j=0; j<4; j++) {
            prob[i][j]=countNum[i][j]/float(countNum[i][0]+countNum[i][1]+countNum[i][2]+countNum[i][3]);
        }
    }
    printProb();
}

void genomeRegions::printProb(){
    cerr<<"prob matrix:"<<"\n"<<"A\t\tC\t\tG\t\tT\n";
    for (int i=0; i<5; i++) {
        for (int j=0 ; j<4 ; j++) {
            cerr<<int(prob[i][j]*1000)/1000.0<<"\t";
        }
        cerr<<"\n";
    }
    cerr<<endl;
}

void printProgress(const int i,const int total,const string& message){
    if (i==0) {
        cerr<<message<<endl;
    }
    if (i%(int((total)/50)+1)==0) {
        cerr<<">";
    }
}

float pow1(float base,int index){
    float pow = 1;
    while (index>0) {
        if (index%2==1) {
            pow = pow * base;
        }
        base = base * base;
        index = int(index/2);
    }
    return pow;
}

string antisense(const string& tempString){
    string output("");
    for (int i=tempString.size()-1; i>=0; i--) {
        switch (tempString[i]) {
            case 'A':
                output.push_back('T');
                break;
            case 'T':
                output.push_back('A');
                break;
            case 'C':
                output.push_back('G');
                break;
            case 'G':
                output.push_back('C');
                break;
            case '#':
                output.push_back('#');
                break;
            case 'N':
                output.push_back('N');
                break;
            default:
                break;
        }
    }
    return output;
}


