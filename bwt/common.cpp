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
#include <stdio.h>
#include <stdlib.h>
#include "FFT.h"
using namespace std;

int K;
map<string,string> option;


/***********************************************************************
 * parseCommandLine
 *
 *
 ***********************************************************************/
void parseCommandLine(int argc,
                      char** argv,
                      map<string, string> &option
                      ) {
    
    // Set default values for command line arguments
    cerr<<"Commandline:";
    
    map<string, bool> optionRequire;
    //t:tag n:normal
    option          ["mode"]        =   "normal";
    optionRequire   ["mode"]        =   true;
    //
    option          ["fastafile"]   =   "";
    optionRequire   ["fastafile"]   =   true;
    //
    option          ["regionfile"]  =   "";
    optionRequire   ["regionfile"]  =   true;
    //
    option          ["tagfile"]     =   "";
    optionRequire   ["tagfile"]     =   false;
    //
    option          ["outdir"]      =   "mosh_out";
    optionRequire   ["outdir"]      =   false;
    //
    option          ["k"]      =   "9";
    optionRequire   ["k"]      =   false;
    //
    option          ["order"]      =   "0";
    optionRequire   ["order"]      =   false;
    //
    option          ["writeloci"]      =  "T";
    optionRequire   ["writeloci"]      =   false;
    //
    option          ["regionfile"]      =  "F";
    optionRequire   ["regionfile"]      =   false;
    //
    option          ["bkgregion"]      =  "region";
    optionRequire   ["bkgregion"]      =   false;
    // Parse the command line.
    string option_name = "";
    string option_value = "";
    
    
    // Read the next option, and break if we're done.
    for (int i=1; i<argc-1; i+=2) {
        option_name = argv[i];
        option_value = argv[i+1];
        cerr<<"("<<option_name<<" "<<option_value<<") ";
        // Assign the parsed value to the appropriate variable
        if (option_name == "-m") {
            option["mode"] = option_value;
            optionRequire["mode"] = false;
            if (option_value=="tag") {
                optionRequire["tagfile"] = true;
            }
            else if (option_value=="normal") {
                
            }
            else if (option_value=="help") {
                
            }
            else {
                cerr<<"wrong mode!"<<endl;
                exit(1);
            }
        } 
        else if (option_name == "-f") {
            option["fastafile"] = option_value;
            optionRequire["fastafile"] = false;
        } 
        else if (option_name == "-r") {
            option["regionfile"] = option_value;
            optionRequire["regionfile"] = false;
        } 
        else if (option_name == "-t") {
            option["tagfile"] = option_value;
            optionRequire["tagfile"] = false;
        } 
        else if (option_name == "-o") {
            option["outdir"] = option_value;
            optionRequire["outdir"] = false;
        } 
        else if (option_name == "-k") {
            option["k"] = option_value;
        } 
        else if (option_name == "-bo") {
            option["order"] = option_value;
        }
        else if (option_name == "-locifile") {
            option["writeloci"] = option_name;
        }
        else if (option_name == "-regionfile") {
            option["regionfile"] = option_name;
        }
        else if (option_name == "-br") {
            option["bkgregion"] = option_name;
        }
        else {
            cerr<<"unrecognized option "<<option_name<<"! skip"<<endl;
            continue;
        } 
    }
    
    // verify required arguments.
    for (map<string, bool>::iterator it=optionRequire.begin();it!=optionRequire.end();it++) {
        if (it->second) {
            cerr<<"option "<<it->first<<" need to be specified in command line"<<endl;
            printUsage();
            exit(1);
        }
    }
    //changing parameters
    K = atoi(option["k"].c_str());
    if (K<1||K>15) {
        cerr<<"K should be integer among [1,15](default 9)"<<endl;
        exit(1);
    }
    int order = atoi(option["order"].c_str());
    if (order!=-1&&order!=0&&order!=1) {
        cerr<<"order should be -1,0,1(default 0)"<<endl;
        exit(1);
    }
    cerr<<endl;
    
}

void printUsage()
{
	string usage =	"*------------------------------------------*\n";
	usage		+=	"        Mosh " + MOSHVERSION+"\n";
	usage		+=	"        Developed by Yizhe Zhang, Chaochun Wei\n";
	usage		+=	"        jeremy071242044@gmail.com\n";
	usage		+=	"*------------------------------------------*\n\n";
	usage		+=	"    mosh [option: <parameter1> <parameter2>]  \n\n";
	usage		+=	"    Required Parameters:\n\n";
	usage		+=	"    -m <normal/tag>\tRuning mode:tagfile should be specified if tag mode is selected\n\n";
	usage		+=	"    -f <DNA sequence file>\tWhole genome DNA sequences\n\n";
    usage		+=	"    -b <BED file>\tregions of interest\n\n";
    usage		+=	"    -t <WIG file>\tTag file for Histone marks or other sources\n\n";
    usage		+=	"    CAVEAT:WIG FILE SHOULD BE SORTED\n\n";
	usage		+=	"    -h help\tPrint help information.\n\n";
    usage		+=	"    Optional Parameters:\n\n";
    usage		+=	"    -o <outputDIR>\tSpecify an output directory.\n\n";
    usage		+=	"    -bo <0,1> \tOrder for background sequence\n\n";
    usage		+=	"    -br <region/genome>\tspecify background sequence\n\n";
    usage		+=	"    -locifile <T/F>\tWhether or not loci file of each cluster should be exported\n\n";
    usage		+=	"    -regionfile <T/F>\tWhether or not loci file of each cluster should be exported\n\n";

    
    
	cerr<<usage<<endl;
}


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
    //already have bed and fasta, slice and paste tag for each chrome for each histone
    ifstream wigFile(filename.c_str());
    //assume sorted wig file
    if (!wigFile) {
        string errorInfo = "Error! Fail to open wig file \"" + filename + "\" !";
		printAndExit(errorInfo);
    }

    string line;
    genomeRegion segment;
    int wigStep=25;
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
            
            
            //check else and put last histone's last chr 
            if (thisHistone!="") {
                getTagBed(thisHistone,currentChr);
                //check histone contains all chrome info that is needed
                for (vector<string>::iterator fastaIt=chromeNames.begin(); fastaIt!=chromeNames.end(); fastaIt++) {
                    
                    if (rawTag[thisHistone].find(*fastaIt)==rawTag[thisHistone].end()) {
                        cerr<<"WARNING:"<<thisHistone<<" in TagFile doesn't contain chrome "<<(*fastaIt)<<" appears in fastaFile"<<endl;
                        rawTag[thisHistone][*fastaIt].assign(genomeLength[*fastaIt], 0);
                        getTagBed(thisHistone,*fastaIt);
                    }
                }
                
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
        //find new chr  push old and detect lack
        if (segment.chr!=currentChr) {
            //push old
            if (currentChr!="") {
                //put last one
                getTagBed(thisHistone,currentChr);
                //check histone contains all chrome info that is needed
                vector<string>::iterator from,to;
                from  = find (chromeNames.begin(), chromeNames.end(), currentChr)+1;
                to  = find (chromeNames.begin(), chromeNames.end(), segment.chr);
                if (from!=to) {
                    for (vector<string>::iterator fastaIt=from; fastaIt!=to; fastaIt++) {
                        if (rawTag[thisHistone].find(*fastaIt)==rawTag[thisHistone].end()) {
                            cerr<<"WARNING:"<<thisHistone<<" in TagFile doesn't contain chrome "<<(*fastaIt)<<" appears in fastaFile"<<endl;
                            rawTag[thisHistone][*fastaIt].assign(genomeLength[*fastaIt], 0);
                            getTagBed(thisHistone,*fastaIt);
                        }
                    }
                }
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
            //cout<<currentChr<<genomeLength[currentChr]<<endl;
        }
        if (score<=0) {
            continue;
        }
        if (segment.startP-1<0||segment.endP>genomeLength[currentChr]) {
            continue;
        }
        //for performance
        vector<short int> &tempMap=rawTag[thisHistone][currentChr];
        for (int j=segment.startP-1; j<segment.endP; j++) {
            
            tempMap[j]+=score;
        //            cout<<j;
        }
        
    }
    //last histone's last chrome
    getTagBed(thisHistone,currentChr);
    //append reverse for antisense,for each histone
    appendReverseTag();
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
    for (vector<string>::iterator it=chromeNames.begin(); it!=chromeNames.end(); it++) {
        genomeLength[*it] = rawGenome[*it].size();
    }
//  initProb(1);
	return true;
}



void genomeRegions::getSeq(const string& outPutDir){
    string fileName = outPutDir+"/allregions.fa";
    ofstream regionFile(fileName.c_str());
    vector<genomeRegion>::iterator it;

    //chr1-chr19-chrX
    for (it = genomes.begin(); it!=genomes.end(); it++) {
        appendSeq(regionFile,it->startP,it->endP,it->chr);
    }
    //clear memory
    sort(chromeNames.begin(), chromeNames.end(),chrCompare);
    for (vector<string>::iterator it=chromeNames.begin(); it!=chromeNames.end(); it++) {
        genomeLength[*it] = rawGenome[*it].size();
        string tempS("");
        rawGenome[*it].swap(tempS);
    }
    
    //output is also pasted sorted chr1-chrn
    
}

//appendSeq to 
void inline genomeRegions::appendSeq(ostream &outFile,int startCut,int endCut,string &currentChr){
    //cout<<currentChr<<" "<<startCut<<" "<<endCut<<" "<<rawGenome[currentChr].substr(startCut,endCut-startCut+1)<<endl;
    segmentCount++;
    if (endCut<genomeLength[currentChr]-1) {
        segmentStartPos.push_back(endCut-startCut+2+segmentStartPos.back());
    }
    else {
        segmentStartPos.push_back(genomeLength[currentChr]+1-startCut+segmentStartPos.back());
    }
    segmentGenomePos.push_back(make_pair(currentChr, startCut));
    
    try {
        if (endCut<=startCut+1) {
            return;
        }
        regionSeqs[currentChr] += rawGenome[currentChr].substr(startCut,endCut-startCut+1)+"#";
        //generate fasta
        if (option["regionfile"]=="T") {
            outFile<<">"<<currentChr<<":"<<startCut<<"-"<<endCut<<"\n";
            for (int i = startCut; i<endCut; i+=50) {
                if (i+50 > endCut) {
                    outFile<<rawGenome[currentChr].substr(i,endCut-i)<<"\n";
                    break;
                }
                outFile<<rawGenome[currentChr].substr(i,50)<<"\n";
            }
            outFile<<endl;
        }
    } catch (exception &e) {
        cerr<<e.what()<<endl;
        cerr<<"Seq insert error: chrome:"<<currentChr<<" start:"<<startCut<<" end:"<<endCut<<" chromesize"<<rawGenome[currentChr].size()<<endl;
        return;
    }

}

int genomeRegions::appendReverseGenome(char* T){
    //extend segmentStartPos
    assert(segmentStartPos.size()==segmentCount+1);
    assert(segmentGenomePos.size()==segmentCount);
    int Halfway = segmentStartPos.back();
    segmentStartPos.back()++;
    for (int i=segmentCount-1; i>=0; i--) {
        segmentStartPos.push_back(Halfway-segmentStartPos[i]+Halfway+1);
    }
    //cout<<segmentStartPos.size()<<segmentStartPos<<endl;
    //extend segmentGenomePos
    for (int i=segmentCount-1; i>=0; i--) {
        segmentGenomePos.push_back(segmentGenomePos[i]);
    } 
    vector<string>::iterator chrit;
    long int offset=0;
    string tempString("");
    for (chrit=chromeNames.begin(); chrit!=chromeNames.end(); chrit++) {
        tempString += regionSeqs[*chrit];
        string tempS("");
        regionSeqs[*chrit].swap(tempS);
    }
    tempString += antisense(tempString)+"#";
    offset = tempString.size();
    cerr<<"regionSize:"<<offset<<endl;
    assert(offset<=MAX_LENGTH);
    strcpy(T, tempString.c_str()); 
    tempString.clear();
    return offset;
}


/*
void genomeRegions::writeRawTag(genomeRegions &tagBed){
    //sort tagBed
    //sort(tagBed.genomes.begin(), tagBed.genomes.end(), compareGenome);
    string currentChr("");
    string thisHistone = "bed";
    map<string, vector<short int> > tempMap;
    rawTag[thisHistone] = tempMap;
    for (int i=0; i<tagBed.genomes.size(); i++) {
        if (tagBed.genomes[i].chr!=currentChr) {
            currentChr = tagBed.genomes[i].chr;
            vector<short int> temp;
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
*/

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
    vector<short int> tempV;
    rawTag[thisHistone][currentChr].swap(tempV);
}

void genomeRegions::appendTag(int a,int b,const string& chr,const string& thisHistone){   
    try {
        for (int i=a-1-EXTENDBOUND; i<b+EXTENDBOUND; i++) {
            //if b exceed fasta's length
            if (i>=genomeLength[chr]-1) {
                if (b>=genomeLength[chr]-1) {
                    for (int j=0; j<EXTENDBOUND; j++) {
                        regionTags[thisHistone].push_back(0);
                    }
                }
                // if b<...
                else {
                    for (int j=0; j<EXTENDBOUND-(genomeLength[chr]-1-b); j++) {
                        regionTags[thisHistone].push_back(0);
                    }
                }
                break;
            }
            else {
                if (i>=rawTag[thisHistone][chr].size()) {
                    regionTags[thisHistone].push_back(0);
                }
                else {
                    regionTags[thisHistone].push_back(rawTag[thisHistone][chr][i]);
                }
            }

        }
        regionTags[thisHistone].push_back(0);
    } catch (exception &e) {
        cerr<<e.what()<<endl;
        cerr<<"Tag insert error: chrome:"<<chr<<" start:"<<a<<" end:"<<b<<" chromesize"<<rawTag[thisHistone][chr].size()<<endl;
    }
}

void genomeRegions::appendReverseTag(){
    //extend genomeTag for each
    for (vector<string>::iterator tagIt=tagName.begin(); tagIt!=tagName.end(); tagIt++) {
        string thisHistone = *tagIt;
        vector<short int> temp(regionTags[thisHistone]);
        reverse(temp.begin(), temp.end());
        copy(temp.begin(), temp.end(),back_inserter(regionTags[thisHistone]));
        regionTags[thisHistone].push_back(0);
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
        for (it=regionSeqs.begin(); it!=regionSeqs.end(); it++ ) {
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
    if (i==total-1) {
        cerr<<endl;
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

char degenerate(char a,char b){
    //assert
    assert(a=='A'||a=='C'||a=='G'||a=='T');
    switch (a) {
        case 'A':
            switch (b) {
                case 'A':
                    return 'A';    
                case 'C':
                    return 'M';
                case 'G':
                    return 'R';
                case 'T':
                    return 'W';
                case 'M':
                    return 'M';
                case 'R':
                    return 'R';
                case 'W':
                    return 'W';
                case 'S':
                    return 'V';
                case 'Y':
                    return 'H';
                case 'K':
                    return 'D';
                case 'B':
                    return 'N';
                case 'D':
                    return 'D';
                case 'H':
                    return 'H';
                case 'V':
                    return 'V';
                case 'N':
                    return 'N';
                default:
                    break;
            }
            break;
        
        case 'C':
            switch (b) {
                case 'A':
                    return 'M';    
                case 'C':
                    return 'C';
                case 'G':
                    return 'S';
                case 'T':
                    return 'Y';
                case 'M':
                    return 'M';
                case 'R':
                    return 'V';
                case 'W':
                    return 'H';
                case 'S':
                    return 'S';
                case 'Y':
                    return 'Y';
                case 'K':
                    return 'B';
                case 'B':
                    return 'B';
                case 'D':
                    return 'N';
                case 'H':
                    return 'H';
                case 'V':
                    return 'V';
                case 'N':
                    return 'N';
                default:
                    break;
            }
            break;
            
        case 'G':
            switch (b) {
                case 'A':
                    return 'R';    
                case 'C':
                    return 'S';
                case 'G':
                    return 'G';
                case 'T':
                    return 'K';
                case 'M':
                    return 'V';
                case 'R':
                    return 'R';
                case 'W':
                    return 'D';
                case 'S':
                    return 'S';
                case 'Y':
                    return 'B';
                case 'K':
                    return 'K';
                case 'B':
                    return 'B';
                case 'D':
                    return 'D';
                case 'H':
                    return 'N';
                case 'V':
                    return 'V';
                case 'N':
                    return 'N';
                default:
                    break;
            }
            break;
            
        case 'T':
            switch (b) {
                case 'A':
                    return 'W';    
                case 'C':
                    return 'Y';
                case 'G':
                    return 'K';
                case 'T':
                    return 'T';
                case 'M':
                    return 'H';
                case 'R':
                    return 'D';
                case 'W':
                    return 'W';
                case 'S':
                    return 'B';
                case 'Y':
                    return 'Y';
                case 'K':
                    return 'K';
                case 'B':
                    return 'B';
                case 'D':
                    return 'D';
                case 'H':
                    return 'H';
                case 'V':
                    return 'N';
                case 'N':
                    return 'N';
                default:
                    break;
            }
            break;
            
        default:
            break;
    }
    return '_';
}

bool chrCompare(const string& chr1,const std::string &chr2){
    int num1,num2;
    if (chr1.find("X")!=chr1.npos) {
        num1='X';
    }
    else if (chr1.find("Y")!=chr1.npos) {
        num1='Y';
    }
    else if (chr1.find("M")!=chr1.npos) {
        num1='M';
    }
    else {
        num1=atoi(chr1.substr(3,chr1.size()-3).c_str());
    }
    
    if (chr2.find("X")!=chr2.npos) {
        num2='X';
    }
    else if (chr2.find("Y")!=chr2.npos) {
        num2='Y';
    }
    else if (chr2.find("M")!=chr2.npos) {
        num2='M';
    }
    else {
        num2=atoi(chr2.substr(3,chr2.size()-3).c_str());
    }
    return num1<num2;
}

pair<int,int> locateSubscript(const vector<int> &listObj, vector<int>::const_iterator begin,vector<int>::const_iterator end, int queryVal){
	//assume sorted
	int span = end - begin;
	int currentVal = *(begin+int(span/2));
	if (queryVal<0||queryVal>=listObj.back()){
		return std::make_pair(-1, -1);
	}
	else if (span == 1){
		return make_pair(begin-listObj.begin(), queryVal-currentVal) ;
	}
	else if (queryVal>currentVal){
		return locateSubscript(listObj, begin+span/2, end, queryVal);
	}
	else {
		return locateSubscript(listObj, begin, end-span/2, queryVal);
	}  
}