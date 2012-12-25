//
//  common.cpp
//  MOST
//
//  Created by zhang yizhe on 12-5-20.
//  Copyright (c) 2012年 SJTU. All rights reserved.
//
//  We thank John D. Cook for his help and his statistical C++ toolkit

#include "common.h"
#include "motif.h"
#include "FFT.h"
#include <bitset>

using namespace std;

int K;
float DELTA;
int MAX_SHIFT;
float MAX_DIST;
float MAX_KL_DIV;
int MAX_CLUSTER_NUM;
int MAX_WORD_NUM;
long int HASH_TABLE_SIZE;
int PEUSUDO_COUNT;

map<string,string> option;


/***********************************************************************
 * parseCommandLine
 *
 *
 ***********************************************************************/
void parseCommandLine(int argc,
                      char** argv
                      ) {
    
    // Set default values for command line arguments
    
    map<string, bool> optionRequire;
    //t:tag n:normal
    option          ["mode"]        =   "normal";
    optionRequire   ["mode"]        =   true;
    //
    option          ["fastafile"]   =   "";
    optionRequire   ["fastafile"]   =   true;
    //
    option          ["regionfastafile"]   =   "";
    optionRequire   ["regionfastafile"]   =   false;
    //
    option          ["regionfile"]  =   "";
    optionRequire   ["regionfile"]  =   true;
    //
    option          ["tagfile"]     =   "";
    optionRequire   ["tagfile"]     =   false;
    //
    option          ["outdir"]      =   "MOST_out";
    optionRequire   ["outdir"]      =   false;
    //
    option          ["k"]      =   "9";
    optionRequire   ["k"]      =   false;
    //
    option          ["order"]      =   "0";
    optionRequire   ["order"]      =   false;
    //
    option          ["writeloci"]      =  "F";
    optionRequire   ["writeloci"]      =   false;
    //
    option          ["writeregion"]      =  "F";
    optionRequire   ["writeregion"]      =   false;
    //if remove control
    option          ["drawdist"]      =  "T";
    optionRequire   ["drawdist"]      =   false;
    //
    option          ["bkgregion"]      =  "region";
    optionRequire   ["bkgregion"]      =   false;
    //
    option          ["FFT"]      =  "T";
    optionRequire   ["FFT"]      =   false;
    //
    option          ["trim"]      =  "50";
    optionRequire   ["trim"]      =   false;
    //
    option          ["clusterlength"]      =  "20";
    optionRequire   ["clusterlength"]      =   false;
    //p-value for over-representation
    option          ["pvalue"]      =  "0.05";
    optionRequire   ["pvalue"]      =   false;
    //if remove repeat
    option          ["rmrepeat"]      =  "2";
    optionRequire   ["rmrepeat"]      =   false;
    //
    option          ["clusterStringency"]      =  "0.2";
    optionRequire   ["clusterStringency"]      =   false;
    //
    option          ["MAX_CLUSTER_NUM"]      =  "1000";
    optionRequire   ["MAX_CLUSTER_NUM"]      =   false;
    //if remove repeat
    option          ["extendmotif"]      =  "T";
    optionRequire   ["extendmotif"]      =   false;
    //
    option          ["MAX_WORD_NUM"]      =  "1000";
    optionRequire   ["MAX_WORD_NUM"]      =   false;
    //
    option          ["flanking"]      =  "0";
    optionRequire   ["flanking"]      =   false;
    //if remove repeat
    option          ["ROC"]      =  "F";
    optionRequire   ["ROC"]      =   false;
    //if remove control
    option          ["control"]      =  "";
    optionRequire   ["control"]      =   false;
    //if remove control
    option          ["peusudo"]      =  "1";
    optionRequire   ["peusudo"]      =   false;

    // Parse the command line.
    string option_name = "";
    string option_value = "";
    
    
    // Read the next option, and break if we're done.
    if (argc==1) {
        printUsage();
        exit(0);
    }
    for (int i=1; i<argc+1; i++) {
        if (i==argc) {
            goto CHECKPARA;
        }
        else if (argv[i][0]!='-') {
            option_value+=argv[i];
            option_value+=" ";
        }
        else {
        CHECKPARA:
            // Assign the parsed value to the appropriate variable
            // and validate each parameter input
            trim(option_value);
            if (option_name=="") {
                //do nothing
            }
            else if (option_name == "-h"||option_name == "--help") {
                printUsage();
                exit(0);
            }
            else if (option_name == "--version") {
                cerr<<MOSTVERSION;
                cerr<<"by Yizhe Zhang";
                exit(0);
            }
            else if (option_name == "-m") {
                option["mode"] = option_value;
                optionRequire["mode"] = false;
                if (option_value=="tag") {
                    optionRequire["tagfile"] = true;
                }
                else if (option_value=="normal") {
                    
                }
                else if (option_value=="help") {
                    printUsage();
                    exit(0);
                }
                else if (option_value=="extract") {
                    
                }
                else {
                    cerr<<"FATAL:"<<"wrong mode!"<<endl;
                    exit(1);
                }
            } 
            else if (option_name == "-f") {
                option["fastafile"] = option_value;
                optionRequire["fastafile"] = false;
            } 
            else if (option_name == "-rf") {
                option["regionfastafile"] = option_value;
                optionRequire["regionfastafile"] = false;
            } 
            else if (option_name == "-r") {
                option["regionfile"] = option_value;
                optionRequire["regionfile"] = false;
            } 
            else if (option_name == "-t") {
                option["tagfile"] = option_value;
                optionRequire["tagfile"] = false;
            } 
            else if (option_name == "-c") {
                option["control"] = option_value;
                optionRequire["control"] = false;
            } 
            else if (option_name == "-o") {
                option["outdir"] = option_value;
            } 
            else if (option_name == "-k") {
                option["k"] = option_value;
                if (atoi(option["k"].c_str())<4||
                    atoi(option["k"].c_str())>15) {
                    printAndExit("K must be within [4,15]");
                }
                
            } 
            else if (option_name == "-bo") {
                option["order"] = option_value;
                if (atoi(option["order"].c_str())<-1||
                    atoi(option["order"].c_str())>1) {
                    printAndExit("order must be within [-1,1]");
                }
            }
            else if (option_name == "-locifile") {
                option["writeloci"] = option_value;
                if (option["writeloci"]!="T"&&option["writeloci"]!="F") {
                    printAndExit("writeloci value should be either T or F");
                }
            }
            else if (option_name == "-writeregion") {
                option["writeregion"] = option_value;
                if (option["writeregion"]!="T"&&option["writeregion"]!="F") {
                    printAndExit("writeregion value should be either T or F");
                }
            }
            else if (option_name == "-drawdist"||option_name == "-dd") {
                option["drawdist"] = option_value;
                if (option["drawdist"]!="T"&&option["drawdist"]!="F") {
                    printAndExit("drawdist value should be either T or F");
                }
            }
            else if (option_name == "-br") {
                option["bkgregion"] = option_value;
            }
            else if (option_name == "-fft") {
                option["FFT"] = option_value;
                if (option["FFT"]=="T"&&option_value!="tag") {
                    printAndExit("Mode must be tag if FFT is chosen");
                }
                if (option["FFT"]!="T"&&option["FFT"]!="F") {
                    printAndExit("FFT value should be either T or F");
                }
            }
            else if (option_name == "-trim") {
                option["trim"] = option_value;
                if (atoi(option["trim"].c_str())<10||
                    atoi(option["trim"].c_str())>100) {
                    printAndExit("trim must be within [10,100]");
                }
            }
            else if (option_name == "-cl") {
                option["clusterlength"] = option_value;
                if (atoi(option["clusterlength"].c_str())<4||
                    atoi(option["clusterlength"].c_str())>30) {
                    printAndExit("clusterlength must be within [4,30]");
                }
            }
            else if (option_name == "-pvalue"||option_name == "-p") {
                option["pvalue"] = option_value;
                if (atof(option["pvalue"].c_str())>=1||
                    atof(option["pvalue"].c_str())<0) {
                    printAndExit("pvalue should be within (0,1), <0.05 recommend");
                }
            }
            else if (option_name == "-rmrepeat") {
                option["rmrepeat"] = option_value;
                if (atoi(option["rmrepeat"].c_str())<0||
                    atoi(option["rmrepeat"].c_str())>2) {
                    printAndExit("rmrepeat value should be either 0,1 or 2");;
                }
            }
            else if (option_name == "-cs") {
                option["clusterStringency"] = option_value;
                if (atof(option["clusterStringency"].c_str())>=1||
                    atof(option["clusterStringency"].c_str())<0) {
                    printAndExit("clusterStringency should be within [0,1), <0.05 recommend");
                }
            }
            else if (option_name == "-cn") {
                option["MAX_CLUSTER_NUM"] = option_value;
                if (atoi(option["MAX_CLUSTER_NUM"].c_str())<2||
                    atoi(option["MAX_CLUSTER_NUM"].c_str())>500) {
                    printAndExit("cluster num must be within [1,500]");
                }
            }
            else if (option_name == "-mn") {
                option["MAX_WORD_NUM"] = option_value;
                if (atoi(option["MAX_WORD_NUM"].c_str())<2||
                    atoi(option["MAX_WORD_NUM"].c_str())>20000) {
                    printAndExit("cluster num must be within [1,20000]");
                }
            }
            else if (option_name == "-extend") {
                option["extendmotif"] = option_value;
                if (option["extendmotif"]!="T"&&option["extendmotif"]!="F") {
                    printAndExit("Whether Extend motif should be either T or F");
                }
            }
            else if (option_name == "-fl") {
                option["flanking"] = option_value;
                if (atoi(option["flanking"].c_str())<2||
                    atoi(option["flanking"].c_str())>20000) {
                    printAndExit("extender must be within [1,20000]");
                }
            }
            else if (option_name == "-roc") {
                option["ROC"] = option_value;
                if (option["ROC"]!="T"&&option["ROC"]!="F") {
                    printAndExit("ROC value should be either T or F");
                }           
            }
            else if (option_name == "-ps") {
                option["peusudo"] = option_value;
                if (atoi(option["peusudo"].c_str())<0||
                    atoi(option["peusudo"].c_str())>200) {
                    printAndExit("Peusudo count must be within [0,200]");
                }          
            }
            else {
                cerr<<"unrecognized option "<<option_name<<"! skip"<<endl;
            } 
            
            //update
            if (i==argc) {
                break;
            }
            option_value = "";
            option_name = argv[i];
            transform(option_name.begin(), option_name.end(), option_name.begin(), ::tolower);
        }
        
        
        
    }
    //output commandlines:
    cerr<<"\nCommandline:";
    for (map<string, string>::iterator it=option.begin();it!=option.end();it++) {
        cerr<<"(--"<<it->first<<" "<<it->second<<") ";
    }
    cerr<<"\n";
    // verify required arguments .
    if (option["regionfastafile"]!=""&&option["mode"]=="tag") {
        option["mode"]="normal";
        cerr<<"WARNING:Automatically choosing normal mode when -rf is specified"<<"\n";
    }
    
    for (map<string, bool>::iterator it=optionRequire.begin();it!=optionRequire.end();it++) {
        if (it->second) {
            cerr<<"FATAL:"<<"option "<<it->first<<" need to be specified in command line!"<<endl;
            exit(1);
        }
    }
    //changing parameters
    if (option["mode"]=="extract")
        option["writeregion"]="T";
    K = atoi(option["k"].c_str());
    float stringency = atof(option["clusterStringency"].c_str());
    MAX_SHIFT = int(K*(1-stringency)/2.5+1);
    MAX_DIST = K/3.0*(1-stringency)+1;
    MAX_KL_DIV = float(0.001+0.004*(1-stringency));
    DELTA = NormalCDFInverse(1.0-atof(option["pvalue"].c_str()))+1.0;
    MAX_WORD_NUM = atoi(option["MAX_WORD_NUM"].c_str());
    MAX_CLUSTER_NUM = atoi(option["MAX_CLUSTER_NUM"].c_str());
    PEUSUDO_COUNT = atoi(option["peusudo"].c_str());

    
}

void printUsage()
{
	string usage =	"***********************************************************************************************\n";
	usage		+=	"                             MOST " + MOSTVERSION+"\n";
	usage		+=	"                             Developed by Yizhe Zhang, Chaochun Wei\n";
	usage		+=	"                             jeremy071242044@gmail.com\n";
	usage		+=	"***********************************************************************************************\n\n";
	usage		+=	"    ./most [option: <parameter1> <parameter2>]  \n\n";
    usage		+=	"    -h help                    Print help information.\n\n";
	usage		+=	"    Required Parameters:\n\n";
	usage		+=	"    -m <normal/tag/help/extract>       Runing mode:tagfile should be specified if tag mode is selected\n\n";
	usage		+=	"    -f <DNA sequence file>     Whole genome DNA sequences\n\n";
    //usage		+=	"    -rf <DNA sequence file>    Regional DNA 
    usage		+=	"    -r <BED file>              Regions of interest\n\n";
    usage		+=	"    -t <WIG file>              Tag file for Histone marks or other sources\n\n";
    //usage		+=	"    CAVEAT:WIG FILE SHOULD BE SORTED\n\n";
    usage		+=	"    Optional Parameters:\n\n";
    usage		+=	"    -c <BEG file>              Control file for ChIP-seq experiment\n\n";
    usage       +=  "    -extend <T/F>              Whether or not to extend motif\n\n";
    usage		+=	"    -o <outputDIR>             Specify an output directory(default MOST_out).\n\n";
    usage		+=	"    -bo <0,1>                  Order for background sequence(default 0)\n\n";
    usage		+=	"    -br <region/genome>        Specify background sequence(default region)\n\n";
    usage		+=	"    -locifile <T/F>            Whether or not loci file of each cluster should be exported(default T)\n\n";
    //usage		+=	"    -writeregion <T/F>          Whether or not region file of each cluster should be exported(default F)\n\n";
    usage		+=	"    -drawdist <T/F>            Whether or not tag distribution file of each cluster should be exported(default T)\n\n";
    usage		+=	"    -fft <T/F>                 Whether or not to cast FFT on each bins(default T)\n\n";
    usage		+=	"    -trim <10-100 integer>     Larger trimer means larger length of cluster(default 50)\n\n";
    usage		+=	"    -cl <4-30 integer>         maximal length of single cluster(default 20)\n\n";
    usage		+=	"    -p(-pvalue) <0-1 float>    P-value for low occurence word filtering(default 0.05)\n\n";
    usage		+=	"    -rmrepeat <0,1,2>          Whether or not to remove repeat words(default 1)\n\n";
    usage		+=	"    -cs <0-1 float>            Stringency of clustering(default 0.2)\n\n";
    usage		+=	"    -cn <1-500 integer>        maximal number of clusters(default 1000)\n\n";
    usage		+=	"    -mn <1-20000 integer>      maximal number of qualified words to be clustered(default 1000)\n\n";
    //usage		+=	"    Testing Parameters:\n\n";
    //usage		+=	"    -roc <T/F>                 Plot roc(default F)\n\n";
    usage		+=	"    -fl <1-20000 integer>  Extract extender(default 0)\n\n";
    usage		+=	"    -ps <0-200 integer>        Peusudo count added on PFM(default 1)\n\n";
    usage               +=      "    Examples:\n\n";
    usage		+=	"    ./most -m normal -r VDR-sti.chr3.bed -f chr3.fa.masked -o VDRnormal >VDR_normal.txt \n\n";
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
    string thistag("");
    //mode 1 (chr1 800 1.25)density    mode 2 (800 1.25)wig   mode 3 (1.25 fixed)fixwig
    int mode = 0;  
    int fixedPos = 1;
    string WigChr;
    int lineCounter=0;
    map<int, string> modelabel;
    modelabel[1]="density";
    modelabel[2]="variableStep wig";
    modelabel[3]="fixedStep wig";
    while (getline(wigFile, line)) {
        lineCounter++;
        //new histone/tag
        if (line[0]=='v'||line[0]=='f') {
            size_t found;
            //if mode is 1 check else and put last histone's last chr 
            if (currentChr!="") {
                getTagBed(thistag,currentChr);
            }
            
            //find info 
            if (lineCounter==1) {
                found = line.find("histone=");
                if (found!=string::npos){
                    size_t found2 = line.substr(found+8,line.npos).find_first_of(" "); 
                    tagName.push_back(line.substr(found+8,found2));
                    
                }
                else {
                    char temptagName[20];
                    sprintf (temptagName, "TagNo.%lu", tagName.size()+1);
                    tagName.push_back(temptagName);
                    cerr<<" You'd better contain histone=**"<<endl;
                }
            }
            thistag = tagName.back();
            
            //find chr info,determine mode
            found = line.find("chrom="); 
            if (found==string::npos||line.substr(found+6,6)=="chrAll") {
                mode = 1;
            }
            else {
                string tempS=line.substr(found+6,5);
                WigChr = trim(tempS);
                rawTag[thistag][WigChr].assign(genomeLength[WigChr], 0); 
                if (line[0]=='f') {
                    mode = 3;
                    size_t found1 = line.find("start="); 
                    if (found1==string::npos) {
                        printAndExit("check wig file,fixedStep should have \"start=\"");
                    }
                    else {
                        string tempS=line.substr(found1+6,2);
                        fixedPos = atoi(trim(tempS).c_str());
                    }
                }
                else {
                    mode = 2;
                }
            }
            if (lineCounter==1) {
                cerr<<"Find Tag "<<tagName.back()<<" mode "<<modelabel[mode]<<endl;
            }
            //span= or step=
            found = line.find("span=");
            size_t found1 = line.find("step=");
            if (found!=string::npos&&found1!=string::npos) {
                printAndExit("check wig file's info line: shouldn't contain both span= and step= ");
            }
            else if (found==string::npos&&found1==string::npos) {
                printAndExit("check wig file's info line: should contain span=**");
            }
            else {
                found = min(found, found1);
                wigStep=atoi(line.substr(found+5,line.size()-found-5).c_str());
            }
            if (mode==1) {
                currentChr = "";
            }
            continue;
        }
		if(line.size() == 0) continue;
		/* Read wig line */
		istringstream ss(line);
        int score;
		if (mode == 1) {
            string temp;
            ss>>segment.chr;
            if (!existChr(segment.chr)) continue;
            ss>>temp;
            segment.startP = atoi(temp.c_str());
            if (segment.startP<0) continue;
            segment.endP = segment.startP + wigStep;
            ss>>temp;
            score = int(atof(temp.c_str())*4);
        }
        else if (mode == 2) {
            string temp;
            currentChr=WigChr;
            ss>>temp;
            segment.startP = atoi(temp.c_str());
            if (segment.startP<0) continue;
            segment.endP = segment.startP + wigStep;
            ss>>temp;
            score = atoi(temp.c_str());
        }
        //mode = 3
        else {
            string temp;
            currentChr=WigChr;
            fixedPos += wigStep;
            segment.startP = fixedPos;
            segment.endP = segment.startP + wigStep;
            ss>>temp;
            score = atoi(temp.c_str());
        }
		
        //find new chr  push old and detect lack
        if (mode==1) {
            if (segment.chr!=currentChr) {
                //push old
                if (currentChr!="") {
                    //put last one
                    getTagBed(thistag,currentChr);
                    
                }
                
                
                
                currentChr = segment.chr;
                //or find  != .end()
                if (rawTag[thistag].count(currentChr)) {
                    string errorInfo = "Error! chrome "+currentChr+" information collides in Tag file\n";
                    cout<<errorInfo;
                    continue;
                }
                
                rawTag[thistag][currentChr].assign(genomeLength[currentChr], 0);  
            }
        }
        
        if (score<=0) {
            continue;
        }
        if (segment.startP-1<0||segment.endP>genomeLength[currentChr]) {
            continue;
        }
        //for performance
        vector<short int> &tempMap=rawTag[thistag][currentChr];
        for (int j=segment.startP-1; j<segment.endP; j++) {
            tempMap[j]+=score;
        //            cout<<j;
        }
        
    }
    //last tag's last chrome
    getTagBed(thistag,currentChr);
    
    //check histone contains all chrome info that is needed
    for (vector<string>::iterator ChrIt=chromeNames.begin(); ChrIt!=chromeNames.end(); ChrIt++) {
        if (rawTagSegments[thistag].find(*ChrIt)==rawTagSegments[thistag].end()) {
            cerr<<"WARNING:"<<thistag<<" in TagFile doesn't contain chrome "<<(*ChrIt)<<" appears in fastaFile"<<endl;
            rawTag[thistag][*ChrIt].assign(genomeLength[*ChrIt], 0);
            getTagBed(thistag,*ChrIt);
        }
    }
    return true;
    
}
bool genomeRegions::catenateTags(){
    
    //catenate tag
    assert(regionTags.size()==0);
    for (int i=0; i<tagName.size(); i++) {
        for (vector<string>::iterator ChrIt=chromeNames.begin(); ChrIt!=chromeNames.end(); ChrIt++) {
            regionTags[tagName[i]].insert(regionTags[tagName[i]].end(),rawTagSegments[tagName[i]][*ChrIt].begin(), rawTagSegments[tagName[i]][*ChrIt].end());
            vector<short int> temp;
            rawTagSegments[tagName[i]][*ChrIt].swap(temp);
        }
    }   
    //append reverse for antisense,for each tag
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

bool genomeRegions::readRegionFasta(const string &filename){
    ifstream fastaFile(filename.c_str());
    if (!fastaFile) {
        string errorInfo = "Error! Fail to open fasta file \"" + filename + "\" !";
		printAndExit(errorInfo);
    }
    
    string chr("region");
    string line;
    while(getline(fastaFile,line))
	{
		trim(line);
		if(line.size() == 0)    continue;
		if(line[0] == '>')      continue;
        //if seq line
        else {
            transform(line.begin(), line.end(), line.begin(), ::toupper);
            regionSeqs[chr] += line +"#";
        }
    }
    genomeLength[chr] = regionSeqs[chr].size();
    chromeNames.push_back(chr);
    //  initProb(1);
	return true;
}


//rm control Peaks
void genomeRegions::rmControlPeaks(const string& filename){
    genomeRegions tempgR(0);
    tempgR.chromeNames=chromeNames;
    tempgR.readBed(filename);
    
    //sort(tempgR.genomes.begin(), tempgR.genomes.end(), compareGenome);
    vector<Interval> intervalLib;
    vector<genomeRegion>::iterator genomeit;
    //origin
    for (genomeit = genomes.begin(); genomeit!=genomes.end(); genomeit++) {
        Interval tempS={genomeit->startP,true,genomeit->chr,'o'};
        intervalLib.push_back(tempS);
        Interval tempE={genomeit->endP,false,genomeit->chr,'o'};
        intervalLib.push_back(tempE);
    }
    //control
    for (genomeit = tempgR.genomes.begin(); genomeit!=tempgR.genomes.end(); genomeit++) {
        Interval tempS={genomeit->startP,true,genomeit->chr,'c'};
        intervalLib.push_back(tempS);
        Interval tempE={genomeit->endP,false,genomeit->chr,'c'};
        intervalLib.push_back(tempE);
    }
    callOverLaps(intervalLib);
    return;

}
//mask region labeled "control"
void genomeRegions::callOverLaps(vector<Interval> &intervalLib) {
    sort(intervalLib.begin(), intervalLib.end(),compareInterval);
    //can deal with overlapped control set and orginal set
    assert(intervalLib.size()>1);
    bool _originOpen=false;
    bool _controlOpen=false;
    int _currentS = 0;
    //initiate chr
    string _currentChr = intervalLib[0].chr;
    vector<genomeRegion> masked_region;
    for (int i=0; i<intervalLib.size(); i++) {
        if (intervalLib[i].chr==_currentChr) {
            if ((intervalLib[i].start&&intervalLib[i].tag=='c')||
                (!intervalLib[i].start&&intervalLib[i].tag=='o')) {
                if (_originOpen&&(!_controlOpen)) {
                    //push back this segment
                    if (intervalLib[i].Point>_currentS+K) {
                        genomeRegion tempG={_currentS,intervalLib[i].Point,_currentChr};
                        masked_region.push_back(tempG);
                    }
                }
            }
            else {
                //update to new start
                _currentS = intervalLib[i].Point;
            }
            //switch open accessibility
            if (intervalLib[i].tag=='o') {
                _originOpen=intervalLib[i].start;
            }
            if (intervalLib[i].tag=='c') {
                _controlOpen=intervalLib[i].start;
            }
        }
        else {
            _currentChr=intervalLib[i].chr;
            _originOpen=false;
            _controlOpen=false;
        }
    }
    genomes.swap(masked_region);
    sort(genomes.begin(), genomes.end(),compareGenome);
}


//for allchrom
void genomeRegions::getSeq(const string& outPutDir){
    string fileName = outPutDir+"/allregions.fa";
    ofstream regionFile(fileName.c_str());
    vector<genomeRegion>::iterator it;

    //chr1-chr19-chrX   iterator changed in it
    for (it = genomes.begin(); it!=genomes.end();) {
        //cerr<<it->startP<<it->chr<<it->endP<<" "<<genomes.size()<<endl;
        appendSeq(regionFile,it);
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
void genomeRegions::appendSeq(ostream &outFile,vector<genomeRegion>::iterator& currentGenomeRegion){
    int startCut=currentGenomeRegion->startP;
    int endCut=currentGenomeRegion->endP;
    string &currentChr=currentGenomeRegion->chr;
    //cout<<currentChr<<" "<<startCut<<" "<<endCut<<endl;
    
    
    try {
        assert(endCut>startCut+1);

        regionSeqs[currentChr] += rawGenome[currentChr].substr(startCut,endCut-startCut+1)+"#";
        //generate fasta
        if (option["writeregion"]=="T") {
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
        //update segmentstart pos
        segmentCount++;
        if (endCut<genomeLength[currentChr]-1) {
            segmentStartPos.push_back(endCut-startCut+2+segmentStartPos.back());
        }
        else {
            segmentStartPos.push_back(genomeLength[currentChr]+1-startCut+segmentStartPos.back());
        }
        segmentGenomePos.push_back(make_pair(currentChr, startCut));
        //next iter
        currentGenomeRegion++;
    } catch (exception &e) {
        currentGenomeRegion = genomes.erase(currentGenomeRegion);
        
        cerr<<e.what()<<endl;
        cerr<<"Seq insert error: chrome:"<<currentChr<<" start:"<<startCut<<" end:"<<endCut<<" chromesize"<<rawGenome[currentChr].size()<<endl;

    }
    return;
}

int genomeRegions::appendReverseGenome(string& T){
    //extend segmentStartPos
    //cerr<<segmentStartPos.size()<<segmentCount+1<<endl;
    assert(segmentStartPos.size()==segmentCount+1);
    assert(segmentGenomePos.size()==segmentCount);
    int Halfway = segmentStartPos.back();
    segmentStartPos.back()++;
    for (int i=segmentCount-1; i>=0; i--) {
        segmentStartPos.push_back(Halfway-segmentStartPos[i]+Halfway+1);
    }
    //extend segmentGenomePos
    for (int i=segmentCount-1; i>=0; i--) {
        segmentGenomePos.push_back(segmentGenomePos[i]);
    } 
    vector<string>::iterator chrit;
    T.clear();
    T.reserve(MAX_LENGTH);
    for (chrit=chromeNames.begin(); chrit!=chromeNames.end(); chrit++) {
        T += regionSeqs[*chrit];
        string tempS("");
        regionSeqs[*chrit].swap(tempS);
    }
    T += antisense(T)+"#";
    cerr<<"regionSize:"<<T.size()<<endl;
    //assert(OFF_SET<=MAX_LENGTH);
    return T.size();
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
    } catch (exception &e) {
        cerr<<e.what()<<endl;
        cerr<<"chrome:"<<currentChr<<" start:"<<startCut<<" end:"<<endCut<<endl;
    }
    genomes.swap(mergedGenome);
}

//for one chromesome
void genomeRegions::getTagBed(const string& thistag,const string& currentChr){   
    vector<genomeRegion>::iterator it;
    bool rear = false;
    for (it = genomes.begin(); it!=genomes.end(); it++) {
        if (it->chr!=currentChr) {
            if (rear) 
                break;
            else 
                continue;
        }
        else {
            appendTag(it->startP,it->endP,currentChr,thistag);
            rear=true;

        }
    }
    //clear memory
    vector<short int> tempV;
    rawTag[thistag][currentChr].swap(tempV);
}

void genomeRegions::appendTag(int a,int b,const string& chr,const string& thistag){   
    try {
        a -= EXTEND_BOUND+1;
        b += EXTEND_BOUND;
        int end=b;
        //assign this segment with 0
        if (b>genomeLength[chr]-1) {
            end=genomeLength[chr]-1;
        }
        rawTagSegments[thistag][chr].insert(rawTagSegments[thistag][chr].end(),rawTag[thistag][chr].begin()+a, rawTag[thistag][chr].begin()+end);
        if (end!=b) {
            vector<int> temp;
            if (b-EXTEND_BOUND>genomeLength[chr]-1) {
                temp.assign(EXTEND_BOUND, 0);
            }
            else {
                temp.assign(b-end,0);
            }
            rawTagSegments[thistag][chr].insert(rawTagSegments[thistag][chr].end(),temp.begin(),temp.end());
        }
        rawTagSegments[thistag][chr].push_back(0);
    } catch (exception &e) {
        cerr<<e.what()<<endl;
        cerr<<"Tag insert error: chrome:"<<chr<<" start:"<<a<<" end:"<<b<<" chromesize"<<rawTag[thistag][chr].size()<<endl;
    }
}

void genomeRegions::appendReverseTag(){
    //extend genomeTag for each
    for (vector<string>::iterator tagIt=tagName.begin(); tagIt!=tagName.end(); tagIt++) {
        string thistag = *tagIt;
        vector<short int> temp(regionTags[thistag]);
        reverse(temp.begin(), temp.end());
        copy(temp.begin(), temp.end(),back_inserter(regionTags[thistag]));
        regionTags[thistag].push_back(0);
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
    static int innerCounter=0;
    if (i!=0&&innerCounter==0) {
        return;
    }
    if (i==0) {
        cerr<<message<<endl;
        cerr<<"|0%_______________50%______________100%|\n";
        
    }
    if (total<40) {
        cerr<<"=";
        innerCounter++;
    }
    else if (i%(int((total)/40))==0) {
        cerr<<"=";
        innerCounter++;
    }
    if (i==total-1&&innerCounter<40) { 
        for (int j=innerCounter; j<40; j++) {
            cerr<<"=";
        }
        innerCounter=40;
    }
    if (innerCounter>=40) {
        cerr<<endl;
        innerCounter=0;
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
            return b;
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



float symKLDiv(const vector<float> &lhs, const vector<float> &rhs){
    assert(lhs.size()==rhs.size());
    float temp = 0;
    for (int t=0; t<lhs.size(); t++) {
        temp += lhs[t]*log10f((lhs[t]+SMALLNUM)/(rhs[t]+SMALLNUM))+rhs[t]*log10f((rhs[t]+SMALLNUM)/(lhs[t]+SMALLNUM));
    }
    return temp;
}

float normalization(vector<float>& fvec){
    float total;
    for (int i=0; i<fvec.size(); i++) {
        total += fvec[i];
    }
    for (int i=0; i<fvec.size(); i++) {
        fvec[i]/=total;
    }
    return total;
}

//statistic snippet implement from http://www.johndcook.com/cpp_phi.html by John D. Cook, revised a little bit 
float calPhi(float x)
{
    // constants
    const float a1 =  0.254829592;
    const float a2 = -0.284496736;
    const float a3 =  1.421413741;
    const float a4 = -1.453152027;
    const float a5 =  1.061405429;
    const float p  =  0.3275911;
    
    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x)/sqrt(2.0);
    
    // A&S formula 7.1.26
    float t = 1.0/(1.0 + p*x);
    float y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);
    
    return 0.5*(1.0 + sign*y);
}

// compute log(1+x) without losing precision for small values of x
float LogOnePlusX(float x)
{   try{
        if (x <= -1.0)
        {
            std::stringstream os;
            os << "Invalid input argument (" << x 
            << "); must be greater than -1.0";
            throw std::invalid_argument( os.str() );
            
        }
        
        if (fabs(x) > 1e-4)
        {
            // x is large enough that the obvious evaluation is OK
            return log(1.0 + x);
        }
        
        // Use Taylor approx. log(1 + x) = x - x^2/2 with error roughly x^3/3
        // Since |x| < 10^-4, |x|^3 < 10^-12, relative error less than 10^-8
        
        return (-0.5*x + 1.0)*x;
    }
    catch (exception &e){
        cerr<<e.what()<<endl;
        return 0;
    }
}

float RationalApproximation(float t)
{
    // Abramowitz and Stegun formula 26.2.23.
    // The absolute value of the error should be less than 4.5 e-4.
    float c[] = {2.515517, 0.802853, 0.010328};
    float d[] = {1.432788, 0.189269, 0.001308};
    return t - ((c[2]*t + c[1])*t + c[0]) / 
    (((d[2]*t + d[1])*t + d[0])*t + 1.0);
}

float NormalCDFInverse(float p)
{
    if (p <= 0.0 || p >= 1.0)
    {
        std::stringstream os;
        os << "Invalid input argument (" << p 
        << "); must be larger than 0 but less than 1.";
        //throw std::invalid_argument( os.str() );
        printAndExit(os.str());
    }
    
    // See article above for explanation of this section.
    if (p < 0.5)
    {
        // F^-1(p) = - G^-1(p)
        return -RationalApproximation( sqrt(-2.0*log(p)) );
    }
    else
    {
        // F^-1(p) = G^-1(1-p)
        return RationalApproximation( sqrt(-2.0*log(1-p)) );
    }
    
}


