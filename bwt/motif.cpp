//
//  motif.cpp
//  bwt
//
//  Created by zhang yizhe on 12-5-20.
//  Copyright (c) 2012年 SJTU. All rights reserved.
//


#include "motif.h"
#include "bwt.h"

using namespace std;

int RegionSize;
/*
int Motif::distance[4][15] = {  {0,1,1,1,0,0,0,1,1,1,1,0,0,0,0},
                                {1,0,1,1,0,1,1,0,0,1,0,1,0,0,0}, 
                                {1,1,0,1,1,0,1,0,1,0,0,0,1,0,0},
                                {1,1,1,0,1,1,0,1,0,0,0,0,0,1,0}};
int Motif::distance[4][15] = {  {0,6,6,6,1,1,1,5,5,5,4,2,2,2,3},
                                {6,0,6,6,1,5,5,1,1,5,2,4,2,2,3}, 
                                {6,6,0,6,5,1,5,1,5,1,2,2,4,2,3},
                                {6,6,6,0,5,5,1,5,1,1,2,2,2,4,3}};
*/
/*
vector<int> Motif::explainMotif(){
    //find all possible explanation(generalization)
    vector<int> tempV;
    expMotifs.swap(tempV);
    if (noWildcard()) {
        expMotifs.push_back(index);
        return expMotifs;
    }
    string tempString("");
    int tempIndex=index;
    int currentNum;
    int traverseIndex = 1;
    //for 0: 1-4 traverse
    vector<int> traversalPoint;
    vector<char> x(query.size(),'_');
    //0:* 1:A 2:C 3:G 4:T 
    for (int i=0; i<query.size(); i++){
        if (query[i]!='A'&&query[i]!='T'&&query[i]!='C'&&query[i]!='G'){
            traversalPoint.push_back(i);
            x[i] = 'A';
        }
    }
    sort(traversalPoint.begin(), traversalPoint.end());
    while (true) {
        for (int i=0; i<query.size(); i++) {
            if (query[i]!='A'&&query[i]!='T'&&query[i]!='C'&&query[i]!='G') {
                tempString.push_back(query[i]);
            }
            else {
                tempString.push_back(x[i]);
            }
            
        }
        Clusterexp.push_back(tempString);
        tempString = "";
        //if over
        
        if(ascending(x,traversalPoint,traverseIndex))
            break;
    }
    return expMotifs;
}
*/

void Motif::locateMotif(const char T[]){
    //implement1 low efficiency
    bool flag = true;
    loci.clear();
    //cout<<RegionSize<<endl;
    for (int i=0; i<RegionSize; i++) {
        for (int index=0; index<query.size(); index++) {
            if((query[index]!='N'&&T[i+index]!=query[index])||T[i+index]=='#'||T[i+index]=='N'){
                flag = false;
                break;
            }
        }
        if (flag) {
            loci.push_back(i);
        }
        else {
            flag = true;
        }
    }     
}


inline int Motif::mapLoci(int genomePos,const vector<int>& StartPs){
    //subscript and dist
    pair<int, int> no_dist = locateSubscript(StartPs, StartPs.begin(), StartPs.end(), genomePos);
    if (no_dist.first==-1) {
        cerr<<genomePos<<endl;
        return 0;
    }
    return EXTENDBOUND*(no_dist.first*2+1)+genomePos;   
}



void Motif::testMotifTag(genomeRegions &gR,bool ifDraw){
    if (option["mode"]!="tag")  return;
    //protocol: input:   loci
    //          outlet: intensity noise bipeak asymmetry sumbin
    if(loci.size()==0){
        initBin(gR);
        return;
    }
    //only need loci and tagSeq
    tagBiPeak.clear();
    tagSymmetry.clear();
    tagNoise.clear();
    //if not cluster
    if (index!=-1) {
        initBin(gR);
        for (int t=0; t<gR.tagName.size(); t++) {
            vector<short int> &tag = gR.regionTags[gR.tagName[t]];
            string tagName = gR.tagName[t];
            //what's wrong with query.size()?  answer:   unsigned int always render the expression to be possitive
            //int counter = loci.size();
            float totalAround=0;
            for (int i=0; i<loci.size(); i++) {
                if (loci[i]<SAMPLESIZE*BINSPAN+SAMPLESIZE*offset+1||loci[i]>RegionSize-(SAMPLESIZE+1)*BINSPAN-SAMPLESIZE*offset+1) {
                    //counter--;
                    continue;
                }
                int thisLociCenter = mapLoci(loci[i],gR.segmentStartPos);
                //cerr<<tag[thisLociCenter]<<" "<<thisLociCenter<<endl;
                int tempsize = gR.segmentStartPos.size()/2;
                //+ direction
                if (loci[i]<gR.segmentStartPos[tempsize]) {
                    for (int l=0; l<2*SAMPLESIZE+1; l++) {
                        for (int j=(l-SAMPLESIZE)*BINSPAN+offset*(l-SAMPLESIZE); j<(l+1-SAMPLESIZE)*BINSPAN+offset*(l-SAMPLESIZE); j++) {
                            sumBin[l][t] += tag[thisLociCenter+j];
                        }
                    }
                }
                //- direction
                else {
                    for (int l=0; l<2*SAMPLESIZE+1; l++) {
                        for (int j=(l-SAMPLESIZE)*BINSPAN+offset*(l-SAMPLESIZE); j<(l+1-SAMPLESIZE)*BINSPAN+offset*(l-SAMPLESIZE); j++) {
                            sumBin[2*SAMPLESIZE-l][t] += tag[thisLociCenter+j];
                        }
                    }
                }
            }
            for (int i=0; i<SAMPLESIZE*2+1; i++) {
                totalAround += sumBin[i][t];
            }
            if (totalAround==0) {
                totalAround+=SMALLNUM;
            }
            for (int i=0; i<SAMPLESIZE*2+1; i++) {
                //normalization for KL divergence
                sumBin[i][t] /= totalAround;
            }
            signalIntensity.push_back(totalAround/loci.size());
        }
    }
    for (int t=0; t<gR.tagName.size(); t++) {
        //vector<short int> &tag = gR.regionTags[gR.tagName[t]];
        //what's wrong with query.size()?  answer:   unsigned int always render the expression to be possitive
        //int counter = loci.size();
        //FFT
        float assymetry=0;
        if (option["FFT"]=="T") {
            vector<float> tempVec;
            for (int i=0; i<SAMPLESIZE*2; i++) {
                tempVec.push_back(sumBin[i][t]);
            }
            FFT tempFFT(tempVec);
            tagNoise.push_back(tempFFT.denoise(int(SAMPLESIZE*5/8))*NOISEWEIGHT*sqrt(loci.size()));
            /* smoothing */
            vector<float> tempBin;
            for (int i=0; i<SAMPLESIZE*2; i++) {
                tempBin.push_back(tempFFT.invTrans[i].real()>0?tempFFT.invTrans[i].real():0);
            }
            tempBin.push_back(sumBin[2*SAMPLESIZE][t]);
            assymetry = testSymmety(tempBin);
            
        }
        else {
            vector<float> tempVec;
            for (int i=0; i<SAMPLESIZE*2+1; i++) {
                tempVec.push_back(sumBin[i][t]);
            }
            assymetry = testSymmety(tempVec);
        }
        float biPeak=fabsf(sumBin[SAMPLESIZE][t]*2);
        for (int i=1; i<PEAKRANGE; i++) {
            float temp = fabsf(sumBin[SAMPLESIZE-i][t]+sumBin[SAMPLESIZE+i][t]);
            if (temp>biPeak) {
                biPeak=temp;
            }
        }
        tagBiPeak.push_back(biPeak*BIPEAKWEIGHT);
        tagSymmetry.push_back(assymetry*SYMMETRYWEIGHT*sqrt(loci.size()));
        //cout<<sumMotif<<sumBin<<endl;
        //cout<<tagName<<" "<<query<<" "<<sumMotif<<" "<<" "<<(2*BINSPAN+query.size()-1)<<" "<<sumMotif/counter/(2*BINSPAN+query.size()-1)<<endl;        
    }
    if (ifDraw) {
        string fileName = option["outdir"] + "/motif.dist";
        ofstream distFile(fileName.c_str(),ios::app);
        if (!distFile) {
            printAndExit("Error! Fail to open distFile for writing!");
        }
        drawDist(gR,distFile);
    }
}

void Motif::drawDist(genomeRegions& gR,ostream &distFile){    
    distFile<<">"<<query<<"_Tag_"<<gR.tagName<<"_Bipeak_"<<tagBiPeak<<"_symmetry_"<<tagSymmetry<<"\n";
    for (int l=0; l<2*SAMPLESIZE+1; l++) {
        int pos = (l-SAMPLESIZE)*BINSPAN+offset*(l-SAMPLESIZE);
        distFile<<pos<<"\t";
        for (int t=0; t<gR.tagName.size(); t++) {
            distFile<<sumBin[l][t]<<"\t";
        }
        distFile<<"\n";
    }
    distFile<<endl;
}

void Motif::initBin(genomeRegions &gR){
    int tagSize = gR.tagName.size();
    for (int l=0; l<2*SAMPLESIZE+1; l++) {
        sumBin[l].assign(tagSize, 0);
    }
}
//score methods

void Motif::initProb(const genomeRegions& gR,int order){
    float probBackgrd[4]={0.25,0.25,0.25,0.25};
    float tempProb = 0;
    if (order==0) {
        tempProb = 1;
        for (int i=0; i<query.size(); i++){
            tempProb *= gR.prob[0][alp2num(query[i])];        
        }
    }
    else if (order==1) {
        tempProb = gR.prob[0][alp2num(query[0])];
        for (int i=1; i<query.size(); i++){
            tempProb *= gR.prob[alp2num(query[i-1])+1][alp2num(query[i])];       
        }
    }
    else if (order==-1) {
        //motifProb = pow1(0.25, K);
        tempProb = probBackgrd[alp2num(query[0])];
        for (int i=1; i<query.size(); i++){
            tempProb *= probBackgrd[alp2num(query[i])];       
        } 
    }
    motifProb = tempProb;
}

void Motif::initPWM(){
    for (int i=0; i<query.size(); i++) {
        int count = loci.size();
        switch (query[i]) {
            case 'A':
                pwm[0][i] += count;
                break;
            case 'C':
                pwm[1][i] += count;
                break;
            case 'G':
                pwm[2][i] += count;
                break;
            case 'T':
                pwm[3][i] += count;
                break;
            default:
                break;
        }
    }
}


bool Motif::ascending(vector<char> &x,const vector<int> &traverseP,int &travIndex){
    int tempIndex = travIndex;
    for (int i=0; i<traverseP.size(); i++) {
        x[traverseP[traverseP.size()-i-1]]=tempIndex%4+1;
        tempIndex = int(tempIndex/4);
    } 
    travIndex++;
    if (travIndex>pow1(4,traverseP.size())){
        //traversal end 
        travIndex=0;
        return true;
    }
    return false;
}

//loci methods

bool Motif::writeLoci(ostream &s,genomeRegions &gR){
    //cerr<<query<<lociScore.size()<<" "<<loci.size()<<endl;
    assert(lociScore.size()==loci.size());
    string strand;
    int startP = 0;
    int endP = 0;
    string chr;
    pair<int,int> sub_dist;
    s<<">"<<query<<"_lociSize"<<loci.size()<<"\n";
    //for roc 
    vector<tempLociScore> lociScores;
    int extender=atoi(option["extend"].c_str());
    try {
        for (int i=0; i<loci.size(); i++) {
            Motif tempM(1);
            tempM.loci.clear();
            tempM.loci.push_back(loci[i]);
            if (option["mode"]=="tag") {
                tempM.testMotifTag(gR ,false);  
            }
            //for roc
            lociScores.push_back(tempLociScore(tempM,lociScore[i]));
            sub_dist = locateSubscript(gR.segmentStartPos, gR.segmentStartPos.begin(), gR.segmentStartPos.end(), loci[i]);
            if (sub_dist.first==-1) {
                continue;
            }
            //for roc
            if (option["ROC"]=="T") {
                if (sub_dist.second>extender-100&&sub_dist.second<extender+200) {
                    lociScores.back().TP=true;
                }
            }
            //
            if (loci[i]<=RegionSize/2) {
                strand = "+";
                chr = gR.segmentGenomePos[sub_dist.first].first;
                startP = gR.segmentGenomePos[sub_dist.first].second+sub_dist.second;
                endP = startP + query.size();
            }
            else if (loci[i]<=RegionSize){
                strand = "-";
                chr = gR.segmentGenomePos[sub_dist.first].first;
                endP = gR.segmentGenomePos[sub_dist.first].second + (gR.segmentStartPos[sub_dist.first]-gR.segmentStartPos[sub_dist.first-1]-sub_dist.second);
                startP = endP - query.size();
            }
            else {
                continue;
            }
            s<<chr<<"\t"<<startP<<"\t"<<endP<<"\t"<<strand<<"\t"<<lociScore[i]<<"\t";
            //put two
            /*
            for (int i=0; i<tempM.tagBiPeak.size(); i++) {
                s<<tempM.tagBiPeak[i]<<"\t";
            }
            for (int i=0; i<tempM.signalIntensity.size(); i++) {
                s<<tempM.signalIntensity[i]<<"\t";
            }
            */ 
            s<<"\n";
            
        }
        
/*        
#ifndef WRITELOCIDIST
#define WRITELOCIDIST
        if (option["mode"]=="tag") {
            for (int i=0; i<loci.size(); i++) {
                ofstream specLociFile((option["outdir"]+"/specloci.dist").c_str(),ios::app);
                ostringstream ss;
                ss<<i;
                lociScores[i].query="loci No."+ss.str();
                lociScores[i].drawDist(gR, specLociFile);
            }
        }
#endif
*/        
        
        //for roc
        if (option["ROC"]=="T") {
            // from big to small
            ofstream rocFile((option["outdir"]+"/roc.txt").c_str(),ios::app);
            for (int i=0; i<lociScores.size(); i++) {
                lociScores[i].CountOnly();
            }
            sort(lociScores.begin(), lociScores.end(), compareLoci());
            int tempIndex=0;
            int tempTP=0;
            int tempTPFP=0;
            int TPFN=gR.segmentCount;
            int TN=gR.segmentCount*5; 
            rocFile<<"only count"<<"\n";
            rocFile<<"thresh\tsensitivity\tprecision\tspecificity\n";
            for (int scorethresh=int(lociScores[0].myScore+1); scorethresh>=int(lociScores.back().myScore-1); scorethresh--) {
                while (lociScores[tempIndex].myScore>(scorethresh-1)) {
                    if (tempIndex>=lociScores.size()) {
                        break;
                    }   
                    tempTPFP++;
                    
                    assert(lociScores[tempIndex].loci.size()==1);
                    if (lociScores[tempIndex].TP) {
                        tempTP++;
                    }
                    tempIndex++;
                }
                rocFile<<scorethresh<<"\t"<<tempTP/float(TPFN)<<"\t"<<(tempTP+SMALLNUM)/(float(tempTPFP)+SMALLNUM)<<"\t"<<TN/float(TN+tempTPFP-tempTP)<<"\n";
            }
            
            if (option["mode"]=="tag") {
                // roc 2
                for (int i=0; i<lociScores.size(); i++) {
                    lociScores[i].CountAndBipeak();
                }
                sort(lociScores.begin(), lociScores.end(), compareLoci());
                tempIndex=0;
                tempTP=0;
                tempTPFP=0;
                rocFile<<"count and bipeak"<<"\n";
                rocFile<<"thresh\tsensitivity\tprecision\tspecificity\n";
                for (int scorethresh=int(lociScores[0].myScore+1); scorethresh>=int(lociScores.back().myScore-1); scorethresh--) {
                    while (lociScores[tempIndex].myScore>(scorethresh-1)) {
                        if (tempIndex>=lociScores.size()) {
                            break;
                        }   
                        tempTPFP++;
                        
                        assert(lociScores[tempIndex].loci.size()==1);
                        if (lociScores[tempIndex].TP) {
                            tempTP++;
                        }
                        tempIndex++;
                    }
                    rocFile<<scorethresh<<"\t"<<tempTP/float(TPFN)<<"\t"<<tempTP/float(tempTPFP)<<"\t"<<TN/float(TN+tempTPFP-tempTP)<<"\n";
                }
                
                
                // roc 3
                
                for (int i=0; i<lociScores.size(); i++) {
                    lociScores[i].CountAndIntensity();
                }
                sort(lociScores.begin(), lociScores.end(), compareLoci());
                tempIndex=0;
                tempTP=0;
                tempTPFP=0;
                rocFile<<"count and intensity"<<"\n";
                rocFile<<"thresh\tsensitivity\tprecision\tspecificity\n";
                for (int scorethresh=int(lociScores[0].myScore+1); scorethresh>=int(lociScores.back().myScore-1); scorethresh--) {
                    while (lociScores[tempIndex].myScore>(scorethresh-1)) {
                        if (tempIndex>=lociScores.size()) {
                            break;
                        }   
                        tempTPFP++;
                        
                        assert(lociScores[tempIndex].loci.size()==1);
                        if (lociScores[tempIndex].TP) {
                            tempTP++;
                        }
                        tempIndex++;
                    }
                    rocFile<<scorethresh<<"\t"<<tempTP/float(TPFN)<<"\t"<<tempTP/float(tempTPFP)<<"\t"<<TN/float(TN+tempTPFP-tempTP)<<"\n";
                }
                

            }
            option["ROC"]="F";
        }
    } catch (const exception &e) {
        cerr<<"Write loci file err:"<<e.what()<<endl;
        return false;
    }
    s<<endl;
    return true;
}


void Motif::initLociScore(){
    //for roc?
    lociScore.clear();
    lociScore.assign(loci.size(), 1);
}

float Motif::PearsonCorrPWM(const Motif &m,int offset,int clustersize,bool strand){
    float innerProduct = 0;
    float sigmaX = 0;
    float sigmaY = 0;
    float sigmaSquareX = 0;
    float sigmaSquareY = 0;
    int samplesize = 0;
    if (strand){
        for (int i=max(-offset, 0); i<min(K, clustersize-offset); i++) {
            for (int nucleotide=0; nucleotide<4; nucleotide++) {
                innerProduct += pwm[nucleotide][i+offset]*m.pwm[nucleotide][i];
                sigmaX += pwm[nucleotide][i+offset];
                sigmaY += m.pwm[nucleotide][i];
                sigmaSquareX += pwm[nucleotide][i+offset]*pwm[nucleotide][i+offset];
                sigmaSquareY += m.pwm[nucleotide][i]*m.pwm[nucleotide][i];
                samplesize++;
            }
        }
    }
    else {
        int pwmsize = m.pwm[0].size()-1;
        for (int i=max(-offset, 0); i<min(K, clustersize-offset); i++) {
            for (int nucleotide=0; nucleotide<4; nucleotide++) {
                innerProduct += pwm[nucleotide][i+offset]*m.pwm[3-nucleotide][pwmsize-i];
                sigmaX += pwm[nucleotide][i+offset];
                sigmaY += m.pwm[3-nucleotide][pwmsize-i];
                sigmaSquareX += pwm[nucleotide][i+offset]*pwm[nucleotide][i+offset];
                sigmaSquareY += m.pwm[3-nucleotide][pwmsize-i]*m.pwm[3-nucleotide][pwmsize-i];
                samplesize++;
            }
        }
        
    }
    
    
    return (samplesize*innerProduct-sigmaX*sigmaY)/sqrtf((samplesize*sigmaSquareX-sigmaX*sigmaX)*(samplesize*sigmaSquareY-sigmaY*sigmaY));
}



void Motif::generateIUPAC(){
    //protocol: in: PWM
    //out: query
    float THRESHOLD = 0.15;
    string tempStr("");
    sumPWM();
    char alp[4]={'A','C','G','T'};
    for (int pos = 0; pos<pwm[0].size(); pos++) {
        tempStr.push_back(' ');
        for (int base = 0; base<4; base++) {
            if (pwm[base][pos]/float(totalPWM[pos])>THRESHOLD) {
                tempStr[tempStr.size()-1] = degenerate(tempStr[tempStr.size()-1], alp[base]);
            }
        }
    }
    query = tempStr;   
}




//clustering methods


pair<float,int> Cluster::editDistance(const Motif& m){
    Cluster &cluster = (*this);
    float score = 1000;
    int optimShift = -100;
    int DISTWEIGHT =10;
    int bound = MAXSHIFT+1+cluster.pwm[0].size()-K;
    for (int i=-MAXSHIFT; i<bound; i++) {
        //cout<<i<<cluster.query<<endl;
        //offset punishment 0.5 per site
        float tempScore = 0;
        if (i<0) {
            tempScore = -0.5*i;
        }
        else if (i>=bound-MAXSHIFT) {
            tempScore = 0.5*(MAXSHIFT-bound+i+1);
            //tempScore = K-cluster.query.size()+i;
        }
        /*
        for (int j=0; j<K; j++) {
            if (j+i<0||j+i>=int(cluster.query.size())) {
                continue;
            }
            
            else {
                //tempScore += distance[alp2num(m.query[j])][alp2num(cluster.query[j+i])];
                tempScore += 1-log2f(cluster.pwm[alp2num(m.query[j])][j+i]/float(cluster.totalPWM[j+i])+1);
                
            
            } 
        }
        */
        tempScore = (1-PearsonCorrPWM(m,i,int(cluster.pwm[0].size()),true))*DISTWEIGHT;
        if (tempScore<score) {
            score = tempScore;
            optimShift = i;
        }
    }
    // discard antisenses
    bool antiBetter = false;
    //string anti(::antisense(m.query));
    for (int i=-MAXSHIFT; i<bound; i++) {
        //cout<<i<<cluster.query<<endl;
        float tempScore = 0;
        if (i<0) {
            tempScore = -0.5*i;
        }
        else if (i>=bound-MAXSHIFT) {
            tempScore = 0.5*(MAXSHIFT-bound+i+1);
            //tempScore = K-cluster.query.size()+i;
        }
        tempScore = (1-PearsonCorrPWM(m,i,int(cluster.pwm[0].size()),false))*DISTWEIGHT;
        if (tempScore<score) {
            score = tempScore;
            optimShift = i;
            antiBetter = true;
        }
    }
    if (antiBetter) {
        score = -score;
    }
    return make_pair(score, optimShift);
}



void Cluster::concatenate(const Motif& m,int index,int optimShift){
    // degenerate query
    //assert this is a cluster
    //assert(m.index!=-1);
    // expMotifs++
    // index is qualified index 
    //expMotifs.push_back(index);
    if (optimShift>=0&&optimShift<query.size()-K+1) {
        for (int j=0; j<K; j++) {
            if (query[j+optimShift]!=m.query[j]) {
                query[j+optimShift]=degenerate(m.query[j],query[j+optimShift]);
            }
        }
    }
    //extending cluster in front
    else if (optimShift<0){
        //marked as cluster
        for (int j=-optimShift; j<K; j++) {
            if (query[j+optimShift]!=m.query[j]) {
                query[j+optimShift]=degenerate(m.query[j],query[j+optimShift]);
            }
        }
        index = -1;
        string temp(m.query.substr(0,-optimShift));
        query = temp+query;
    }
    //extending cluster in back
    else {
        //marked as cluster
        //offset = 1,2,3
        int offset = -(query.size()-K-optimShift);
        for (int j=0; j<K-offset; j++) {
            if (query[j+optimShift]!=m.query[j]) {
                query[j+optimShift]=degenerate(m.query[j],query[j+optimShift]);
            }
        }
        index = -1;
        string temp(m.query.substr(K-offset,offset));
        query = query+temp;
    }
    //    cout<<query<<endl;
}
//calculate KL divergence: motif to cluster
float Cluster::tagDistrDistance(const Motif& m){
    //major distrib:cluster  3 tag's sum
    float temp = 0;
    for (int t=0; t<sumBin[0].size(); t++) {
        for (int l=0; l<SAMPLESIZE*2+1; l++) {
            temp += sumBin[l][t]*log10f((sumBin[l][t]+SMALLNUM)/(m.sumBin[l][t]+SMALLNUM));
        }
    }
    //shrinker for small loci.size
    float shrinker=600/(loci.size()+1)+600/(m.loci.size()+1)+1;
    return temp/sumBin[0].size()/shrinker;
}




void Cluster::mergeLoci(){
    // change score implicitly
    
    lociScore.clear();
    sort(loci.begin(),loci.end());
    int prevSize = loci.size();
    vector<int> newloci(loci.begin(),loci.begin()+1);
    lociScore.push_back(1);
    for (int i=1; i<loci.size(); i++) {
        assert(lociScore.size()==newloci.size());
        if (loci[i]-newloci.back()<=query.size()-K) {
            lociScore.back()++;
        }
        else {
            newloci.push_back(loci[i]);
            lociScore.push_back(1);
        }
    }
    loci.swap(newloci);
    score = score*loci.size()/float(prevSize);
}

void Cluster::trim(){
    int start = 0;
    int end = query.size()-1;
    int maxPWM = sumPWM();
    int trimMultiper = 100 - atoi(option["trim"].c_str());
    while (trivial(start)||
           totalPWM[start]<maxPWM/trimMultiper||
           ((totalPWM[start]<maxPWM*2/trimMultiper)&&oligo(start))) {
        start++;
    }
    while (trivial(end)||
           totalPWM[end]<maxPWM/trimMultiper||
           ((totalPWM[end]<maxPWM*2/trimMultiper)&&oligo(start))) {
        end--;
    }
    query = query.substr(start,end-start+1);
    for (int j=0; j<4; j++) {
        vector<int> tempVec(pwm[j].begin()+start,pwm[j].begin()+end+1);
        // cout<<*this<<endl<<tempVec<<endl;
        pwm[j].swap(tempVec);
    }
}

void Cluster::reCalSumBin(const Motif& m,const genomeRegions &gR){
    //reCalsumBin and signalIntensity
    for (int t=0; t<gR.tagName.size(); t++) {
        for (int l=0; l<2*SAMPLESIZE+1; l++) {
            sumBin[l][t]=(sumBin[l][t]*loci.size()+m.sumBin[l][t]*m.loci.size());
        }
        /*
        float totalAround = 0;
        for (int i=0; i<SAMPLESIZE*2+1; i++) {
            totalAround += sumBin[i][t];
        }
        if (totalAround==0) {
            totalAround+=SMALLNUM;
        }
        for (int i=0; i<SAMPLESIZE*2+1; i++) {
            //normalization
            sumBin[i][t] /= totalAround;
        }
        */
    }
}

bool Cluster::trivial(int pos){
    /* 
    //only one
    for (int i=0; i<4; i++) {
        if (pwm[i][pos]==totalPWM[pos])
            return true;
    }
    */
    //all same
    for (int i=0; i<4; i++) {
        //cerr<<*this<<endl<<PEUSUDOCOUNT<<endl;
        assert(totalPWM[pos]!=0);
        if (pwm[i][pos]/float(totalPWM[pos])>0.4)
            return false;
    }
    return true;
}
bool Cluster::oligo(int pos){
    //only one
    for (int i=0; i<4; i++) {
        assert(totalPWM[pos]!=0);
        if (pwm[i][pos]==totalPWM[pos])
            return true;
    }
    return false;
}

void Cluster::getExtended(const Motif &m, genomeRegions &gR, Suffix & active){
    //protocol:
    //         outlet: pwm,loci,sumbin,4*score
    // extend word
    initBin(gR);
    string originWord = m.query;
    query = m.query;
    char alps[4] = {'A','C','G','T'};
    for (int pos=0; pos<K; pos++) {
        string tempWord = originWord;
        //permutation in first places
        for (int alp=0; alp<4; alp++) {
            tempWord[pos] = alps[alp];
            Motif thisMotif(tempWord);
            thisMotif.initProb(gR, atoi(option["order"].c_str()));
            thisMotif.loci = active.locateMotif(thisMotif);
            thisMotif.calConscore(RegionSize);
            if (thisMotif.score<=0) {
                continue;
            }
            pwm[alp][pos] = thisMotif.loci.size();
            appendLoci(thisMotif);
            addProb(thisMotif, K);
        }        
    }
    
    //need normalization sumbin
    for (int t=0; t<gR.tagName.size(); t++) {
        float totalAround = 0;
        for (int i=0; i<SAMPLESIZE*2+1; i++) {
            totalAround += sumBin[i][t];
        }
        if (totalAround==0) {
            totalAround+=SMALLNUM;
        }
        for (int i=0; i<SAMPLESIZE*2+1; i++) {
            //normalization
            sumBin[i][t] /= totalAround;
        }
    }
    calConscore(RegionSize);
    testMotifTag(gR,false);
    sumOverallScore();
    
    //IUPAC annotation

}


int Motif::sumPWM(){
    totalPWM.clear();
    for (int i=0; i<pwm[0].size(); i++) {
        totalPWM.push_back(pwm[0][i]+pwm[1][i]+pwm[2][i]+pwm[3][i]);
    }
    return *max_element(totalPWM.begin(), totalPWM.end());
}

void Cluster::calPWM(const Motif& m,int optimShift){
    //    explainMotif();
    //assert this is a cluster
    //not query.size() because query.size()!=pwm[0].size()now!!

    if (optimShift>=0&&optimShift<pwm[0].size()-K+1) {
        for (int i=0; i<K; i++) {
            for (int j=0; j<4; j++) {
                pwm[j][i+optimShift] += m.pwm[j][i];
            }
            //              cout<<*this<<endl;
            //    cout<<m<<endl;
        }
    }

    //extending cluster in front
    else if (optimShift<0){
        //marked as cluster
        for (int i=-optimShift; i<K; i++) {
            for (int j=0; j<4; j++) {
                pwm[j][i+optimShift] += m.pwm[j][i];
            }
        }
        for (int i=-optimShift-1; i>=0; i--) {
            for (int j=0; j<4; j++) {
                pwm[j].insert(pwm[j].begin(),m.pwm[j][i]) ;
            }
        }
    }
    //extending cluster in back
    else {
        //marked as cluster
        //offset = 1,2,3..
        int offset = -(pwm[0].size()-K-optimShift);
        for (int i=0; i<K-offset; i++) {
            for (int j=0; j<4; j++) {
                pwm[j][i+optimShift] += m.pwm[j][i];
            }
        }
        for (int i=K-offset; i<K; i++) {
            for (int j=0; j<4; j++) {
                pwm[j].push_back(m.pwm[j][i]) ;
            }
        }
        
    }
    sumPWM();
    //    cout<<query<<endl;
}




void Motif::sumOverallScore(){
    float STDScore = 0;
    if (tagBiPeak.size()>0) {
        for(int i=0;i<tagBiPeak.size();i++){
            STDScore += tagBiPeak[i];
        }
        STDScore /= tagBiPeak.size();
    }
    float noiseScore = 0;
    if (tagNoise.size()>0) {
        for(int i=0;i<tagNoise.size();i++){
            noiseScore += tagNoise[i];
        }
        noiseScore /= tagNoise.size();
    }
    float symmetryScore = 0;
    if (tagSymmetry.size()>0) {
        for(int i=0;i<tagSymmetry.size();i++){
            symmetryScore += tagSymmetry[ i];
        }
        symmetryScore /= tagSymmetry.size();
    }
    overallScore = score+STDScore-noiseScore-symmetryScore+0.0001;
}

void Cluster::mergeProb(const Motif& m){
    score =(score*loci.size()+m.score*m.loci.size())/(loci.size()+m.loci.size());
    motifProb = (loci.size()+m.loci.size())/RegionSize/score;
    return;
}


float Motif::sumTagScore(){
    float tempScore = 0;
    if (tagBiPeak.size()==0) {
        return 0;
    }
    for(int i=0;i<tagBiPeak.size();i++){
        tempScore += tagBiPeak[i];
    }
    return tempScore;
}

void Motif::printMotif(ostream &s){
    s.precision(4);
    s<<std::setw(25)<<setiosflags(std::ios::left)<<query<<"\t"<<pvalue()<<"\t"<<score<<"\t";
    for (int i=0; i<tagBiPeak.size(); i++) {
        s<<tagBiPeak[i]<<"\t";
    }
    for (int i=0; i<tagSymmetry.size(); i++) {
        s<<tagSymmetry[i]<<"\t";
    }
    if (option["FFT"]=="T") {
        for (int i=0; i<tagNoise.size(); i++) {
            s<<tagNoise[i]<<"\t";
        }
    }
    /*
    for (int i=0; i<signalIntensity.size(); i++) {
        s<<signalIntensity[i]<<"\t";
    }
    */
    s<<loci.size()<<endl;
}

ostream &operator<<( ostream &s, Motif &motif ){
    s<<">"<<motif.query<<"_Conscore_"<<motif.score;
    if (motif.tagBiPeak.size()!=0) {
        s<<"_tagBiPeak_"<<motif.tagBiPeak;
    }
    if (motif.tagSymmetry.size()!=0) {
        s<<"_tagSymmetry_"<<motif.tagSymmetry;
    }
    if (motif.tagNoise.size()!=0) {
        s<<"_tagNoise_"<<motif.tagNoise;
    }
    s<<"\n";
    for (int i=0; i<4; i++) {
        for (int j=0; j<motif.pwm[i].size(); j++) {
            s<<motif.pwm[i][j]/float(motif.pwm[0][j]+motif.pwm[1][j]+motif.pwm[2][j]+motif.pwm[3][j])<<'\t';
        }
        s<<endl;
    }
    return s;
}

