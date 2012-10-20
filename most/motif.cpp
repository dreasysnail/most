//
//  motif.cpp
//  MOST
//
//  Created by zhang yizhe on 12-5-20.
//  Copyright (c) 2012年 SJTU. All rights reserved.
//


#include "motif.h"
#include "bwt.h"

using namespace std;

int RegionSize;

void Motif::locateMotif(const char T[]){
    //implement1 low efficiency
    bool flag = true;
    loci.clear();
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
    return EXTEND_BOUND*(no_dist.first*2+1)+genomePos;   
}

void Motif::CalSumBin(genomeRegions &gR){
    //protocol: in: loci
    //          out:normalized sumbin
    initBin(gR);
    assert(loci.size()!=0);
    for (int t=0; t<gR.tagName.size(); t++) {
        vector<short int> &tag = gR.regionTags[gR.tagName[t]];
        string tagName = gR.tagName[t];
        float totalAround=0;
        for (int i=0; i<loci.size(); i++) {
            if (loci[i]<SAMPLESIZE*BIN_SPAN+SAMPLESIZE*OFF_SET+1||loci[i]>RegionSize-(SAMPLESIZE+1)*BIN_SPAN-SAMPLESIZE*OFF_SET+1) {
                //counter--;
                continue;
            }
            int thisLociCenter = mapLoci(loci[i],gR.segmentStartPos);
            //cerr<<tag[thisLociCenter]<<" "<<thisLociCenter<<endl;
            int tempsize = gR.segmentStartPos.size()/2;
            //+ direction
            if (loci[i]<gR.segmentStartPos[tempsize]) {
                for (int l=0; l<2*SAMPLESIZE+1; l++) {
                    for (int j=(l-SAMPLESIZE)*BIN_SPAN+OFF_SET*(l-SAMPLESIZE); j<(l+1-SAMPLESIZE)*BIN_SPAN+OFF_SET*(l-SAMPLESIZE); j++) {
                        sumBin[l][t] += tag[thisLociCenter+j];
                    }
                }
            }
            //- direction
            else {
                for (int l=0; l<2*SAMPLESIZE+1; l++) {
                    for (int j=(l-SAMPLESIZE)*BIN_SPAN+OFF_SET*(l-SAMPLESIZE); j<(l+1-SAMPLESIZE)*BIN_SPAN+OFF_SET*(l-SAMPLESIZE); j++) {
                        sumBin[2*SAMPLESIZE-l][t] += tag[thisLociCenter+j];
                    }
                }
            }
        }
        //normalized sumbin
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
        //signalIntensity.push_back(totalAround/loci.size());
    }
}


void Motif::testMotifTag(genomeRegions &gR,bool ifDraw){
    if (option["mode"]!="tag")  return;
    //protocol: input:   sumbin
    //          outlet: noise bipeak asymmetry

    //only need loci and tagSeq
    tagBiPeak.clear();
    tagSymmetry.clear();
    tagNoise.clear();
    
    
    //if not initiated
    if (index != -1) {
        CalSumBin(gR);
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
            tagNoise.push_back(tempFFT.denoise(int(SAMPLESIZE*5/8))*NOISE_WEIGHT*sqrt(loci.size()));
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
        for (int i=1; i<PEAK_RANGE; i++) {
            float temp = fabsf(sumBin[SAMPLESIZE-i][t]+sumBin[SAMPLESIZE+i][t]);
            if (temp>biPeak) {
                biPeak=temp;
            }
        }
        tagBiPeak.push_back(biPeak*BIPEAK_WEIGHT);
        tagSymmetry.push_back(assymetry*SYMMETRY_WEIGHT*sqrt(loci.size()));
        //cout<<sumMotif<<sumBin<<endl;
        //cout<<tagName<<" "<<query<<" "<<sumMotif<<" "<<" "<<(2*BIN_SPAN+query.size()-1)<<" "<<sumMotif/counter/(2*BIN_SPAN+query.size()-1)<<endl;        
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
        int pos = (l-SAMPLESIZE)*BIN_SPAN+OFF_SET*(l-SAMPLESIZE);
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

void Motif::initPFM(){
    //protocol : in  loci
    //           out pfm
    for (int i=0; i<query.size(); i++) {
        int count = loci.size();
        switch (query[i]) {
            case 'A':
                pfm[0][i] += count;
                break;
            case 'C':
                pfm[1][i] += count;
                break;
            case 'G':
                pfm[2][i] += count;
                break;
            case 'T':
                pfm[3][i] += count;
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
    assert(lociScore.size()==loci.size());
    string strand;
    int startP = 0;
    int endP = 0;
    string chr;
    pair<int,int> sub_dist;
    s<<">"<<query<<"_lociSize"<<loci.size()<<"\n";
    //for roc 
    vector<tempLociScore> lociScores;
    int extender=atoi(option["flanking"].c_str());
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
            s<<"\n";
            
        }
        
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

float Motif::PearsonCorrPFM(const Motif &m,int offset,int clustersize,bool strand){
    float innerProduct = 0;
    float sigmaX = 0;
    float sigmaY = 0;
    float sigmaSquareX = 0;
    float sigmaSquareY = 0;
    int samplesize = 0;
    if (strand){
        for (int i=max(-offset, 0); i<min(K, clustersize-offset); i++) {
            for (int nucleotide=0; nucleotide<4; nucleotide++) {
                innerProduct += pfm[nucleotide][i+offset]*m.pfm[nucleotide][i];
                sigmaX += pfm[nucleotide][i+offset];
                sigmaY += m.pfm[nucleotide][i];
                sigmaSquareX += pfm[nucleotide][i+offset]*pfm[nucleotide][i+offset];
                sigmaSquareY += m.pfm[nucleotide][i]*m.pfm[nucleotide][i];
                samplesize++;
            }
        }
    }
    else {
        int pfmsize = m.pfm[0].size()-1;
        for (int i=max(-offset, 0); i<min(K, clustersize-offset); i++) {
            for (int nucleotide=0; nucleotide<4; nucleotide++) {
                innerProduct += pfm[nucleotide][i+offset]*m.pfm[3-nucleotide][pfmsize-i];
                sigmaX += pfm[nucleotide][i+offset];
                sigmaY += m.pfm[3-nucleotide][pfmsize-i];
                sigmaSquareX += pfm[nucleotide][i+offset]*pfm[nucleotide][i+offset];
                sigmaSquareY += m.pfm[3-nucleotide][pfmsize-i]*m.pfm[3-nucleotide][pfmsize-i];
                samplesize++;
            }
        }
        
    }
    
    
    return (samplesize*innerProduct-sigmaX*sigmaY)/sqrtf((samplesize*sigmaSquareX-sigmaX*sigmaX)*(samplesize*sigmaSquareY-sigmaY*sigmaY));
}



void Motif::generateIUPAC(){
    //protocol: in: PFM
    //out: query
    float THRESHOLD = 0.15;
    string tempStr("");
    sumPFM();
    char alp[4]={'A','C','G','T'};
    for (int pos = 0; pos<pfm[0].size(); pos++) {
        tempStr.push_back(' ');
        for (int base = 0; base<4; base++) {
            if (pfm[base][pos]/float(totalPFM[pos])>THRESHOLD) {
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
    int bound = MAX_SHIFT+1+cluster.pfm[0].size()-K;
    for (int i=-MAX_SHIFT; i<bound; i++) {
        //cout<<i<<cluster.query<<endl;
        //offset punishment 0.5 per site
        float tempScore = 0;
        if (i<0) {
            tempScore = -0.5*i;
        }
        else if (i>=bound-MAX_SHIFT) {
            tempScore = 0.5*(MAX_SHIFT-bound+i+1);
            //tempScore = K-cluster.query.size()+i;
        }
        tempScore += (1-PearsonCorrPFM(m,i,int(cluster.pfm[0].size()),true))*DISTWEIGHT;
        if (tempScore<score) {
            score = tempScore;
            optimShift = i;
        }
    }
    // discard antisenses
    bool antiBetter = false;
    //string anti(::antisense(m.query));
    for (int i=-MAX_SHIFT; i<bound; i++) {

        //cout<<i<<cluster.query<<endl;
        float tempScore = 0;
        if (i<0) {
            tempScore = -0.5*i;
        }
        else if (i>=bound-MAX_SHIFT) {
            tempScore = 0.5*(MAX_SHIFT-bound+i+1);
            //tempScore = K-cluster.query.size()+i;
        }
        tempScore += (1-PearsonCorrPFM(m,i,int(cluster.pfm[0].size()),false))*DISTWEIGHT;
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
        //OFF_SET = 1,2,3
        int OFF_SET = -(query.size()-K-optimShift);
        for (int j=0; j<K-OFF_SET; j++) {
            if (query[j+optimShift]!=m.query[j]) {
                query[j+optimShift]=degenerate(m.query[j],query[j+optimShift]);
            }
        }
        index = -1;
        string temp(m.query.substr(K-OFF_SET,OFF_SET));
        query = query+temp;
    }
    //    cout<<query<<endl;
}
//calculate KL divergence: motif to cluster
float Cluster::tagDistrDistance(const Motif& m){
    //major distrib:cluster  3 tag's sum

    float temp = 0;
        
    for (int t=0; t<sumBin[0].size(); t++) {
        for (int l=0; l<2*SAMPLESIZE+1; l++) {
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
    //int prevSize = loci.size();
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
    //score = score*loci.size()/float(prevSize);
}

void Cluster::trim(){
    int start = 0;
    int end = query.size()-1;
    int maxPFM = sumPFM();
    int trimMultiper = 100 - atoi(option["trim"].c_str());
    while (trivial(start)||
           totalPFM[start]<maxPFM/trimMultiper||
           ((totalPFM[start]<maxPFM*2/trimMultiper)&&oligo(start))) {
        start++;
    }
    while (trivial(end)||
           totalPFM[end]<maxPFM/trimMultiper||
           ((totalPFM[end]<maxPFM*2/trimMultiper)&&oligo(start))) {
        end--;
    }
    query = query.substr(start,end-start+1);
    for (int j=0; j<4; j++) {
        vector<int> tempVec(pfm[j].begin()+start,pfm[j].begin()+end+1);
        // cout<<*this<<endl<<tempVec<<endl;
        pfm[j].swap(tempVec);
    }
}

void Cluster::reCalSumBin(const Motif& m,const genomeRegions &gR){
    //protocol: in loci
    //reCalsumBin
    for (int t=0; t<gR.tagName.size(); t++) {
        for (int l=0; l<2*SAMPLESIZE+1; l++) {
            sumBin[l][t]=(sumBin[l][t]*loci.size()+m.sumBin[l][t]*m.loci.size())/(loci.size()+m.loci.size());
        }
    }
}

bool Cluster::trivial(int pos){
    return query[pos]=='N';
}

bool Cluster::oligo(int pos){
    //only one
    for (int i=0; i<4; i++) {
        assert(totalPFM[pos]!=0);
        if (pfm[i][pos]>=totalPFM[pos]*0.95)
            return true;
    }
    return false;
}

void Cluster::getExtended(const Motif &m, genomeRegions &gR, Suffix & active){
    //protocol:
    //         outlet: pfm,loci,sumbin,4*score
    // extend word
    initBin(gR);
    string originWord = m.query;
    query = m.query;
    Members[0].query = m.query;
    char alps[4] = {'A','C','G','T'};
    for (int pos=0; pos<K; pos++) {
        string tempWord = originWord;
        //permutation in first places
        for (int alp=0; alp<4; alp++) {
            //for each word
            tempWord[pos] = alps[alp];
            Motif thisMotif(tempWord);
            if (originWord[pos]==alps[alp]) {
                thisMotif = m;
            }
            else {
                thisMotif.initProb(gR, atoi(option["order"].c_str()));
                thisMotif.loci = active.locateMotif(thisMotif);
                thisMotif.calConscore(RegionSize);
                if (thisMotif.score<=0) {
                    continue;
                }                
            }
            pfm[alp][pos] = thisMotif.loci.size();
            appendLoci(thisMotif);
            addProb(thisMotif, K);
        }        
    }
    
    //need normalization sumbin
    CalSumBin(gR);
    calConscore(RegionSize);
    testMotifTag(gR,false);
    sumOverallScore();
    
    //IUPAC annotation

}

void Cluster::printMember(ostream &clusterFile){
    
    vector<int> offset1,offset2;
    
    int tempOffset = 0;
    for (int i=0; i<Members.size(); i++) {
        tempOffset+=(Members[i].shift<0?(-Members[i].shift):0);
        offset1.push_back(tempOffset);
    }
    offset2.push_back(offset1.back());
    for (int i=0; i<Members.size(); i++) {
        offset2.push_back(offset1.back()-offset1[i]);
    }
    
    clusterFile<<" "<<query<<endl;
    for (int i=0; i<Members.size(); i++) {
        clusterFile<<setw(8)<<"MEMBER\t"<<setw(3)<<i<<"\t";
        for (int counter=0; counter<Members[i].shift+offset2[i]; counter++) {
            clusterFile<<" ";
        }
        clusterFile<<" "<<Members[i].query;
        clusterFile<<setw(8)<<"Bipeak\t"<<Members[i].tagBiPeak<<setw(8)<<"symmetry\t"<<Members[i].tagSymmetry<<setw(8)<<"dist\t"<<Members[i].dist<<setw(8)<<"shift\t"<<Members[i].shift<<"\n";
    }
    clusterFile<<"\n";
}


int Motif::sumPFM(){
    totalPFM.clear();
    for (int i=0; i<pfm[0].size(); i++) {
        totalPFM.push_back(pfm[0][i]+pfm[1][i]+pfm[2][i]+pfm[3][i]);
    }
    return *max_element(totalPFM.begin(), totalPFM.end());
}

void Cluster::calPFM(const Motif& m,int optimShift){
    //    explainMotif();
    //assert this is a cluster
    //not query.size() because query.size()!=pfm[0].size()now!!

    if (optimShift>=0&&optimShift<pfm[0].size()-K+1) {
        for (int i=0; i<K; i++) {
            for (int j=0; j<4; j++) {
                pfm[j][i+optimShift] += m.pfm[j][i];
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
                pfm[j][i+optimShift] += m.pfm[j][i];
            }
        }
        for (int i=-optimShift-1; i>=0; i--) {
            for (int j=0; j<4; j++) {
                pfm[j].insert(pfm[j].begin(),m.pfm[j][i]) ;
            }
        }
    }
    //extending cluster in back
    else {
        //marked as cluster
        //offset = 1,2,3..
        int offset = -(pfm[0].size()-K-optimShift);
        for (int i=0; i<K-offset; i++) {
            for (int j=0; j<4; j++) {
                pfm[j][i+optimShift] += m.pfm[j][i];
            }
        }
        for (int i=K-offset; i<K; i++) {
            for (int j=0; j<4; j++) {
                pfm[j].push_back(m.pfm[j][i]) ;
            }
        }
        
    }
    sumPFM();
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
    //merging score
    //score = avg.score of members weighted by loci.size()
    score =(score*loci.size()+m.score*m.loci.size())/(loci.size()+m.loci.size());
    //motifProb = (loci.size()+m.loci.size())/RegionSize/score;
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
    s<<loci.size()<<endl;
}

ostream &operator<<( ostream &s, Motif &motif ){
    s<<">"<<motif.query;
    /*
    s<<"_Conscore_"<<motif.score;
    if (motif.tagBiPeak.size()!=0) {
        s<<"_tagBiPeak_"<<motif.tagBiPeak;
    }
    if (motif.tagSymmetry.size()!=0) {
        s<<"_tagSymmetry_"<<motif.tagSymmetry;
    }
    if (motif.tagNoise.size()!=0) {
        s<<"_tagNoise_"<<motif.tagNoise;
    }
    */
    s<<"\n";
    for (int i=0; i<4; i++) {
        for (int j=0; j<motif.pfm[i].size(); j++) {
            s<<motif.pfm[i][j]/float(motif.pfm[0][j]+motif.pfm[1][j]+motif.pfm[2][j]+motif.pfm[3][j])<<'\t';

        }
        s<<endl;
    }
    return s;
}

