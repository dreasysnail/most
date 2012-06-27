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

int GenomeSize;
/*
int Motif::distance[4][15] = {  {0,1,1,1,0,0,0,1,1,1,1,0,0,0,0},
                                {1,0,1,1,0,1,1,0,0,1,0,1,0,0,0}, 
                                {1,1,0,1,1,0,1,0,1,0,0,0,1,0,0},
                                {1,1,1,0,1,1,0,1,0,0,0,0,0,1,0}};
*/
///*
int Motif::distance[4][15] = {  {0,6,6,6,1,1,1,5,5,5,4,2,2,2,3},
                                {6,0,6,6,1,5,5,1,1,5,2,4,2,2,3}, 
                                {6,6,0,6,5,1,5,1,5,1,2,2,4,2,3},
                                {6,6,6,0,5,5,1,5,1,1,2,2,2,4,3}};
//*/
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
    //cout<<GenomeSize<<endl;
    for (int i=0; i<GenomeSize; i++) {
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



void Motif::testMotifTag(genomeRegions &gR,const string& outPutDir,bool ifDraw){
    assert(loci.size()!=0);
    //only need loci and tagSeq
    signif.clear();
    for (int j=0; j<gR.tagName.size(); j++) {
        vector<short int> &tag = gR.regionTags[gR.tagName[j]];
        string tagName = gR.tagName[j];
        float sumMotif=0;
        float sumAround[2*SAMPLESIZE]={0};
        //what's wrong with query.size()?  answer:   unsigned int always render the expression to be possitive
        int queryLength=query.size();
        int counter = loci.size();
        for (int i=0; i<loci.size(); i++) {
            if (loci[i]<(2*SAMPLESIZE+1)*dist+SAMPLESIZE*queryLength+offset||loci[i]>GenomeSize-(2*SAMPLESIZE+1)*dist-(SAMPLESIZE+1)*queryLength-1-offset) {
                counter--;
                continue;
            }
            int thisLociCenter = mapLoci(loci[i],gR.segmentStartPos);
            
            for (int j=-dist; j<dist+queryLength; j++) {
                sumMotif += tag[thisLociCenter+j];
                
            }
            //     for (int j=-3*dist-K-offset1; j<-dist-offset1; j++) sumShoulder += tag[thisLociCenter+j];
            //     for (int j=dist+K+offset1; j<3*dist+2*K+offset1; j++) sumShoulder += tag[thisLociCenter+j];
            for (int l=0; l<SAMPLESIZE; l++) {
                for (int j=-(2*l+3)*dist-(l+1)*queryLength-offset*(l+1); j<-(2*l+1)*dist-l*queryLength-offset*(l+1); j++) {
                    sumAround[l] += tag[thisLociCenter+j];
                }
                for (int j=(2*l+1)*dist+(l+1)*queryLength+offset*(l+1); j<(2*l+3)*dist+(l+2)*queryLength+offset*(l+1); j++) {
                    sumAround[l+SAMPLESIZE] += tag[thisLociCenter+j];
                }
            }
        }
        float totalAround=0;
        float variance = 0;
        for (int i=0; i<SAMPLESIZE*2; i++) {
            totalAround += sumAround[i];
        }
        float average = (totalAround+sumMotif)/(SAMPLESIZE*2+1);
        for (int i=0; i<SAMPLESIZE*2; i++) {
            variance += (sumAround[i]-average)*(sumAround[i]-average);
        }
        variance += (sumMotif-average)*(sumMotif-average);
        variance = sqrt(variance/SAMPLESIZE*2)+0.000001;
        //vector<float> temp(sumAround,sumAround+2*SAMPLESIZE);
        //cout<<sumMotif<<temp<<average<<"var "<<variance<<" Counter"<<counter<<" height"<<logf(sumMotif/float(counter))<<endl;
        //signif.push_back((sumMotif-average)/variance*logf(sumMotif/float(counter))) ;
        signif.push_back((sumMotif-average)/variance);
        //cout<<sumMotif<<sumAround<<endl;
        //cout<<tagName<<" "<<query<<" "<<sumMotif<<" "<<" "<<(2*dist+query.size()-1)<<" "<<sumMotif/counter/(2*dist+query.size()-1)<<endl;
        if (ifDraw) {
            string fileName = outPutDir + "/motif.dist";
            ofstream distFile(fileName.c_str(),ios::app);
            if (!distFile) {
                string errorInfo = "Error! Fail to open distFile for writing!";
                printAndExit(errorInfo);
            }
            distFile<<">"<<query<<"_Conscore_"<<score<<"_Tag_"<<tagName<<"_AVG_"<<sumMotif/counter/(2*dist+query.size()-1)<<"_Sig_"<<signif<<"\n";
            for (int l=SAMPLESIZE-1;l>=0 ;l--) {
                int pos = -(2*l+2)*dist-(l+1)*query.size()-offset;
                //cout<<pos<<endl;
                distFile<<pos<<"\t"<<sumAround[l]/4/(dist*2+queryLength)<<"\n";
            }
            distFile<<0<<"\t"<<sumMotif/4/(dist*2+queryLength)<<"\n";
            for (int l=0; l<SAMPLESIZE; l++) {
                distFile<<(2*l+2)*dist+(l+1)*query.size()+offset<<"\t"<<sumAround[l+SAMPLESIZE]/4/(dist*2+queryLength)<<"\n";
            }
            
            distFile<<endl;
        }        
    }
    
        

}


void Motif::initProb(const genomeRegions& gR,int order){
    float probBackgrd[4]={0.2,0.3,0.3,0.2};
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

//clustering methods



pair<int,int> Cluster::editDistance(const Motif& m){
    Cluster &cluster = (*this);
    int score = 1000;
    int optimShift = -100;
    int bound = SHIFT+1+cluster.query.size()-K;
    for (int i=-SHIFT; i<bound; i++) {
        //cout<<i<<cluster.query<<endl;
        //gap punish 2
        int tempScore = 0;
        if (i<0) {
            tempScore = -2*i;
        }
        else if (i>=bound-SHIFT) {
            tempScore = 2*(SHIFT-bound+i+1);
            //tempScore = K-cluster.query.size()+i;
        }
        for (int j=0; j<K; j++) {
            if (j+i<0||j+i>=int(cluster.query.size())) {
                continue;
            }
            else {
                tempScore += distance[alp2num(m.query[j])][alp2num(cluster.query[j+i])];
            } 
        }
        if (tempScore<score) {
            score = tempScore;
            optimShift = i;
        }
    }
    // discard antisenses
    bool antiBetter = false;
    string anti(::antisense(m.query));
    for (int i=-SHIFT; i<bound; i++) {
        //cout<<i<<cluster.query<<endl;
        int tempScore = 0;
        if (i<0) {
            tempScore = -2*i;
        }
        else if (i>=bound-SHIFT) {
            tempScore = 2*(SHIFT-bound+i+1);
            //tempScore = K-cluster.query.size()+i;
        }
        for (int j=0; j<K; j++) {
            if (j+i<0||j+i>=cluster.query.size()) {
                continue;
            }
            else {
                tempScore += distance[alp2num(anti[j])][alp2num(cluster.query[j+i])];
            } 
        }
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

float Cluster::histoneDistrDistance(const Motif& m){
    Cluster &cluster = (*this);
    float signifDist = 0;
    for (int i=0; i<signif.size(); i++) {
        float temp = (m.signif[i]-cluster.signif[i])/(cluster.signif[i]+m.signif[i]);
        if (temp>0) {
            signifDist+=temp;
        }
        else {
            signifDist-=temp;
        }
    }
    return signifDist;
}

void Cluster::trim(){
    int start = 0;
    int end = query.size()-1;
    int maxPWM = sumPWM();
    while (unif(start)||totalPWM[start]<maxPWM/100) {
        start++;
    }
    while (unif(start)||totalPWM[end]<maxPWM/100) {
        end--;
    }
    query = query.substr(start,end-start+1);
    for (int j=0; j<4; j++) {
        vector<int> tempVec(pwm[j].begin()+start,pwm[j].begin()+end+1);
        // cout<<*this<<endl<<tempVec<<endl;
        pwm[j].swap(tempVec);
    }
}
bool Cluster::unif(int pos){
    for (int i=0; i<4; i++) {
        if (pwm[i][pos]/totalPWM[pos]<0.2||pwm[i][pos]/totalPWM[pos]>0.3)
            return false;
    }
    return true;
}

int Cluster::sumPWM(){
    totalPWM.clear();
    for (int i=0; i<pwm[0].size(); i++) {
        totalPWM.push_back(pwm[0][i]+pwm[0][i]+pwm[0][i]+pwm[0][i]);
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

    //    cout<<query<<endl;
}






void Motif::sumOverallScore(){
    overallScore = score*score;
    float tempScore = 0;
    if (signif.size()==0) {
        return;
    }
    for(int i=0;i<signif.size();i++){
        tempScore += fabs(signif[i]);
    }
    overallScore *= tempScore+0.0001;
}

ostream &operator<<( ostream &s, const Motif &motif ){
    if (motif.signif.size()!=0) {
        s<<">"<<motif.query<<"_Conscore_"<<motif.score<<"_Signif_"<<motif.signif;
    }
    else {
        s<<">"<<motif.query<<"_Conscore_"<<motif.score<<"\n";
    }
    for (int i=0; i<4; i++) {
        for (int j=0; j<motif.pwm[i].size(); j++) {
            s<<motif.pwm[i][j]<<'\t';
        }
        s<<endl;
    }
    return s;
}