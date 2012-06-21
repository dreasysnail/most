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
int Motif::motif[K_5] = {0};
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
void Motif::fillMotif(const vector<Motif>& allmotifs){
    //explainMotif();
    vector<int>::iterator it;
    //float motifProb;
    probThresh = 0;
    for (it=expMotifs.begin(); it!=expMotifs.end();it++ ) {
        //add generalization's prob
        probThresh += allmotifs[(*it)].probThresh;
    }
    calConscore(GenomeSize);
    
}
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






void Motif::testMotifTag(genomeRegions &gR,const string& outPutDir){
    signif.clear();
    for (int j=0; j<gR.tagName.size(); j++) {
        vector<short int> &tag = gR.genomeTags[gR.tagName[j]];
        string tagName = gR.tagName[j];
        float sumMotif=0;
        float sumAround[2*SAMPLESIZE]={0};
        //what's wrong with query.size()?  answer:   unsigned int always render the expression to be possitive
        int queryLength=query.size();
        //float sumShoulder=0;
        int counter = loci.size();
        for (int i=0; i<loci.size(); i++) {
            if (loci[i]<(2*SAMPLESIZE+1)*dist+SAMPLESIZE*queryLength+offset||loci[i]>tag.size()-(2*SAMPLESIZE+1)*dist-(SAMPLESIZE+1)*queryLength-1-offset) {
                counter--;
                continue;
            }
            
            
            for (int j=-dist; j<dist+queryLength; j++) {
                sumMotif += tag[loci[i]+j];
                
            }
            //     for (int j=-3*dist-K-offset1; j<-dist-offset1; j++) sumShoulder += tag[loci[i]+j];
            //     for (int j=dist+K+offset1; j<3*dist+2*K+offset1; j++) sumShoulder += tag[loci[i]+j];
            for (int l=0; l<SAMPLESIZE; l++) {
                for (int j=-(2*l+3)*dist-(l+1)*queryLength-offset; j<-(2*l+1)*dist-l*queryLength-offset; j++) {
                    sumAround[l] += tag[loci[i]+j];
                }
                for (int j=(2*l+1)*dist+(l+1)*queryLength+offset; j<(2*l+3)*dist+(l+2)*queryLength+offset; j++) {
                    sumAround[l+SAMPLESIZE] += tag[loci[i]+j];
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
        variance = sqrt(variance/SAMPLESIZE*2);
        //vector<float> temp(sumAround,sumAround+2*SAMPLESIZE);
        //cout<<sumMotif<<temp<<average<<"var "<<variance<<" Counter"<<counter<<" height"<<logf(sumMotif/float(counter))<<endl;
        //signif.push_back((sumMotif-average)/variance*logf(sumMotif/float(counter))) ;
        signif.push_back((sumMotif-average)/variance);
        if (signif.back()<0) {
            signif.back()=-signif.back();
        }
        //cout<<sumMotif<<sumAround<<endl;
        //cout<<tagName<<" "<<query<<" "<<sumMotif<<" "<<" "<<(2*dist+query.size()-1)<<" "<<sumMotif/counter/(2*dist+query.size()-1)<<endl;
        
#ifdef DRAW_MOTIF
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
            distFile<<pos<<"\t"<<sumAround[l]<<"\n";
        }
        distFile<<0<<"\t"<<sumMotif<<"\n";
        for (int l=0; l<SAMPLESIZE; l++) {
            distFile<<(2*l+2)*dist+(l+1)*query.size()+offset<<"\t"<<sumAround[l+SAMPLESIZE]<<"\n";
        }
        
        distFile<<endl;
#endif

        
    }
    
        

}


void Motif::initProb(const genomeRegions& gR,int order){
    float probtemp[4]={0.2,0.3,0.3,0.2};
    if (order==0) {
        float motifProb = 1;
        for (int i=0; i<query.size(); i++){
            motifProb *= gR.prob[0][alp2num(query[i])];        
        }
        probThresh = motifProb;  
    }
    else if (order==1) {
        float motifProb = gR.prob[0][alp2num(query[0])];
        for (int i=1; i<query.size(); i++){
            motifProb *= gR.prob[alp2num(query[i-1])+1][alp2num(query[i])];       
        }
        probThresh = motifProb;  
    }
    else if (order==-1) {
        //probThresh = pow1(0.25, K);
        float motifProb = probtemp[alp2num(query[0])];
        for (int i=1; i<query.size(); i++){
            motifProb *= probtemp[alp2num(query[i])];       
        }
        probThresh = motifProb;  
    }
      
}

void Motif::initPWM(){
    for (int i=0; i<query.size(); i++) {            
        switch (query[i]) {
            case 'A':
                pwm[0][i] += motif[index];
                break;
            case 'C':
                pwm[1][i] += motif[index];
                break;
            case 'G':
                pwm[2][i] += motif[index];
                break;
            case 'T':
                pwm[3][i] += motif[index];
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
pair<int,int> Motif::editDistance(const Motif& cluster){
    int score = 1000;
    int optimShift = -100;
    int bound = SHIFT+1+cluster.query.size()-K;
    for (int i=-SHIFT; i<bound; i++) {
        //cout<<i<<cluster.query<<endl;
        //gap punish 1
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
                tempScore += distance[alp2num(query[j])][alp2num(cluster.query[j+i])];
            } 
        }
        if (tempScore<score) {
            score = tempScore;
            optimShift = i;
        }
    }
    // discard antisenses
    bool antiBetter = false;
    string anti(::antisense(query));
    for (int i=-SHIFT; i<bound; i++) {
        //cout<<i<<cluster.query<<endl;
        int tempScore = 0;
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
            antiBetter = true;
        }
    }
    if (antiBetter) {
        score = - score;
    }
    return make_pair(score, optimShift);
}



void Motif::concatenate(const Motif& m,int index,int optimShift){
    //assert this is a cluster
    assert(m.index!=-1);
    // expMotifs++
    // index is qualified index 
    expMotifs.push_back(index);
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

void Motif::trim(){
    int start = 0;
    int end = query.size()-1;
    while (query[start]=='N') {
        start++;
    }
    while (query[end]=='N') {
        end--;
    }
    query = query.substr(start,end-start+1);
    for (int j=0; j<4; j++) {
        vector<int> tempVec(pwm[j].begin()+start,pwm[j].begin()+end+1);
        // cout<<*this<<endl<<tempVec<<endl;
        pwm[j].swap(tempVec);
    }
}

void Motif::calPWM(const Motif& m,int optimShift){
    //    explainMotif();
    //assert this is a cluster
    //not query.size() because query.size()!=pwm[0].size()now!!
    assert(m.index!=-1);
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
        //offset = 1,2,3
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






void Motif::sumScore(){
    overallScore = score*score;
    float tempScore = 0;
    for(int i=0;i<signif.size();i++){
        tempScore += signif[i];
    }
    overallScore *= tempScore/10.0;
}

ostream &operator<<( ostream &s, const Motif &motif ){
    s<<">"<<motif.query<<"_Conscore_"<<motif.score<<"_Signif_"<<motif.signif;
    for (int i=0; i<4; i++) {
        for (int j=0; j<motif.pwm[i].size(); j++) {
            s<<motif.pwm[i][j]<<'\t';
        }
        s<<endl;
    }
    return s;
}