//
//  main.cpp
//  bwt
//
//  Created by icarus on 12-5-10.
//  Copyright 2012年 sjtu. All rights reserved.
//

#include <iostream>
#include <string>
#include <queue>
#include <algorithm>
#include <cctype>
#include "common.h"
#include "bwt.h"
#include "motif.h"
#include "FFT.h"
//for directory
#include <cstdlib>

using namespace std;

//global var's declarision

extern char T[ MAX_LENGTH ];
extern int N;

extern int GenomeSize;
extern Edge* Edges;
extern Node* Nodes;
extern map<string,string> option;

void test();
int temp;

int main(int argc, char **argv)
{
    clock_t tStart,t1,t2,t2_1,t3,t4,tEnd;
    tStart=clock();
    cerr<<"START MOTIF FINDING"<<endl;
    Edges = new Edge[ HASH_TABLE_SIZE ];
    Nodes = new Node[ MAX_LENGTH * 2 ];
    genomeRegions *gR = new genomeRegions(0);
    
    //for test
    //test();

    parseCommandLine(argc, argv);
    
    Suffix active( 0, 0, -1 );  // The initial active prefix
    
    while (true) {
        
        string outPutDir(option["outdir"]); 
        system(("rm -rf "+outPutDir).c_str());
        if (system(("mkdir "+outPutDir).c_str()) != 0)
            printAndExit("cannot make directory");
        
        //parsing fasta and region
        gR->readFasta(option["fastafile"]);
        gR->readBed(option["regionfile"]);
        gR->mergeOverlap();
        gR->getSeq(outPutDir);
        //region-wide
        switch (option["bkgregion"][0]) {
            case 'r':
            case 'R':
                gR->initProb(2);
                break;
            default:
                gR->initProb(1);
                break;
        } 
        
        //initiate T 
        cerr<<"regionFile（bedFormat）:"<<option["regionfile"]<<endl;
        GenomeSize = gR->appendReverseGenome(T);
        t1=clock();
        cerr<<"Parsing Fasta:"<<double((t1-tStart)/1e6)<<endl;
        // tag mode
        if (option["mode"]=="tag") {

            
            string tagFileName(option["tagfile"]);
            //if tag file is wigfile
            if (tagFileName.substr(tagFileName.size()-3,3)=="wig") {
                cerr<<"tagfile（WigFormat）:"<<tagFileName<<endl;
                gR->readWig(tagFileName);
            }
            else {
                //problematic
                //genomeRegions tag(0);
                //cerr<<"tagfile（BedFormat）:"<<tagFileName<<endl;
                //tag.readBed(tagFileName);
                //gR->writeRawTag(tag);
            }
            //gR->getTagBed();
            t2_1=clock();
            cerr<<"Parsing Wig:"<<double((t2_1-t1)/1e6)<<endl;
            t1=t2_1;
            
            //cout<<gR->regionTags<<endl;
            for (int i=0; i<gR->tagName.size(); i++) {
                cerr<<gR->tagName[i]<<"size:"<<gR->regionTags[gR->tagName[i]].size()<<"\t";
                //cerr<<gR->regionTags[gR->tagName[i]]<<endl;
                cerr<<long(gR->regionTags[gR->tagName[i]].size())-EXTENDBOUND*4*gR->segmentCount<<endl;
                assert(GenomeSize==long(gR->regionTags[gR->tagName[i]].size())-EXTENDBOUND*4*gR->segmentCount);
            }
        }
        cerr<<"genomesize:"<<GenomeSize<<endl;
        assert(GenomeSize<=MAX_LENGTH);
        assert(GenomeSize*2<=HASH_TABLE_SIZE);
        assert(gR->segmentStartPos.back()==GenomeSize);
        N = strlen(T) - 1;
        
        for ( int i = 0 ; i <= N ; i++ )
            active.AddPrefix(i);
        t2=clock();
    
        cerr<<"built suffix tree:"<<double((t2-t2)/1e6)<<endl;
            
        //count words
        string fileName = outPutDir + "/cluster.log";
        ofstream clusterFile(fileName.c_str(),ios::app);
        //pop minimal element    heap 
        priority_queue<Motif,vector<Motif>,compareMotif> MotifHeap;
        
        const int K_4=pow1(4, K);

        int counter1 = 0;
        for (int i = 0; i <(K_4); i++) {
            printProgress(i,K_4,"Qualify kmer from suffix tree:");
            Motif thisMotif(i);
            //if (!thisMotif.noWildcard()) continue;
            //order null
            thisMotif.initProb((*gR),atoi(option["order"].c_str()));
            thisMotif.loci = active.locateMotif(thisMotif);
            //cerr<<thisMotif.loci.size()<<endl;
            temp += thisMotif.loci.size();
            thisMotif.calConscore(GenomeSize);
            if (thisMotif.score&&!thisMotif.isRepeat()) {
                if (counter1==0) {
                    MotifHeap.push(thisMotif);
                }
                counter1++;
                thisMotif.initPWM();
                thisMotif.initLociScore();
                if (option["mode"]=="tag"){
                    thisMotif.testMotifTag(*gR, outPutDir, false);
                }
                thisMotif.sumOverallScore();
                if (MotifHeap.size()<MAXMOTIFNUM||thisMotif.overallScore>MotifHeap.top().overallScore) {
                    //has pwm loci lociscore sign noise conscore motifProb overallscore.
                    MotifHeap.push(thisMotif);
                    if (MotifHeap.size()>MAXMOTIFNUM) {
                        MotifHeap.pop();
                    }
                }
                /*
                if (option["mode"]=="tag"){
                    thisMotif.testMotifTag(*gR, outPutDir, false);
                    float TagscoreThresh=MINTAGSCORE;
                    if (option["FFT"]=="T") {
                        TagscoreThresh-=MAXNOISE;
                    }
                    if (thisMotif.sumTagScore()-thisMotif.tagNoise>TagscoreThresh){
                        
                        qualifiedMotifs.push_back(thisMotif);
                    }
                    
                    else {
                        //filtered motif
                        clusterFile<<"Filtered word:\t"<<thisMotif.query<<"\t"<<thisMotif.score<<"\t"<<thisMotif.tagSTD<<"\t"<<thisMotif.tagNoise<<endl;
                        // for test
                        thisMotif.testMotifTag(*gR, outPutDir, true);
                    }
                    
                }
                else {
                    qualifiedMotifs.push_back(thisMotif);
                }
                */
            }
        }
        cerr<<"clustering result:";
                   cerr<<"\n"<<"STAGE1(filter words with low frequence):"<<K_4-counter1<<" kmer was filtered"<<"\n";
        if (option["mode"]=="tag"){
            cerr<<"STAGE2:"<<counter1-MotifHeap.size()<<" kmer was filtered"<<"\n";
        }
        cerr<<"Left:"<<MotifHeap.size()<<endl;
        cerr<<"total motifs'loci size:"<<temp<<" approximate "<<GenomeSize<<endl;
        //need to eliminate allmotif
        if (MotifHeap.size()==0) {
            printAndExit("too strict parameters!");
        }
        t3=clock();
        cerr<<"count words:"<<double((t3-t2)/1e6)<<endl;
    
        //clustering
   
        //cout<<qualifiedMotifs.size()<<qualifiedMotifs[0]<<qualifiedMotifs[1]<<endl;
        vector<Cluster> clusters;
        
        int maxMotifSize=min(MAXMOTIFNUM, int(MotifHeap.size()));
        vector<Motif> qualifiedMotifs;
        qualifiedMotifs.reserve(MAXMOTIFNUM);
        while (!MotifHeap.empty())
        {
            qualifiedMotifs.push_back(MotifHeap.top());
            MotifHeap.pop();
        }
        //qualifiedMotifs[0].printMotif();
        reverse(qualifiedMotifs.begin(), qualifiedMotifs.end());
        //qualifiedMotifs[0].printMotif();

        Cluster temp0(qualifiedMotifs[0]);
        clusters.push_back(temp0);
        
        // for (int i=0; i<10; i++) qualifiedMotifs[i].testMotifTag(*gR,outPutDir,true);
        // pair<int,int> tempa = qualifiedMotifs[2].editDistance(qualifiedMotifs[4]);
        // cout<<tempa.first<<" "<<tempa.second<<endl;
        cerr<<"clustering:"<<endl;
        for (int i=1; i<maxMotifSize; i++) { 
            float KLDiv;
            printProgress(i,maxMotifSize,"clustering:");
            int bound = clusters.size();
            for (int j=0; j<bound; j++) {
                bool aligned = false;
                pair<int,int> dist_shift = clusters[j].editDistance(qualifiedMotifs[i]);
                if (option["mode"]=="tag"){
                    KLDiv = clusters[j].tagDistrDistance(qualifiedMotifs[i]);
                    if (fabs(dist_shift.first)<=MAXDISTANCE&&KLDiv<=MAXKLDIV) {
                        aligned = true;
                    }
                }
                else {
                    if (fabs(dist_shift.first)<=MAXDISTANCE) {
                        aligned = true;
                    }
                }
                if (dist_shift.first<=0) {
                    if (aligned) {
                    //if highest mark is antisense and matched current cluster,discard
                        goto nextMotif;
                    }
                    //try next cluster
                    continue;
                }
                else if (aligned) {
                    int queryLength = clusters[j].query.size();
                    //if cluster size exceed Max cluster size after this motif appended,discard
                    if (queryLength>=atoi(option["clusterlength"].c_str())&&(dist_shift.second<0||dist_shift.second+K>queryLength)) {
                        goto nextMotif;
                    }
                    //write cluster file
                    if (option["mode"]=="tag"){
                        clusterFile<<i<<"\t";
                        if (dist_shift.second>=0) {
                            for (int counter=0; counter<dist_shift.second; counter++) {
                                clusterFile<<" ";
                            }
                        }
                        clusterFile<<" "<<qualifiedMotifs[i].query;
                        clusterFile<<qualifiedMotifs[i].tagSTD<<"dist"<<dist_shift.first<<"shift"<<dist_shift.second<<"\n";
                        clusterFile<<j<<"\t";
                        if (dist_shift.second<0) {
                            for (int counter=0; counter<abs(dist_shift.second); counter++) {
                                clusterFile<<" ";
                            }
                        }
                        clusterFile<<" "<<clusters[j].query<<clusters[j].tagSTD<<"KL div"<<KLDiv<<endl;
                    } 
                    else {
                        clusterFile<<i<<"\t";
                        if (dist_shift.second>=0) {
                            for (int counter=0; counter<dist_shift.second; counter++) {
                                clusterFile<<" ";
                            }
                        }
                        clusterFile<<" "<<qualifiedMotifs[i].query;
                        clusterFile<<"\tdist"<<dist_shift.first<<"\tshift"<<dist_shift.second<<"\n";
                        clusterFile<<j<<"\t";
                        if (dist_shift.second<0) {
                            for (int counter=0; counter<abs(dist_shift.second); counter++) {
                                clusterFile<<" ";
                            }
                        }
                        clusterFile<<" "<<clusters[j].query<<endl;
                    }
                    int prevSize = clusters[j].query.size();
                    clusters[j].concatenate(qualifiedMotifs[i],i,dist_shift.second);
                    //recalculate four attributes
                    clusters[j].calPWM(qualifiedMotifs[i], dist_shift.second);
                    clusters[j].appendLoci(qualifiedMotifs[i]);
                    clusters[j].addProb(qualifiedMotifs[i],prevSize);
                    clusters[j].calConscore(GenomeSize);
                    if (option["mode"]=="tag"){
                        clusters[j].testMotifTag(*gR,outPutDir,true);
                    }
                    clusters[j].sumOverallScore();
                    sort(clusters.begin(),clusters.end(),compareMotif());
                    goto nextMotif;
                }
                //for test: cluster system
                else if (dist_shift.first>0&&dist_shift.first<=MAXDISTANCE&&KLDiv>MAXKLDIV){
                    if (option["mode"]=="tag"){
                        clusterFile<<"KLDIV not consistent:"<<"\n"<<i<<" "<<qualifiedMotifs[i].query<<qualifiedMotifs[i].tagSTD<<"dist"<<dist_shift.first<<"shift"<<dist_shift.second<<"\n"<<j<<" "<<clusters[j].query<<clusters[j].tagSTD<<endl;
                    }
                }
                else {
                    //test next cluster
                    continue;
                }

            }
            //if not aligned to any cluster
            if (clusters.size()<MAXCLUSTERNUM){
                qualifiedMotifs[i].index = -1;
                Cluster temp0(qualifiedMotifs[i]);
                clusters.push_back(temp0);
            }
        nextMotif:
            continue;
        }
        
        //calculating
        for (int j=0; j<clusters.size(); j++){
            clusters[j].trim(); 
            clusters[j].mergeLoci();
            clusters[j].calConscore(GenomeSize);
            clusters[j].sumOverallScore();
        }
        sort(clusters.begin(),clusters.end(),compareMotif());
        for (int j=0; j<clusters.size(); j++){
            if (option["writeloci"]=="T") {
                fileName = outPutDir + "/clustersLoci.bed";
                ofstream lociFile(fileName.c_str(),ios::app);
                clusters[j].writeLoci(lociFile, *gR);
            }
        }
        t4=clock();
        cerr<<"clustering:"<<double((t4-t3)/1e6)<<endl;

        
        //write score info to cout
        cout<<"Motif\tP-value\tConscore\t";
        for (int i=0; i<gR->tagName.size(); i++) {
            cout<<gR->tagName[i]<<"_STD\t";
            
        }
        if (option["FFT"]=="T") {
            for (int i=0; i<gR->tagName.size(); i++) {
                cout<<gR->tagName[i]<<"_tagNoise\t";
                
            }
        }
        cout<<"lociSize:"<<endl;
        for (int i=0; i<clusters.size(); i++){
            clusters[i].printMotif();
        }
        
        //write pwm file;
        fileName = outPutDir + "/allmotif.pwm";
        ofstream pwmFile(fileName.c_str());
        if (!pwmFile) {
            printAndExit("Error! Fail to open pwmFile for writing!");
        }
        for (int i = 0; i<clusters.size(); i++) {
            printProgress(i,clusters.size(), "Generate PWM file");
            pwmFile<<clusters[i];
            
            
        }            
        //write pwm and dist
        
        tEnd=clock();
        cerr<<"total time eclapse:"<<double((tEnd-tStart)/1e6)<<endl;
        
        exit(0);
    }
    return 0;
}






void test(){
    //test locateSub
    //float a[64]={222196,223883,224365,225399,226522,228458,230182,232294,234638,235955,237100,238276,240537,241488,244776,247591,250682,256220,261218,266507,273588,280507,288908,296303,305834,315449,326204,338612,351108,365119,380065,392915,407475,423412,438997,444292,434680,413637,384791,362230,363998,396845,451613,503894,542050,560474,556176,536605,515048,494379,477460,460567,442164,425085,411363,396023,381178,368042,356459,344080,334199,323655,314402,307031};
    float a[64]={102.75,114.75,124.275,124.275,128,127.35,109.975,103.375,109.525,114.725,124.375,127.275,130.975,126.5,130.15,125.4,141.025,150.475,145.325,143.7,129.65,118.825,104.55,102.05,103.95,111.325,125.825,146.375,152.575,153.125,152.4,155.925,146.1,148.8,163.925,157.6,158.525,137.9,131.125,117.9,116.375,113.1,118.9,121.85,143.275,135.925,155.75,148.65,143.875,146.075,132,122.45,123.525,115.9,114.5,109.65,99.575,93.85,98.525,92.95,103.85,114.275,117.4,119.475};
    vector<float> b(a,a+64);
    //cerr<<b<<locateSubscript(b, b.begin(), b.end(), 120)<<endl;
    //test
    float total=0;
    float variance=0;
    for (int i=0; i<b.size(); i++) {
        total += b[i];
    }
    for (int i=0; i<b.size(); i++) {
        //normalization for KL divergence
        b[i] /= total;
    }
    for (int i=0; i<b.size(); i++) {
        variance += pow1(b[i]-1/float(b.size()),2);
    }

    FFT tempFFT(b);
    for (int i=0; i<b.size(); i++) {
        cerr<<tempFFT.origin[i].real()<<"\t";
    }
    cerr<<"\n"<<"tagNoise"<<tempFFT.denoise(10);
    for (int i=0; i<b.size(); i++) {
        variance += pow1(b[i]-1/float(b.size()),2);
    }
    cerr<<" std"<<sqrtf(variance)<<"\n\n";
    for (int i=0; i<b.size(); i++) {
        cerr<<norm(tempFFT.transformed[i])<<"\t";
    }
    cerr<<"\n\n";
    for (int i=0; i<b.size(); i++) {
        cerr<<tempFFT.invTrans[i].real()<<"\t";
    }
    cerr<<endl;
    assert(0==1);
}