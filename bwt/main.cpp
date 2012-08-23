
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
#include <sys/stat.h>
#include <sys/types.h>
#include "common.h"
#include "bwt.h"
#include "motif.h"
#include "FFT.h"
#include <cstdlib>

using namespace std;

//global var's declarision

extern string T;
extern int N;

extern int RegionSize;
extern Edge* Edges;
extern Node* Nodes;
extern map<string,string> option;

void test();
int temp;

int main(int argc, char **argv)
{
    clock_t tStart,t1,t2,t2_1,t3,t4,tEnd;
    tStart=clock();

    parseCommandLine(argc, argv);

    genomeRegions *gR = new genomeRegions(atoi(option["extend"].c_str()));

    cerr<<"START MOTIF FINDING"<<endl;
    
    while (true) {
        
        string outPutDir(option["outdir"]); 
        system(("rm -rf "+outPutDir).c_str());
        
        if (mkdir(outPutDir.c_str(),S_IRWXU|S_IRWXG|S_IRWXO)!= 0)
            printAndExit("cannot make directory");
        
        //parsing fasta and region
        if (option["regionfastafile"]=="") {
            gR->readFasta(option["fastafile"]);
            gR->readBed(option["regionfile"]);
            gR->mergeOverlap();
            // rm control peaks
            if (option["control"]!="") {
                cerr<<"control File（bedFormat）:"<<option["control"]<<endl;
                gR->rmControlPeaks(option["control"]);
            }
            gR->getSeq(outPutDir);
        }
        else {
            gR->readRegionFasta(option["regionfastafile"]);
        }
        
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
        
        //initiate whole sequence 
        cerr<<"regionFile（bedFormat）:"<<option["regionfile"]<<endl;
        RegionSize = gR->appendReverseGenome(T);
        t1=clock();
        cerr<<"Parsing Fasta:"<<double((t1-tStart)/1e6)<<endl;

    
        // tag mode
        if (option["mode"]=="tag") {

            string fileName("");
            string tagFileNames(option["tagfile"]);
            istringstream ss(tagFileNames);
            while (ss>>fileName) {
                cerr<<"tagfile:"<<fileName<<endl;
                gR->readWig(fileName);
            }
            gR->catenateTags();
            t2_1=clock();
            cerr<<"Parsing Wig:"<<double((t2_1-t1)/1e6)<<endl;
            t1=t2_1;
            
            for (int i=0; i<gR->tagName.size(); i++) {
                cerr<<gR->tagName[i]<<" size:"<<gR->regionTags[gR->tagName[i]].size()<<"\t";
                cerr<<long(gR->regionTags[gR->tagName[i]].size())-EXTENDBOUND*4*gR->segmentCount<<endl;
                assert(RegionSize==long(gR->regionTags[gR->tagName[i]].size())-EXTENDBOUND*4*gR->segmentCount);
            }
        }
        cerr<<"RegionSize:"<<RegionSize<<endl;
        assert(gR->segmentStartPos.back()==RegionSize);
        
        Edges = new Edge[long(RegionSize*2*HASH_TABLE_MULTIPLER)];
        Nodes = new Node[RegionSize*2+1];
        HASH_TABLE_SIZE =long(RegionSize*2*HASH_TABLE_MULTIPLER);
        
        N = T.size() - 1;
        
        Suffix active( 0, 0, -1 );  // The initial active prefix
        for ( int i = 0 ; i <= N ; i++ )
            active.AddPrefix(i);
        t2=clock();
    
        cerr<<"built suffix tree:"<<double((t2-t2)/1e6)<<endl;
            
        //count words
        string fileName = outPutDir + "/cluster.log";
        ofstream clusterFile(fileName.c_str(),ios::app);
        
        //pop minimal element from heap 
        priority_queue<Motif,vector<Motif>,compareMotif> MotifHeap;
        
        const int K_4=int(pow1(4, K));
        int counter1 = 0;
        bool rmrepeat = (option["rmrepeat"]!="0");
        for (int i = 0; i <(K_4); i++) {
            printProgress(i,K_4,"Qualify kmer from suffix tree:");
            Motif thisMotif(i);
            //if (!thisMotif.noWildcard()) continue;
            //order null
            thisMotif.initProb((*gR),atoi(option["order"].c_str()));
            thisMotif.loci = active.locateMotif(thisMotif);
            //cerr<<thisMotif.loci.size()<<endl;
            temp += thisMotif.loci.size();
            thisMotif.calConscore(RegionSize);
            if (thisMotif.score) {
                if (rmrepeat&&thisMotif.isRepeat()) {
                    continue;
                }
                counter1++;
                thisMotif.initPWM();
                thisMotif.initLociScore();
                if (option["mode"]=="tag"){
                    thisMotif.testMotifTag(*gR, false);
                }
                thisMotif.sumOverallScore();
                if (counter1==1) {
                    MotifHeap.push(thisMotif);
                }
                if (MotifHeap.size()<MAXMOTIFNUM||thisMotif.overallScore>MotifHeap.top().overallScore) {
                    //has pwm loci lociscore sign noise conscore motifProb overallscore.
                    MotifHeap.push(thisMotif);
                    if (MotifHeap.size()>MAXMOTIFNUM) {
                        MotifHeap.pop();
                    }
                }
            }
        }
        cerr<<"clustering result:";
        cerr<<"\n"<<"STAGE1(filter words with low frequence):"<<K_4-counter1<<" kmer was filtered"<<"\n";
        if (option["mode"]=="tag"){
            cerr<<"STAGE2:"<<counter1-MotifHeap.size()<<" kmer was filtered"<<"\n";
        }
        cerr<<"Left:"<<MotifHeap.size()<<endl;
        cerr<<"total motifs'loci size:"<<temp<<" approximate "<<RegionSize<<endl;
        //need to eliminate allmotif
        if (MotifHeap.size()==0) {
            printAndExit("too strict parameters!");
        }
        t3=clock();
        
        cerr<<"count words:"<<double((t3-t2)/1e6)<<endl;
        
        delete [] Nodes;
        delete [] Edges;
        
    
        //clustering
   
        vector<Cluster> clusters;
        int exceedLengthCount = 0;
        int inverseAlignedCount = 0;
        int normalAligned = 0;
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
#ifdef QUALIFIED
        ofstream wordFile((option["outdir"]+"/qualified.log").c_str());
        ofstream wordDist((option["outdir"]+"/qualified.dist").c_str());
        for (int i=0; i<qualifiedMotifs.size(); i++) {
            qualifiedMotifs[i].printMotif(wordFile);
        }
        for (int i=0; i<qualifiedMotifs.size(); i++) {
            if (qualifiedMotifs[i].implicit()) {
                continue;
            }
            qualifiedMotifs[i].drawDist(*gR, wordDist);
        } 
#endif
#ifdef CHIPEDPEAKDIST
        //output dist centered by chiped peaks
        ofstream CHIPDist((option["outdir"]+"/ChIPed_Peaks.dist").c_str());
        Motif chipPeaks(1);
        chipPeaks.loci.clear();
        for (int i=0; i<gR->segmentStartPos.size()-1 ;i++) {
            chipPeaks.loci.push_back(gR->segmentStartPos[i]+50);
        }
        chipPeaks.testMotifTag(*gR, false);
        chipPeaks.drawDist(*gR, CHIPDist);
#endif
        
        Cluster temp0(qualifiedMotifs[0]);
        
        clusters.push_back(temp0);
        
        for (int i=1; i<maxMotifSize; i++) { 
            float KLDiv=0;
            printProgress(i-1,maxMotifSize-1,"clustering:");
            int bound = clusters.size();
            for (int j=0; j<bound; j++) {
                bool aligned = false;
                pair<float,int> dist_shift = clusters[j].editDistance(qualifiedMotifs[i]);
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
                        inverseAlignedCount++;
                        goto nextMotif;
                    }
                    //try next cluster
                    continue;
                }
                else if (aligned) {
                    int queryLength = clusters[j].query.size();
                    //if cluster size exceed Max cluster size after this motif appended,discard
                    if (queryLength>=atoi(option["clusterlength"].c_str())&&(dist_shift.second<0||dist_shift.second+K>queryLength)) {
                        exceedLengthCount++;
                        goto nextMotif;
                    }
               
#ifdef CLUSTERLOG
                    //write cluster file
                    if (option["mode"]=="tag"){
                        clusterFile<<i<<setw(8)<<"\t";
                        if (dist_shift.second>=0) {
                            for (int counter=0; counter<dist_shift.second; counter++) {
                                clusterFile<<" ";
                            }
                        }
                        clusterFile<<" "<<qualifiedMotifs[i].query;
                        clusterFile<<qualifiedMotifs[i].tagBiPeak<<qualifiedMotifs[i].tagSymmetry<<clusters[j].signalIntensity<<"dist"<<dist_shift.first<<"shift"<<dist_shift.second<<"\n";
                        clusterFile<<j<<setw(8)<<"\t";
                        if (dist_shift.second<0) {
                            for (int counter=0; counter<abs(dist_shift.second); counter++) {
                                clusterFile<<" ";
                            }
                        }
                        clusterFile<<" "<<clusters[j].query<<clusters[j].tagBiPeak<<clusters[j].tagSymmetry<<clusters[j].signalIntensity<<"KL div"<<KLDiv<<endl;
                    } 
                    else {
                        clusterFile<<i<<setw(8)<<"\t";
                        if (dist_shift.second>=0) {
                            for (int counter=0; counter<dist_shift.second; counter++) {
                                clusterFile<<" ";
                            }
                        }
                        clusterFile<<" "<<qualifiedMotifs[i].query;
                        clusterFile<<"\tdist"<<dist_shift.first<<"\tshift"<<dist_shift.second<<"\n";
                        clusterFile<<j<<setw(8)<<"\t";
                        if (dist_shift.second<0) {
                            for (int counter=0; counter<abs(dist_shift.second); counter++) {
                                clusterFile<<" ";
                            }
                        }
                        clusterFile<<" "<<clusters[j].query<<endl;
                    }
#endif               
                    normalAligned++;
                    int prevSize = clusters[j].query.size();
                    clusters[j].concatenate(qualifiedMotifs[i],i,dist_shift.second);
                    //recalculate four attributes
                    clusters[j].calPWM(qualifiedMotifs[i], dist_shift.second);
                    if (option["mode"]=="tag"){
                        clusters[j].reCalSumBin(qualifiedMotifs[i], *gR);
                    }
                    clusters[j].appendLoci(qualifiedMotifs[i]);
                    clusters[j].addProb(qualifiedMotifs[i],prevSize);
                    clusters[j].calConscore(RegionSize);
                    if (option["mode"]=="tag"){
                        //update noise and tagscore
                        clusters[j].testMotifTag(*gR,false);
                    }
                    clusters[j].sumOverallScore();
                    sort(clusters.begin(),clusters.end(),compareMotif());
                    goto nextMotif;
                }
                //for test: cluster system not aligned by tag
#ifdef CLUSTERLOG   
                else if (dist_shift.first>0&&dist_shift.first<=MAXDISTANCE&&KLDiv>MAXKLDIV){
                    if (option["mode"]=="tag"){
                        clusterFile<<"KLDIV not consistent:"<<"\n";
                        clusterFile<<i<<setw(8)<<"\t";
                        if (dist_shift.second>=0) {
                            for (int counter=0; counter<dist_shift.second; counter++) {
                                clusterFile<<" ";
                            }
                        }
                        clusterFile<<" "<<qualifiedMotifs[i].query;
                        clusterFile<<qualifiedMotifs[i].tagBiPeak<<qualifiedMotifs[i].tagSymmetry<<"dist"<<dist_shift.first<<"shift"<<dist_shift.second<<"\n";
                        clusterFile<<j<<setw(8)<<"\t";
                        if (dist_shift.second<0) {
                            for (int counter=0; counter<abs(dist_shift.second); counter++) {
                                clusterFile<<" ";
                            }
                        }
                        clusterFile<<" "<<clusters[j].query<<clusters[j].tagBiPeak<<clusters[j].tagSymmetry<<"KL div"<<KLDiv<<endl;
                    }
                }
#endif
                else {
                    //test next cluster
                    continue;
                }

            }
            //if not aligned to any cluster, new cluster
            if (clusters.size()<MAXCLUSTERNUM){
                qualifiedMotifs[i].index = -1;
                clusterFile<<"+new cluster: "<<i<<" "<<qualifiedMotifs[i].query<<endl;
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
            clusters[j].testMotifTag(*gR,true);
            clusters[j].calConscore(RegionSize);
            clusters[j].sumOverallScore();
        }
        sort(clusters.begin(),clusters.end(),compareMotif());
        t4=clock();
        cerr<<"clustering:"<<double((t4-t3)/1e6)<<endl;
        cerr<<"skipped:"<<exceedLengthCount<<"(exceedLength)+"<<inverseAlignedCount<<"(inverseAligned)="<<exceedLengthCount+inverseAlignedCount<<"\tAligned:"<<normalAligned<<"\tClusters:"<<clusters.size()<<"\ttotal:"<<exceedLengthCount+normalAligned+inverseAlignedCount+clusters.size()<<endl;
        if (option["writeloci"]=="T") {
            clock_t t5=clock();
            for (int j=0; j<clusters.size(); j++){
                fileName = outPutDir + "/clustersLoci.bed";
                ofstream lociFile(fileName.c_str(),ios::app);
                clusters[j].writeLoci(lociFile, *gR);
            }
            cerr<<"writeLoci:"<<double((t5-t4)/1e6)<<endl;
            t4=t5;
        }
        
        
        
        //write score info to cout
        cout<<setw(20)<<setiosflags(std::ios::left)<<"Motif"<<"\tP-value\tConscore\t";
        for (int i=0; i<gR->tagName.size(); i++) {
            cout<<gR->tagName[i]<<"_Bipeaks\t";
            
        }
        for (int i=0; i<gR->tagName.size(); i++) {
            cout<<gR->tagName[i]<<"_assymmetry\t";
            
        }
        if (option["FFT"]=="T") {
            for (int i=0; i<gR->tagName.size(); i++) {
                cout<<gR->tagName[i]<<"_tagNoise\t";
                
            }
        }
        cout<<"lociSize:"<<endl;
        for (int i=0; i<clusters.size(); i++){
            clusters[i].printMotif(cout);
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




