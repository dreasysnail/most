//
//  main.cpp
//  bwt
//
//  Created by icarus on 11-10-4.
//  Copyright 2011年 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include <string>
#include <algorithm>
#include <cctype>
#include "bwt.h"
#include "common.h"
#include "motif.h"
//for directory
#include <cstdlib>

using namespace std;

//global var's declarision

extern char T[ MAX_LENGTH ];
extern int N;
void printUsage();
extern int GenomeSize;
extern Edge* Edges;
extern Node* Nodes;
void preInitialization(Edge* Edges,Node* Nodes);
void test();
int temp;

int main(int argc, char **argv)
{
    clock_t tStart,t1,t2,t2_1,t3,t4,tEnd;
    tStart=clock();
    cerr<<"START MOTIF FINDING"<<endl;
    Edge *EdgesTemp = new Edge[ HASH_TABLE_SIZE ];
    Node *NodesTemp = new Node[ MAX_LENGTH * 2 ];
    genomeRegions *gR = new genomeRegions(0);
    
    //for test
    //test();
    map<string,string> option;
    parseCommandLine(argc, argv,option);
    
        
    preInitialization(EdgesTemp,NodesTemp);
    
    
    Suffix active( 0, 0, -1 );  // The initial active prefix
    
    while (true) {

        // tag mode
        if (option["mode"]=="tag") {
            /************************
            2:tss.bed
            3:fasta
            4:histone.wig/histone.bed
            6:output directory
            ************************/
            string outPutDir(option["outdir"]); 
            system(("rm -rf "+outPutDir).c_str());
            if (system(("mkdir "+outPutDir).c_str()) != 0)
            {
                cerr<<"cannot make directory"<<endl;
                exit(1);
            }
            
            

            gR->readFasta(option["fastafile"]);
            gR->readBed(option["regionfile"]);
            gR->mergeOverlap();
            gR->getSeq(outPutDir);
            //region-wide
            //gR->initProb(2);
            //initiate T 
            cerr<<"regionFile（bedFormat）:"<<option["regionfile"]<<endl;
            GenomeSize = gR->catSeq(T);
            t1=clock();
            cerr<<"Parsing Fasta:"<<double((t1-tStart)/1e6)<<endl;
            
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
            t2=clock();
            cerr<<"Parsing Wig:"<<double((t2-t1)/1e6)<<endl;


            
            
            //ofstream tfile("aa.fa");
            //tfile<<T<<endl;
            //cerr<<T<<endl;
            //cout<<gR->regionTags<<endl;
            for (int i=0; i<gR->tagName.size(); i++) {
                cerr<<gR->tagName[i]<<"size:"<<gR->regionTags[gR->tagName[i]].size()<<"\t";
                //cerr<<gR->regionTags[gR->tagName[i]]<<endl;
                cerr<<long(gR->regionTags[gR->tagName[i]].size())-EXTENDBOUND*4*gR->segmentCount<<endl;
                assert(GenomeSize==long(gR->regionTags[gR->tagName[i]].size())-EXTENDBOUND*4*gR->segmentCount);
            }
            cerr<<"genomesize:"<<GenomeSize<<endl;
            assert(GenomeSize<=MAX_LENGTH);
            assert(GenomeSize*2<=HASH_TABLE_SIZE);
            assert(gR->segmentStartPos.back()==GenomeSize);
            N = strlen(T) - 1;
            
            for ( int i = 0 ; i <= N ; i++ )
                active.AddPrefix(i);
            t2_1=clock();
        
            cerr<<"built suffix tree:"<<double((t2_1-t2)/1e6)<<endl;
            
            //count words
            string fileName = outPutDir + "/cluster_signif";
            ofstream clusterFile(fileName.c_str(),ios::app);
            vector<Motif> qualifiedMotifs;
            const int K_4=pow1(4, K);
            qualifiedMotifs.reserve(K_4/1000);
            int counter1 = 0;
            for (int i = 0; i <(K_4); i++) {
                printProgress(i,K_4,"Qualify kmer from suffix tree:");
                Motif thisMotif(i);
                //if (!thisMotif.noWildcard()) continue;
                //order null
                thisMotif.initProb((*gR),-1);
                thisMotif.loci = active.locateMotif(thisMotif);
                
                //cout<<allmotifs[i].loci.size()<<" "<<active.countString(allmotifs[i].query)<<endl;
                temp += thisMotif.loci.size();
                thisMotif.calConscore(GenomeSize);
                if (thisMotif.score&&!thisMotif.isRepeat()) {
                    counter1++;
                    thisMotif.initPWM();
                    thisMotif.testMotifTag(*gR, outPutDir, false);
                    if (thisMotif.sumTagScore()>MINTAGSCORE){
                        //has pwm loci sign conscore motifProb.
                        qualifiedMotifs.push_back(thisMotif);
                    }
                    else {
                        clusterFile<<"filtered:\t"<<thisMotif.query<<"\t"<<thisMotif.score<<"\t"<<thisMotif.signif<<endl;
                    }
                }
            }
            cerr<<"clustering result:";
            cerr<<"\n"<<"STAGE1:"<<K_4-counter1<<" kmer was filtered"<<"\n";
            cerr<<"STAGE2:"<<counter1-qualifiedMotifs.size()<<" kmer was filtered"<<"\n";
            cerr<<"Left:"<<qualifiedMotifs.size()<<endl;
            cerr<<"total motifs'loci size:"<<temp<<" approximate "<<GenomeSize<<endl;
            //need to eliminate allmotif
            if (qualifiedMotifs.size()==0) {
                printAndExit("too strict parameters!");
            }
            t3=clock();
            cerr<<"count words:"<<double((t3-t2_1)/1e6)<<endl;
            //clustering
            sort(qualifiedMotifs.begin(), qualifiedMotifs.end(),compareMotif());
       
            //cout<<qualifiedMotifs.size()<<qualifiedMotifs[0]<<qualifiedMotifs[1]<<endl;
            vector<Cluster> clusters;
            int maxMotifSize=min(MOTIFMAX, int(qualifiedMotifs.size()));
            
            Cluster temp0(qualifiedMotifs[0]);
            clusters.push_back(temp0);
            
            for (int i=0; i<10; i++) qualifiedMotifs[i].testMotifTag(*gR,outPutDir,true);
            // pair<int,int> tempa = qualifiedMotifs[2].editBINSPANance(qualifiedMotifs[4]);
            // cout<<tempa.first<<" "<<tempa.second<<endl;
            cerr<<"clustering:"<<endl;
            for (int i=1; i<maxMotifSize; i++) { 
                printProgress(i,maxMotifSize,"clustering:");
                int bound = clusters.size();
                for (int j=0; j<bound; j++) {
                    pair<int,int> dist_shift = clusters[j].editDistance(qualifiedMotifs[i]);
                    float signifDist = clusters[j].histoneDistrDistance(qualifiedMotifs[i]);
                    //if highest mark is antisense and matched current cluster,discard
                    if (dist_shift.first<=0&&(-dist_shift.first)<=MAXDISTANCE) {
                        goto nextMotif;
                    }
                    //try next cluster
                    else if (dist_shift.first<=0&&(-dist_shift.first)>MAXDISTANCE){
                        continue;
                    }
                    //else if (dist_shift.first<=MAXDISTANCE&&signifDist<=MAXKLDIV) {
                    else if (dist_shift.first<=MAXDISTANCE) {
                        int queryLength = clusters[j].query.size();
                        //if cluster size exceed Max cluster size after this motif appended,discard
                        if (queryLength>=MAXCLUSTERSIZE&&(dist_shift.second<0||dist_shift.second+K>queryLength)) {
                            goto nextMotif;
                        }
                        clusterFile<<i<<" "<<qualifiedMotifs[i].query<<qualifiedMotifs[i].signif<<"dist"<<dist_shift.first<<"shift"<<dist_shift.second<<"\n"<<j<<"KL div"<<signifDist<<clusters[j].query<<clusters[j].signif<<endl;
                        int prevSize = clusters[j].query.size();
                        clusters[j].concatenate(qualifiedMotifs[i],i,dist_shift.second);
                        //recalculate four attributes
                        clusters[j].calPWM(qualifiedMotifs[i], dist_shift.second);
                        clusters[j].appendLoci(qualifiedMotifs[i]);
                        clusters[j].testMotifTag(*gR,outPutDir,true);
                        clusters[j].addProb(qualifiedMotifs[i],prevSize);
                        clusters[j].calConscore(GenomeSize);
                        clusters[j].sumOverallScore();
                        sort(clusters.begin(),clusters.end(),compareCluster());
                        goto nextMotif;
                    }
                    //for test: cluster system 
                    else if (dist_shift.first>0&&dist_shift.first<=MAXDISTANCE&&signifDist>MAXKLDIV){
                        clusterFile<<"signif not consistent:"<<"\n"<<i<<" "<<qualifiedMotifs[i].query<<qualifiedMotifs[i].signif<<"dist"<<dist_shift.first<<"shift"<<dist_shift.second<<"\n"<<j<<" "<<clusters[j].query<<clusters[j].signif<<endl;
                    }
                    else {
                        //test next cluster
                        continue;
                    }

                }
                //if not aligned to any cluster
                if (clusters.size()<CLUSTERMAX){
                    if (qualifiedMotifs[i].isRepeat()) {
                        continue;
                    }
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
                //cout<<clusters[i].loci.size()<<endl;
                //cout<<clusters[i].loci.size()<<endl;
                //calculate score for each cluster
            }
            t4=clock();
            cerr<<"clustering:"<<double((t4-t3)/1e6)<<endl;

            
            //write score info to cout
            cout<<"Motif\tCONscore\t";
            for (int i=0; i<gR->tagName.size(); i++) {
                cout<<gR->tagName[i]<<"\t";
                
            }
            cout<<"lociSize:"<<endl;
            for (int i=0; i<clusters.size(); i++){
                clusters[i].printMotif();
            }
            
            //write pwm file;
            fileName = outPutDir + "/allmotif.pwm";
            ofstream pwmFile(fileName.c_str());
            if (!pwmFile) {
                string errorInfo = "Error! Fail to open pwmFile for writing!";
                printAndExit(errorInfo);
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
        else if (option["mode"]=="normal") {
    
            string outPutDir(option["outdir"]); 
            system(("rm -rf "+outPutDir).c_str());
            if (system(("mkdir "+outPutDir).c_str()) != 0)
            {
                cerr<<"cannot make directory"<<endl;
                exit(1);
            }
            
            
            
            gR->readFasta(option["fastafile"]);
            gR->readBed(option["regionfile"]);
            gR->mergeOverlap();
            gR->getSeq(outPutDir);
            //region-wide
            //gR->initProb(2);
            //initiate T 
            cerr<<"regionFile（bedFormat）:"<<option["regionfile"]<<endl;
            GenomeSize = gR->catSeq(T);
            t1=clock();
            cerr<<"Parsing Fasta:"<<double((t1-tStart)/1e6)<<endl;
            //ofstream tfile("aa.fa");
            //tfile<<T<<endl;
            //cerr<<T<<endl;
            cerr<<"genomesize:"<<GenomeSize<<endl;
            assert(GenomeSize<=MAX_LENGTH);
            assert(GenomeSize*2<=HASH_TABLE_SIZE);
            N = strlen(T) - 1;
            
            for ( int i = 0 ; i <= N ; i++ )
                active.AddPrefix(i);
            t2=clock();
            cerr<<"built suffix tree:"<<double((t2-t1)/1e6)<<endl;
            
            //count words
            vector<Motif> qualifiedMotifs;
            const int K_4=pow1(4, K);
            qualifiedMotifs.reserve(K_4/1000);
            int counter1 = 0;
            for (int i = 0; i <(K_4); i++) {
                printProgress(i,K_4,"Qualify kmer from suffix tree:");
                Motif thisMotif(i);
                //if (!thisMotif.noWildcard()) continue;
                //order null
                thisMotif.initProb((*gR),atoi(option["order"].c_str()));
                thisMotif.loci = active.locateMotif(thisMotif);
                
                //cout<<allmotifs[i].loci.size()<<" "<<active.countString(allmotifs[i].query)<<endl;
                temp += thisMotif.loci.size();
                thisMotif.calConscore(GenomeSize);
                if (thisMotif.score&&!thisMotif.isRepeat()) {
                    counter1++;
                    thisMotif.initPWM();
                    qualifiedMotifs.push_back(thisMotif);
                }
            }
            cerr<<"\n"<<"STAGE1(filter words with low frequence):"<<K_4-counter1<<" kmer was filtered"<<"\n";
            cerr<<"Left:"<<qualifiedMotifs.size()<<endl;
            cerr<<"total motifs'loci size:"<<temp<<" approximate "<<GenomeSize<<endl;
            //need to eliminate allmotif
            if (qualifiedMotifs.size()==0) {
                cerr<<"too strict parameters!";
                exit(1);
            }
            t3=clock();
            cerr<<"count words:"<<double((t3-t2)/1e6)<<endl;
            //clustering
            sort(qualifiedMotifs.begin(), qualifiedMotifs.end(),compareMotif());
            vector<Cluster> clusters;
            int maxMotifSize=min(MOTIFMAX, int(qualifiedMotifs.size()));
            Cluster temp0(qualifiedMotifs[0]);
            clusters.push_back(temp0);
            // pair<int,int> tempa = qualifiedMotifs[2].editDistance(qualifiedMotifs[4]);
            // cout<<tempa.first<<" "<<tempa.second<<endl;
            cerr<<"clustering:"<<endl;
            for (int i=1; i<maxMotifSize; i++) { 
                printProgress(i,maxMotifSize,"clustering:");
                int bound = clusters.size();
                for (int j=0; j<bound; j++) {
                    pair<int,int> dist_shift = clusters[j].editDistance(qualifiedMotifs[i]);
                    //if highest mark is antisense and matched current cluster,discard
                    if (dist_shift.first<=0) {
                        goto nextMotif2;
                    }
                    //try next cluster
                    else if (dist_shift.first<=0&&(-dist_shift.first)>MAXDISTANCE){
                        continue;
                    }
                    //else if (dist_shift.first<=MAXDISTANCE&&signifDist<=MAXKLDIV) {
                    else if (dist_shift.first<=MAXDISTANCE) {
                        int queryLength = clusters[j].query.size();
                        //if cluster size exceed Max cluster size after this motif appended,discard
                        if (queryLength>=MAXCLUSTERSIZE&&(dist_shift.second<0||dist_shift.second+K>queryLength)) {
                            goto nextMotif2;
                        }
                        int prevSize = clusters[j].query.size();
                        clusters[j].concatenate(qualifiedMotifs[i],i,dist_shift.second);
                        //recalculate four attributes
                        clusters[j].calPWM(qualifiedMotifs[i], dist_shift.second);
                        clusters[j].appendLoci(qualifiedMotifs[i]);
                        clusters[j].addProb(qualifiedMotifs[i],prevSize);
                        clusters[j].calConscore(GenomeSize);
                        sort(clusters.begin(),clusters.end(),compareMotif());
                        goto nextMotif2;
                    }
                    //for test: cluster system 
                    else {
                        //test next cluster
                        continue;
                    }
                    
                }
                //if not aligned to any cluster
                if (clusters.size()<CLUSTERMAX){
                    if (qualifiedMotifs[i].isRepeat()) {
                        continue;
                    }
                    qualifiedMotifs[i].index = -1;
                    Cluster temp0(qualifiedMotifs[i]);
                    clusters.push_back(temp0);
                }
            nextMotif2:
                continue;
            }
            //calculating
            for (int j=0; j<clusters.size(); j++){
                clusters[j].trim();            
                //cout<<clusters[i].loci.size()<<endl;
                //cout<<clusters[i].loci.size()<<endl;
                //calculate score for each cluster
            }
            t4=clock();
            cerr<<"clustering:"<<double((t4-t3)/1e6)<<endl;
            
            
            //write score info to cout
            cout<<"Motif\tCONscore\t";
            cout<<"lociSize:"<<endl;
            for (int i=0; i<clusters.size(); i++){
                clusters[i].printMotif();
            }
            
            //write pwm file;
            string fileName = outPutDir + "/allmotif.pwm";
            ofstream pwmFile(fileName.c_str());
            if (!pwmFile) {
                string errorInfo = "Error! Fail to open pwmFile for writing!";
                printAndExit(errorInfo);
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
        else if (option["mode"]=="help") {
            printUsage();
            exit(0);
        }
        printUsage();
        exit(1);     
    }
};

void printUsage()
{
	string usage =	"----------------------------------------------------------------------\n";
	usage		+=	"        Motif discovery system \n";
	usage		+=	"        Developed by \n";
	usage		+=	"        jeremy071242044@gmail.com\n";
	usage		+=	"----------------------------------------------------------------------\n";
	usage		+=	"    motif [option] <parameter1> <parameter2> \n";
	usage		+=	"  Options:\n";
	usage		+=	"    -m --manually input\t<DNA sequence file>\t<Annotation file>\n\t\tTrain CTF on dataset given in parameters\n\n";
	usage		+=	"    -i --file\t<bed file>\t<DNA sequence file>\n\t\tPredict on DNA sequences\n\n";
    usage		+=	"    -t --file\t<bed file>\t<DNA sequence file>\n\t<tag bed file>\tPredict on DNA sequences\n\n";
	usage		+=	"    -h --help\n\t\tPrint help information.\n";
    usage		+=	"    -d --debug\n\t\trun some test.\n";
    usage		+=	"    ATTENTION:LAST FILE SHOULD BE SORTED";
	cerr<<usage<<endl;
    
    cerr << "Normally, suffix trees require that the last\n"
    << "character in the input string be unique.  If\n"
    << "you don't do this, your tree will contain\n"
    << "suffixes that don't end in leaf nodes.  This is\n"
    << "often a useful requirement. You can build a tree\n"
    << "in this program without meeting this requirement,\n"
    << "but the validation code will flag it as being an\n"
    << "invalid tree\n\n"<<endl;
}


void preInitialization(Edge* EdgesTemp,Node* NodesTemp){
    Nodes = NodesTemp;
    Edges = EdgesTemp;

}


void test(){
    //test locateSub
    int a[10]={0,3,10,14,20,30,40,60,90,120};
    vector<int> b(a,a+10);
    cerr<<b<<locateSubscript(b, b.begin(), b.end(), 120)<<endl;
    //test
    
    
    
    
    assert(0==1);
}
