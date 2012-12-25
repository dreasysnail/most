
//
//  most.cpp
//  MOST
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
    
    /********************* parse files ***********************/

    parseCommandLine(argc, argv);
    genomeRegions *gR = new genomeRegions(atoi(option["flanking"].c_str()));
    
    //for test
    //test();

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
            if (option["mode"]=="extract")
                exit(0);
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
                cerr<<long(gR->regionTags[gR->tagName[i]].size())-EXTEND_BOUND*4*gR->segmentCount<<endl;
                assert(RegionSize==long(gR->regionTags[gR->tagName[i]].size())-EXTEND_BOUND*4*gR->segmentCount);
            }
        }
        cerr<<"RegionSize:"<<RegionSize<<endl;
        assert(gR->segmentStartPos.back()==RegionSize);
        
        
        /********************* Count occurence ***********************/
        
        
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
                thisMotif.initPFM();
                //thisMotif.initLociScore();
                if (option["mode"]=="tag"){
                    thisMotif.testMotifTag(*gR, false);
                }
                thisMotif.sumOverallScore();
                if (counter1==1) {
                    MotifHeap.push(thisMotif);
                }
                if (MotifHeap.size()<MAX_WORD_NUM||thisMotif.overallScore>MotifHeap.top().overallScore) {
                   //protocol:has pfm loci lociscore sign noise conscore motifProb overallscore  bins  sumbin.
                    MotifHeap.push(thisMotif);
                    if (MotifHeap.size()>MAX_WORD_NUM) {
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
        

        
    
        /********************* Clustering ***********************/
        
        vector<Cluster> clusters;
        int exceedLengthCount = 0;
        int inverseAlignedCount = 0;
        int normalAligned = 0;
        int maxMotifSize=min(MAX_WORD_NUM, int(MotifHeap.size()));
        vector<Cluster> qualifiedMotifs;
        qualifiedMotifs.reserve(MAX_WORD_NUM);
        
        
        //extend motif
        if (option["extendmotif"]=="T") {
            while (!MotifHeap.empty())
            {
                Cluster temp0(Motif (1));
                temp0.getExtended(MotifHeap.top(), *gR, active);
                qualifiedMotifs.push_back(temp0);
                MotifHeap.pop();
            }
        }
        else {
            while (!MotifHeap.empty())
            {
                Cluster temp0(MotifHeap.top());
                qualifiedMotifs.push_back(temp0);
                MotifHeap.pop();
            }
        }
        
        
        
        //??
        delete [] Nodes;
        delete [] Edges;
        
        //qualifiedMotifs[0].printMotif();
        sort(qualifiedMotifs.begin(),qualifiedMotifs.end(),compareMotif());
        
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
#ifdef QUALIFIEDPFM
        //write pfm file;
        ofstream pfmWordFile((option["outdir"]+"/qualified.pfm").c_str());
        if (!pfmWordFile) {
            printAndExit("Error! Fail to open pfmFile for writing!");
        }
        for (int i = 0; i<qualifiedMotifs.size(); i++) {
            printProgress(i,qualifiedMotifs.size(), "Generate PFM file for qualified words");
            pfmWordFile<<qualifiedMotifs[i];
        }
#endif
#ifdef CHIPEDPEAKDIST
        //output dist centered by chiped peaks
        ofstream CHIPDist((option["outdir"]+"/ChIPed_Peaks.dist").c_str());
        Motif chipPeaks(1);
        chipPeaks.loci.clear();
        for (int i=0; i<gR->segmentStartPos.size()-2 ;i++) {
            chipPeaks.loci.push_back((gR->segmentStartPos[i]+gR->segmentStartPos[i+1])/2);
        }
        for (int j=0; j<clusters.size(); j++){
            fileName = outPutDir + "/chippeaks.bed";
            ofstream lociFile(fileName.c_str(),ios::app);
            chipPeaks.writeLoci(lociFile, *gR);
        }

        chipPeaks.testMotifTag(*gR, false);
        chipPeaks.drawDist(*gR, CHIPDist);
#endif
        //qualifiedMotifs[0].printMotif();
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
                    if (fabs(dist_shift.first)<=MAX_DIST&&KLDiv<=MAX_KL_DIV) {
                        aligned = true;
                    }
                }
                else {
                    if (fabs(dist_shift.first)<=MAX_DIST) {
                        aligned = true;
                    }
                }
                //reverse aligned
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
                    int queryLength = clusters[j].pfm[0].size();
                    //if cluster size exceed Max cluster size after this motif appended,discard
                    if (queryLength>=atoi(option["clusterlength"].c_str())&&(dist_shift.second<0||dist_shift.second+K>queryLength)) {
                        exceedLengthCount++;
                        goto nextMotif;
                    }
                    Member tempM ={dist_shift.second,qualifiedMotifs[i].tagBiPeak,qualifiedMotifs[i].tagSymmetry,dist_shift.first,qualifiedMotifs[i].query};
                    clusters[j].Members.push_back(tempM);
              
                    normalAligned++;
                   
                    //recalculate four attributes
                    clusters[j].calPFM(qualifiedMotifs[i], dist_shift.second);
                    if (option["mode"]=="tag"){
                        clusters[j].reCalSumBin(qualifiedMotifs[i], *gR);
                    }
                    clusters[j].appendLoci(qualifiedMotifs[i]);
                    //merge score and prob
                    clusters[j].mergeProb(qualifiedMotifs[i]);
                    //update noise and tagscore
                    clusters[j].testMotifTag(*gR,false);
                    clusters[j].sumOverallScore();
                    sort(clusters.begin(),clusters.end(),compareMotif());
                    goto nextMotif;
                }
                //for test: cluster system not aligned by tag
#ifdef CLUSTERLOG   
                else if (dist_shift.first>0&&dist_shift.first<=MAX_DIST&&KLDiv>MAX_KL_DIV){
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
            if (clusters.size()<MAX_CLUSTER_NUM){
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
            
            clusters[j].generateIUPAC();
            clusters[j].trim();
            clusters[j].mergeLoci();
            clusters[j].testMotifTag(*gR,option["drawdist"]=="T");
            //clusters[j].calConscore(RegionSize);
            clusters[j].sumOverallScore();
        }
        
//write cluster file
#ifdef CLUSTERLOG
        clusterFile<<"\n";
        for (int j=0; j<clusters.size(); j++){
            clusterFile<<setw(8)<<"CLUSTER\t"<<setw(3)<<j<<"\t";
            clusters[j].printMember(clusterFile);
        }
#endif
        
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
        
        //write pfm file;
        fileName = outPutDir + "/allmotif.pfm";
        ofstream pfmFile(fileName.c_str());
        if (!pfmFile) {
            printAndExit("Error! Fail to open pfmFile for writing!");
        }
        for (int i = 0; i<clusters.size(); i++) {
            printProgress(i,clusters.size(), "Generate PFM file");
            pfmFile<<clusters[i];
            
            
        }            
        //write pfm and dist
        
        tEnd=clock();
        cerr<<"total time eclapse:"<<double((tEnd-tStart)/1e6)<<endl;
        
        exit(0);
    }
    return 0;
}




void test(){
    //test locateSub
    /*
    float a[64]={222196,223883,224365,225399,226522,228458,230182,232294,234638,235955,237100,238276,240537,241488,244776,247591,250682,256220,261218,266507,273588,280507,288908,296303,305834,315449,326204,338612,351108,365119,380065,392915,407475,423412,438997,444292,434680,413637,384791,362230,363998,396845,451613,503894,542050,560474,556176,536605,515048,494379,477460,460567,442164,425085,411363,396023,381178,368042,356459,344080,334199,323655,314402,307031};
    //float a[64]={102.75,114.75,124.275,124.275,128,127.35,109.975,103.375,109.525,114.725,124.375,127.275,130.975,126.5,130.15,125.4,141.025,150.475,145.325,143.7,129.65,118.825,104.55,102.05,103.95,111.325,125.825,146.375,152.575,153.125,152.4,155.925,146.1,148.8,163.925,157.6,158.525,137.9,131.125,117.9,116.375,113.1,118.9,121.85,143.275,135.925,155.75,148.65,143.875,146.075,132,122.45,123.525,115.9,114.5,109.65,99.575,93.85,98.525,92.95,103.85,114.275,117.4,119.475};
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
    cerr<<"\n"<<"tagNoise"<<tempFFT.denoise(20)*NOISEWEIGHT;
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
     */
    
    //float fvec[65]={1.37328e-05,6.53941e-06,1.30788e-05,1.30788e-05,0.00022234,0.000268116,0.000251113,0.000245228,0.000288061,0.000306208,0.00027956,0.000291331,0.000235419,0.000249805,0.00024719,0.000232476,0.000251113,0.000255691,0.000256835,0.000305554,0.000314873,0.000518085,0.000557975,0.000603261,0.000775574,0.000794211,0.000812685,0.000797645,0.000743694,0.000814156,0.00076413,0.000725221,0.000744675,0.000657865,0.000679935,0.00073176,0.000633669,0.00060604,0.000524951,0.000410675,0.000354109,0.000292639,0.000260759,0.000248988,0.000327461,0.000362283,0.000376016,0.000482608,0.000475415,0.00044517,0.000505823,0.000432255,0.000471164,0.000398904,0.000324845,0.000314709,0.000279723,0.000241141,0.000207626,0.000196182,0.000188825,0.000121143,0.000122123,8.33775e-05,4.31601e-05};
    //vector<float> tempVec(fvec,fvec+4);
    /*
    FFT tempFFT(tempVec);
    cerr<<(tempFFT.denoise(int(SAMPLESIZE*5/8))*NOISEWEIGHT)<<endl;
    vector<float> tempBin;
    for (int i=0; i<SAMPLESIZE*2; i++) {
        tempBin.push_back(tempFFT.invTrans[i].real()>0?tempFFT.invTrans[i].real():0);
    }
    tempBin.push_back(fvec[64]);
    cerr<<tempBin<<endl;
    cerr<<testSymmety(tempBin)<<endl;
    
    for (vector<float>::iterator i=tempVec.begin();i!=tempVec.end();) {
        cerr<<*i<<endl;
        i++;
    }
    */
    
    //test pearson
    K=6;
    /*
    int a[4][6]={{1,2,3,4,5,6},
                 {2,3,4,5,6,7},
                 {3,2,3,4,5,6},
                 {4,2,3,4,5,6}};
     int b[4][6]={        {6,5,4,3,2,1},
     {6,5,4,3,2,1},
     {6,5,4,3,2,1},
     {6,5,4,3,2,1}};
    */
    int a[4][6]={{1,2,3,4,5,6},
        {1,2,3,4,5,6},
        {1,2,3,4,5,6},
        {1,2,3,4,5,8}};
    

    int b[4][6]={{8,5,4,3,2,1},
        {6,5,4,3,2,1},
        {6,5,4,3,2,1},
        {6,5,4,3,2,1}};
    
    
    Motif m1(1),m2(2);
    for (int j=0; j<6; j++) {
        for (int i=0; i<4; i++) {
            m1.pfm[i][j]=a[i][j];
            m2.pfm[i][j]=b[i][j];
        }
    }
    cerr<<m1.PearsonCorrPFM(m2, 0, 6, false)<<endl;
    assert(0==1);
}


