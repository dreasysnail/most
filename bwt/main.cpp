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
int temp;

int main(int argc, char **argv)
{
    clock_t tStart,t1,t2,t2_1,t3,t4,tEnd;
    tStart=clock();
    cerr<<"START MOTIF FINDING"<<endl;
    Edge *EdgesTemp = new Edge[ HASH_TABLE_SIZE ];
    Node *NodesTemp = new Node[ MAX_LENGTH * 2 ];
    genomeRegions *gR = new genomeRegions(0);
    
    preInitialization(EdgesTemp,NodesTemp);
    
    
    Suffix active( 0, 0, -1 );  // The initial active prefix
    while (true) {
        if (argc<2||(argv[1][1]!='m'&&argc<4)) {
            printUsage();
            exit(1);
        }

        // tag mode
        if (argv[1][1]=='t'&&argc==7) {
            /************************
            2:tss.bed
            3:fasta
            4:histone.wig/histone.bed
            6:output directory
            ************************/
            string outPutDir(argv[6]); 
            system(("rm -rf "+outPutDir).c_str());
            if (system(("mkdir "+outPutDir).c_str()) != 0)
            {
                cerr<<"cannot make directory"<<endl;
                exit(1);
            }
            //extend 300bp
            

            gR->readFasta(argv[3]);
            gR->readBed(argv[2]);
            gR->mergeOverlap();
            gR->getSeq(outPutDir);
            //region-wide
            //gR->initProb(2);
            //initiate T 
            cerr<<"regionFile（bedFormat）:"<<argv[2]<<endl;
            GenomeSize = gR->catSeq(T);
            t1=clock();
            cerr<<"Parsing Fasta:"<<double((t1-tStart)/1e6)<<endl;
            
            string tagFileName(argv[4]);
            //if tag file is wigfile
            if (tagFileName.substr(tagFileName.size()-3,3)=="wig") {
                cerr<<"tagfile（WigFormat）:"<<tagFileName<<endl;
                gR->readWig(tagFileName);
            }
            else {
                //problematic
                genomeRegions tag(0);
                cerr<<"tagfile（BedFormat）:"<<tagFileName<<endl;
                tag.readBed(tagFileName);
                gR->writeRawTag(tag);
            }
            //gR->getTagBed();
            t2=clock();
            cerr<<"Parsing Wig:"<<double((t2-t1)/1e6)<<endl;


            
            
            //ofstream tfile("aa.fa");
            //tfile<<T<<endl;
            //cerr<<T<<endl;
            //cout<<gR->genomeTags<<endl;
            for (int i=0; i<gR->tagName.size(); i++) {
                cerr<<gR->tagName[i]<<"size:"<<gR->genomeTags[gR->tagName[i]].size()<<"\t";
                //cerr<<gR->genomeTags[gR->tagName[i]]<<endl;
                assert(GenomeSize==gR->genomeTags[gR->tagName[i]].size());
            }
            cerr<<"genomesize:"<<GenomeSize<<endl;
            assert(GenomeSize<=MAX_LENGTH);
            assert(GenomeSize*2<=HASH_TABLE_SIZE);
            N = strlen(T) - 1;
            
            for ( int i = 0 ; i <= N ; i++ )
                active.AddPrefix(i);
            t2_1=clock();
            cerr<<"built suffix tree:"<<double((t2_1-t2)/1e6)<<endl;
            
            //count words
            vector<Motif> qualifiedMotifs;
            qualifiedMotifs.reserve(K_5/1000);
            int counter1 = 0;
            for (int i = 0; i <(K_5); i++) {
                printProgress(i,K_5,"Qualify kmer from suffix tree:");
                Motif thisMotif(i);
                if (!thisMotif.noWildcard()) {
                    continue;
                }
                
                //order null
                thisMotif.initProb((*gR),-1);
                thisMotif.loci = active.locateMotif(thisMotif,qualifiedMotifs);
                
                //cout<<allmotifs[i].loci.size()<<" "<<active.countString(allmotifs[i].query)<<endl;
                temp += thisMotif.loci.size();
                thisMotif.calConscore(GenomeSize);
                if (thisMotif.score&&!thisMotif.isRepeat()) {
                    counter1++;
                    thisMotif.initPWM();
                    thisMotif.expMotifs.push_back(qualifiedMotifs.size());
                    thisMotif.testMotifTag(*gR, outPutDir, false);
                    thisMotif.sumScore();
                    if (thisMotif.overallScore>MINOVERALLSCORE)
                        qualifiedMotifs.push_back(thisMotif);
                    else {
                        thisMotif.printMotif();
                    }
                }
            }
            cout<<"clustering result:";
            cerr<<"\n"<<"STAGE1:"<<pow1(4, K)-counter1<<" kmer was filtered"<<"\n";
            cerr<<"STAGE2:"<<counter1-qualifiedMotifs.size()<<" kmer was filtered"<<"\n";
            cerr<<"Left:"<<qualifiedMotifs.size()<<endl;
            cerr<<"total motifs'loci size:"<<temp<<" approximate "<<GenomeSize<<endl;
            //need to eliminate allmotif
            
            
            t3=clock();
            cerr<<"count words:"<<double((t3-t2_1)/1e6)<<endl;
            
            //clustering
            vector<Motif> originalQualified(qualifiedMotifs);
            sort(qualifiedMotifs.begin(), qualifiedMotifs.end(),compareMotif());
            
            
            //cout<<qualifiedMotifs.size()<<qualifiedMotifs[0]<<qualifiedMotifs[1]<<endl;
            vector<Motif> clusters;
            int maxMotifSize=(MOTIFMAX<qualifiedMotifs.size()?MOTIFMAX:qualifiedMotifs.size());
            clusters.push_back(qualifiedMotifs[0]);
            //for (int i=0; i<10; i++) qualifiedMotifs[i].printMotif();
            // pair<int,int> tempa = qualifiedMotifs[2].editDistance(qualifiedMotifs[4]);
            // cout<<tempa.first<<" "<<tempa.second<<endl;
            cerr<<"clustering:"<<endl;
            for (int i=1; i<maxMotifSize; i++) { 
                printProgress(i,maxMotifSize,"clustering:");
                int bound = clusters.size();
                for (int j=0; j<bound; j++) {
                    pair<int,int> dist_shift = qualifiedMotifs[i].editDistance(clusters[j]);
                    float signifDist = qualifiedMotifs[i].histoneDistrDistance(clusters[j]);
                    //if highest mark is antisense and matched current cluster
                    if (dist_shift.first<0&&(-dist_shift.first)<=MAXDISTANCE) {
                        goto nextMotif;
                    }
                    else if (dist_shift.first<=MAXDISTANCE&&signifDist<=MAXSIGNDIST) {
                        //cluster size exceed Max cluster size
                        int queryLength = clusters[j].query.size();
                        if (queryLength>=MAXCLUSTERSIZE&&(dist_shift.second<0||dist_shift.second+K>queryLength)) {
                            goto nextMotif;
                        }
                        //cout<<qualifiedMotifs[i].query<<qualifiedMotifs[i].probThresh<<" "<<qualifiedMotifs[i].loci.size()<<endl<<j<<" "<<clusters[j].query<<endl;
                        clusters[j].concatenate(qualifiedMotifs[i],i,dist_shift.second);
                        clusters[j].calPWM(qualifiedMotifs[i], dist_shift.second);
                        clusters[j].loci = active.locateMotif(clusters[j],originalQualified);
                        clusters[j].testMotifTag(*gR,outPutDir,true);
                        goto nextMotif;
                    }

                }
                //if not aligned to any cluster
                if (clusters.size()<CLUSTERMAX){
                    qualifiedMotifs[i].index = -1;
                    clusters.push_back(qualifiedMotifs[i]);
                }
            nextMotif:
                continue;
            }
            //calculating
            for (int i=0; i<clusters.size(); i++){
                /*
                if (clusters[i].implicit()) {
                    clusters[i].score = 0;
                    continue;
                }
                */
                clusters[i].trim();
                //cout<<clusters[i].loci.size()<<endl;
                clusters[i].loci = active.locateMotif(clusters[i],originalQualified);
                //cout<<clusters[i].loci.size()<<endl;
                //calculate score for each cluster
                clusters[i].fillMotif(originalQualified);
                clusters[i].testMotifTag(*gR,outPutDir,true);
                clusters[i].sumScore();
            }
            sort(clusters.begin(),clusters.end(),compareCluster());
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
            
            
            /*
            
            //pwm file;
            delete [] EdgesTemp;
            delete [] NodesTemp;
            
            cerr<<endl<<"count words:"<<double((t3-t2)/1e6)<<endl;
            vector<Motif> significantMotifs;
            
            //start calculating wildcards
            for (int i =0; i <(K_5); i++) {
                printProgress(i,K_5,"calculate wildcard motifs:");
                if (Motif::motif[i]!=0||i%5==0||allmotifs[i].noWildcard()||allmotifs[i].wildcardNum()>K/2||allmotifs[i].implicit()||allmotifs[i].isRepeat()) {
                    continue;
                }
                allmotifs[i].fillMotif(allmotifs);
                if (allmotifs[i].score) {
                    //loci for this motif
                    
                    
                    //allmotifs[i].locateMotif(T);
                    allmotifs[i].loci = active.locateMotif(allmotifs[i],allmotifs);
                    for (int j=0; j<gR->tagName.size(); j++) {
                        allmotifs[i].testMotifTag(gR->genomeTags[gR->tagName[j]]);
                    }
                    
                    allmotifs[i].printMotif();
                    if (allmotifs[i].signif[0]>SINGIFTHRESH&&allmotifs[i].score>SCORETHRESH) {
                        //calPWM
                        allmotifs[i].sumScore();
                        allmotifs[i].calPWM(allmotifs);
                        significantMotifs.push_back(allmotifs[i]);
                    }
                }
                
            }
            t4 = clock();
            cerr<<endl<<"add words:"<<double((t4-t3)/1e6)<<endl;
            sort(significantMotifs.begin(), significantMotifs.end(),compareMotif());
            fileName = outPutDir + "/motif.dist";
            ofstream distFile(fileName.c_str());
            cerr<<endl<<"Generate PWM and DIST file with top "<<significantMotifs.size()<<" motifs:"<<endl;
            
            delete gR;
           
             
            */
            

            return 1;
        }

        
        
        cout << "Enter string: " << flush;
        cin.getline( T, MAX_LENGTH - 1 );
        cout<<T<<endl;
        
        if (T[0]=='.')
            break;
        else {
            N = strlen( T ) - 1;
            //cout<<N<<endl;
            //
            // The active point is the first non-leaf suffix in the
            // tree.  We start by setting this to be the empty string
            // at node 0.  The AddPrefix() function will update this
            // value after every new prefix is added.
            // suffix is a tree
            //
            
            for ( int i = 0 ; i <= N ; i++ )
                active.AddPrefix( i );
            
            
            //display
            for(int i=0;i<Node::Count;i++){
                cout<<i<<"   "<<Nodes[i].father<<"   "<<Nodes[i].leaf_index<<"   "<<Nodes[i].suffix_node<<"    "<<endl;
                //         Nodes[i].showNodeString();
                
                
                
                
                
            }
            active.initialize();
            
            
            
            
            // leaf index means it's a leaf and its no. 
            // Once all N prefixes have been added, the resulting table
            // of edges is printed out, and a validation step is
            // optionally performed.
            //
            //   dump_edges( N );
            //   cout << "Would you like to validate the tree?"
            //        << flush;
            //   std::string s;
            //    std::getline( cin, s ); 
            //   if ( s.size() > 0 && s[ 0 ] == 'Y' || s[ 0 ] == 'y' )
            //       validate();

        }
    }
    //print info
    
    
    
   
       return 1;
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
	cout<<usage<<endl;
    
    cout << "Normally, suffix trees require that the last\n"
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
