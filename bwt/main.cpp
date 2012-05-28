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

using namespace std;

//global var's declarision

extern char T[ MAX_LENGTH ];
extern int N;
void printUsage();
void printProgress(const int i,const string& message);
extern int GenomeSize;
extern Edge* Edges;
extern Node* Nodes;
void preInitialization(Edge* Edges,Node* Nodes);
int temp;

int main(int argc, char **argv)
{
    clock_t tStart,t1,t2,t3,tEnd;
    tStart=clock();
    Edge *EdgesTemp = new Edge[ HASH_TABLE_SIZE ];
    Node *NodesTemp = new Node[ MAX_LENGTH * 2 ];
    
    preInitialization(EdgesTemp,NodesTemp);
    t1=clock();
    
    Suffix active( 0, 0, -1 );  // The initial active prefix
    while (true) {
        if (argc<2||(argv[1][1]!='m'&&argc<4)) {
            printUsage();
            exit(1);
        }
        if (argv[1][1]=='i'&&argc==4) {
/*
            genomeRegions gR(0);
            gR.readBed(argv[2]);
            gR.readFasta(argv[3]);
            vector<string>::iterator it;
            string tempString;
            for (it=gR.genomeSeqs.begin(); it!=gR.genomeSeqs.end(); it++) {
                if ((*it)=="") {
                    continue;
                }
                tempString += (*it)+"#";
                transform(tempString.begin(), tempString.end(), tempString.begin(), ::toupper);
                GenomeSize += (*it).size();
                (*it).clear();
            }
           
            //tempString = trimN(tempString);
            strcpy(T,(tempString).c_str());
            tempString.clear();
            //cout<<T<<endl;
            N = strlen(T) - 1;
            
            for ( int i = 0 ; i <= N ; i++ )
                active.AddPrefix(i);
            for (int i = 0; i <(K_5); i++) {
                string query=translate(i);
                motif[i]=active.countString(query);
                 temp+=motif[i];
                printProgress(i,"find kmer from suffix tree:");
            }
            cout<<temp<<endl;
            for (int i =0; i <(K_5); i++) {
                if (motif[i]!=0||i%5==0) {
                    continue;
                }
                
                float score = fillMotif(i);
                if (score) {
                    cout<<"motif:"<<translate(i)<<"\tscore:"<<score<<endl;
                }
            // printProgress(i,"calculate motifs");
            }
            delete [] EdgesTemp;
            delete [] NodesTemp;
            return 1;
*/
        }
        
        
        
        
        // tag mode
        if (argv[1][1]=='t'&&argc==5) {
            //extend 300bp
            cerr<<"initialize:"<<double((t1-tStart)/1e6)<<endl;
            genomeRegions gR(0);
            gR.readBed(argv[2]);
            gR.readFasta(argv[3]);
            genomeRegions tag(0);
            tag.readBed(argv[4]);
            gR.writeRawTag(tag);
            gR.getTagBed();
            
            //initiate T and 
            vector<string>::iterator it;
            string tempString;
            for (it=gR.genomeSeqs.begin(); it!=gR.genomeSeqs.end(); it++) {
                if ((*it)=="") {
                    continue;
                }
                tempString += (*it)+"#";
                transform(tempString.begin(), tempString.end(), tempString.begin(), ::toupper);
                (*it).clear();
            }
            tempString += antisense(tempString)+"#";
            strcpy(T,(tempString).c_str());
            GenomeSize=tempString.size();
            tempString.clear();
            //cout<<T<<endl;
            //cout<<gR.genomeTags<<endl;
            
            int tempGenomeSize=gR.genomeTags.size();
            cout<<"tagsize:"<<tempGenomeSize<<"\tgenomesize:"<<GenomeSize<<endl;
            assert(GenomeSize==tempGenomeSize);
            assert(GenomeSize<=MAX_LENGTH);
            assert(GenomeSize*2<=HASH_TABLE_SIZE);
            N = strlen(T) - 1;
            t2=clock();
            
            cerr<<"built suffix tree:"<<double((t2-t1)/1e6)<<endl;
            vector<Motif> allmotifs;
            for ( int i = 0 ; i <= N ; i++ )
                active.AddPrefix(i);
            for (int i = 0; i <(K_5); i++) {
                Motif thisMotif(i);
                allmotifs.push_back(thisMotif);
                Motif::motif[i]=active.countString(thisMotif.query);
                temp += Motif::motif[i];
                printProgress(i,"find kmer from suffix tree:");
            }
            cout<<temp<<endl;
            
            
            t3=clock();
            cerr<<"count words:"<<double((t3-t2)/1e6)<<endl;
            cout<<"motif"<<"\tCONscore"<<"\tTAGscore"<<endl;
            for (int i =0; i <(K_5); i++) {
                if (Motif::motif[i]!=0||i%5==0) {
                    continue;
                }
                allmotifs[i].fillMotif();
                if (allmotifs[i].score) {
                    //loci for this motif
                    
                    //allmotifs[i].locateMotif(T);
                    allmotifs[i].loci = active.locateMotif(allmotifs[i]);
                    allmotifs[i].testMotifTag(gR.genomeTags);
                    allmotifs[i].printMotif();
                }
                // printProgress(i,"calculate motifs");
            }
            tEnd=clock();
            cerr<<"initialize:"<<double((t1-tStart)/1e6)<<endl;
            cerr<<"built suffix tree:"<<double((t2-t1)/1e6)<<endl;
            cerr<<"count words:"<<double((t3-t2)/1e6)<<endl;
            cerr<<"add words:"<<double((tEnd-t3)/1e6)<<endl;
            
            delete [] EdgesTemp;
            delete [] NodesTemp;
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
                cout<<i<<"   "<<Nodes[i].father<<"   "<<Nodes[i].leaf_index<<"   "<<Nodes[i].suffix_node<<"    "<<Nodes[i].leaf_count_beneath<<endl;
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
