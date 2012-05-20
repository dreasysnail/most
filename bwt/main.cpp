//
//  main.cpp
//  bwt
//
//  Created by icarus on 11-10-4.
//  Copyright 2011年 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include "bwt.h"
#include "common.h"

using namespace std;

//global var's declarision
extern Edge Edges[ HASH_TABLE_SIZE ];
extern Node Nodes[ MAX_LENGTH * 2 ];
extern char T[ MAX_LENGTH ];
extern int N;
void printUsage();

int main(int argc, char **argv)
{
    Suffix active( 0, 0, -1 );  // The initial active prefix
    
    while (true) {
        if (argv[1][1]!='m'&&argc<4) {
            printUsage();
            exit(1);
        }
        if (argv[1][1]=='i'&&argc==4) {
            genomeRegions gR;
            gR.readBed(argv[2]);
            gR.readFasta(argv[3]);
            vector<string>::iterator it;
            string tempString;
            for (it=gR.genomeSeqs.begin(); it!=gR.genomeSeqs.end(); it++) {
                tempString += (*it)+"#";
                (*it).clear();
            }
            cout<<tempString<<endl;
            strcpy(T, tempString.c_str());
            tempString.clear();
            N = strlen( T ) - 1;
            for ( int i = 0 ; i <= N ; i++ )
                AddPrefix( active, i );
            cout<<active.countString("tcc")<<endl;
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
                AddPrefix( active, i );
            
            
            //display
            for(int i=0;i<Node::Count;i++){
                cout<<i<<"   "<<Nodes[i].father<<"   "<<Nodes[i].leaf_index<<"   "<<Nodes[i].suffix_node<<"    "<<Nodes[i].leaf_count_beneath<<endl;
                //         Nodes[i].showNodeString();
                
                
                
                
                
            }
            cout<<active.countString("MAD")<<endl;
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

