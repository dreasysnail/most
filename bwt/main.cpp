//
//  main.cpp
//  bwt
//
//  Created by icarus on 11-10-4.
//  Copyright 2011年 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include "bwt.h"

using namespace std;

//global var's declarision
extern Edge Edges[ HASH_TABLE_SIZE ];
extern Node Nodes[ MAX_LENGTH * 2 ];
extern char T[ MAX_LENGTH ];
extern int N;

int main()
{
    
    cout << "Normally, suffix trees require that the last\n"
    << "character in the input string be unique.  If\n"
    << "you don't do this, your tree will contain\n"
    << "suffixes that don't end in leaf nodes.  This is\n"
    << "often a useful requirement. You can build a tree\n"
    << "in this program without meeting this requirement,\n"
    << "but the validation code will flag it as being an\n"
    << "invalid tree\n\n";
    
    
    Suffix active( 0, 0, -1 );  // The initial active prefix
    
    while (true) {
        cout << "Enter string: " << flush;
        cin.getline( T, MAX_LENGTH - 1 );
        if (T[0]=='.')
            break;
        else {
            N = strlen( T ) - 1;
            //
            // The active point is the first non-leaf suffix in the
            // tree.  We start by setting this to be the empty string
            // at node 0.  The AddPrefix() function will update this
            // value after every new prefix is added.
            // suffix is a tree
            //
            
            for ( int i = 0 ; i <= N ; i++ )
                AddPrefix( active, i );
            
            
            
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
    for(int i=0;i<Node::Count;i++){
        cout<<i<<"   "<<Nodes[i].father<<"   "<<Nodes[i].leaf_index<<"   "<<Nodes[i].first_char_index<<"   "<<Nodes[i].last_char_index<<"  "<<Nodes[i].leaf_count_beneath<<endl;
        //  Nodes[i].showNodeString();
        
        
        
        
        
    }
    
    cout<<active.countString("hah")<<endl;
   
       return 1;
};

