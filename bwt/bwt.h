//
//  bwt.h
//  bwt
//
//  Created by zhang yizhe on 12-5-15.
//  Copyright (c) 2012年 SJTU. All rights reserved.
//

#pragma once
#ifndef BWT_T
#define BWT_T

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <string>
#include <algorithm>
#include "common.h"
#include "motif.h"
#include "FFT.h"
using std::cout;
using std::cin;
using std::cerr;
using std::setw;
using std::flush;
using std::endl;
using std::string;
using std::getline;
using std::istream;
using std::ostream;
using std::min;


class Node;
class Edge;
class Motif;

class Suffix {

public :
    int origin_node;
    int first_char_index;
    int last_char_index;
    int stringcount;
    Suffix( int node, int start, int stop)
    : origin_node( node ),
    first_char_index( start ),
    last_char_index( stop )
    {stringcount=1;};
    int Explicit(){ return last_char_index < first_char_index; }
    int Implicit(){ return last_char_index >= first_char_index; }
    void Canonize();

    //my custom function
    //restore searched motif-loci
    std::vector<int> loci;
    
    void AddPrefix(int last_char_index );
    int countString(const string &query );
    bool isExistString(const string &query);
    bool initialize();
    //static Node (*Nodes);
    //static Edge (*Edges);
    void AddSuffixLink( int &last_parent, int parent); 
    vector<int> locateMotif(Motif& currentMotif);
    void traverseLoci(int offset, int nodeIndex);

};


class Edge {
public :
    int first_char_index;
    int last_char_index;
    int end_node;
    int start_node;
    void Insert();
    void Remove();
    Edge():start_node(-1){};
    Edge( int init_first_char_index,
         int init_last_char_index,
         int parent_node );
    int SplitEdge( Suffix &s );
    static Edge Find( int node, int c );
    static long int Hash( int node, int c );
    

    //   static Node (*Nodes);
    //static Edge (*Edges);
};

class Node {
    
public :
    int suffix_node;
    int father;
    int leaf_index;
    Node() { suffix_node = -1;
        father=-1;
        leaf_index=-1;}
    void showNodeString();
    static int Count;
    static int Leaf;
    
    //my custome variance

    

    
    //    static Node (*Nodes);
    //static Edge (*Edges);
    
};



class Aux {
    public :
    int i;
    Aux( int rhs ){ i = rhs; }
    operator int(){ return i; }
};


class Buffer {
    public :
    char data[ MAX_LENGTH ];
    int N;
    Aux operator[]( int size ) const;
};

inline Aux Buffer::operator[]( int i ) const
{
    if ( i >= N )
        return Aux( 256 );
    else
        return Aux( data[ i ] );
}





//
// Necessary forward references
//
void validate();
int walk_tree( int start_node, int last_char_so_far );




//
void dump_edges( int current_n );
void AddPrefix( Suffix &active, int last_char_index );
ostream &operator<<( ostream &s, const Suffix &str );
ostream &operator<<( ostream &s, Aux &a );
ostream &operator<<( ostream &s, const Edge &edge );
ostream &operator<<( ostream &s, const Node &node );
istream &operator>>( istream &s, Buffer &b );
void print_parents( ostream &s, int node);

#endif
