//
//  bwt.h
//  bwt
//
//  Created by zhang yizhe on 12-5-15.
//  Copyright (c) 2012年 SJTU. All rights reserved.
//

#pragma once

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <string>
#include <algorithm>
#include "common.h"
#include "motif.h"
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


//
// The maximum input string length this program
// will handle is defined here.  A suffix tree
// can have as many as 2N edges/nodes.  The edges
// are stored in a hash table, whose size is also
// defined here.
//



//
// When a new tree is added to the table, we step
// through all the currently defined suffixes from
// the active point to the end point.  This structure
// defines a Suffix by its final character.
// In the canonical representation, we define that last
// character by starting at a node in the tree, and
// following a string of characters, represented by
// first_char_index and last_char_index.  The two indices
// point into the input string.  Note that if a suffix
// ends at a node, there are no additional characters
// needed to characterize its last character position.
// When this is the case, we say the node is Explicit,
// and set first_char_index > last_char_index to flag
// that.
// 
// from first_c_i to last_c_i :current implicit one expect to be inserted.

class Node;
class Edge;

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

//
// The suffix tree is made up of edges connecting nodes.
// Each edge represents a string of characters starting
// at first_char_index and ending at last_char_index.
// Edges can be inserted and removed from a hash table,
// based on the Hash() function defined here.  The hash
// table indicates an unused slot by setting the
// start_node value to -1.
//

class Edge {
public :
    int first_char_index;
    int last_char_index;
    int end_node;
    int start_node;
    void Insert();
    void Remove();
    Edge();
    Edge( int init_first_char_index,
         int init_last_char_index,
         int parent_node );
    int SplitEdge( Suffix &s );
    static Edge Find( int node, int c );
    static long int Hash( int node, int c );
    

    //   static Node (*Nodes);
    //static Edge (*Edges);
};

//
//  The only information contained in a node is the
//  suffix link. Each suffix in the tree that ends
//  at a particular node can find the next smaller suffix
//  by following the suffix_node link to a new node.  Nodes
//  are stored in a simple array.
//
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
#ifdef display
      int my_node_index;    
      int first_char_index;
      int last_char_index;
      int above_edge_first_char_index;
#endif
    int leaf_count_beneath;
    
    

    
    //    static Node (*Nodes);
    //static Edge (*Edges);
    
};


// The Buffer class exists purely to overload operator[],
// which allows me to return 256 for T[ N ].  Note also
// that operator[] doesn't exactly return an int or a char
// like you might think.  Instead, it returns an Aux object,
// which is just really a wrapper around an integer.  This
// lets me write an operator<<() for Aux and ostream so that
// outputting 256 will actually print "<EOF>"
//

//
// Since Aux is just a wrapper around an integer, all I
// need is a holder for the integer, a constructor,
// and a casting operator.
//

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
// The default ctor for Edge just sets start_node
// to the invalid value.  This is done to guarantee
// that the hash table is initially fildled with unused
// edges.
//


//
void dump_edges( int current_n );
void AddPrefix( Suffix &active, int last_char_index );
ostream &operator<<( ostream &s, const Suffix &str );
ostream &operator<<( ostream &s, Aux &a );
ostream &operator<<( ostream &s, const Edge &edge );
ostream &operator<<( ostream &s, const Node &node );
istream &operator>>( istream &s, Buffer &b );
void print_parents( ostream &s, int node);


