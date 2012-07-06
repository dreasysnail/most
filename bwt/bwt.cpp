//
//  bwt.cpp
//  bwt
//
//  Created by zhang yizhe on 12-5-20.
//  Copyright (c) 2012年 SJTU. All rights reserved.
//
//
// Suffix tree creation
//
// we thank Mark Nelson for his open source
// implement of Ukekknon algorithm
//
//

#include "bwt.h"
extern int GenomeSize;

//
// This is the hash table where all the currently
// defined edges are stored.  
//

Edge *Edges;

// This is all nodes
Node *Nodes;

//
// The array of defined nodes.  The count is 1 at the
// start because the initial tree has the root node
// defined, with no children.
//

int Node::Count = 1;
int Node::Leaf = 1;


//
// The input buffer and character count.  Please note that N
// is the length of the input string -1, which means it
// denotes the maximum index in the input buffer.
//

// This is the concatenated seq
char T[ MAX_LENGTH ];
int N;




Edge::Edge( int init_first, int init_last, int parent_node )
{
    first_char_index = init_first;
    last_char_index = init_last;
    start_node = parent_node;
    end_node = Node::Count++;
    Nodes[end_node].father=start_node;
    
    // Nodes[end_node].leaf_count_beneath=1;
    
    //Nodes[end_node].my_node_index=end_node;
}

//
// Edges are inserted into the hash table using this hashing
// function.
//

long int Edge::Hash( int node, int c )
{   int offset;
    switch (c) {
        case 'A':
            offset=1;
            break;
        case 'C':
            offset=2;
            break;
        case 'G':
            offset=3;
            break;
        case 'T':
            offset=4;
            break;
        case 'N':
            offset=5;
            break;
        case '#':
            offset=6;
            break;
        default:
            offset=c;
            break;
    }
    return ( ( node << 3 ) + offset  ) % HASH_TABLE_SIZE;
}

//
// A given edge gets a copy of itself inserted into the table
// with this function. 
//


void Edge::Insert()
{
    //truncate # and N
    //   if (start_node==0&&(T[ first_char_index ]=='#'||T[ first_char_index ]=='N'))        return;
    long int i = Hash( start_node, T[ first_char_index ] );
    while ( Edges[ i ].start_node != -1 )
        i = ++i % HASH_TABLE_SIZE;
    Edges[ i ] = *this;
    
}

//
// Removing an edge from the hash table
// Knuth, Sorting and Searching, Algorithm R, p. 527
//

void Edge::Remove()
{
    long int i = Hash( start_node, T[ first_char_index ] );
    while ( Edges[ i ].start_node != start_node ||
           Edges[ i ].first_char_index != first_char_index )
        i = ++i % HASH_TABLE_SIZE;
    for ( ; ; ) {
        Edges[ i ].start_node = -1;
        int j = i;
        for ( ; ; ) {
            i = ++i % HASH_TABLE_SIZE;
            if ( Edges[ i ].start_node == -1 )
                return;
            long int r = Hash( Edges[ i ].start_node, T[ Edges[ i ].first_char_index ] );
            if ( i >= r && r > j )
                continue;
            if ( r > j && j > i )
                continue;
            if ( j > i && i >= r )
                continue;
            break;
        }
        Edges[ j ] = Edges[ i ];
    }
}

//
// find a particular edge leading a particular node
//

Edge Edge::Find( int node, int c )
{
    long int i = Hash( node, c );
    for ( ; ; ) {
        if ( Edges[ i ].start_node == node ) 
//poisonous when string changed, just take for granted that no collision occurs
            if ( c == T[ Edges[ i ].first_char_index ] )
                return Edges[ i ];
        if ( Edges[ i ].start_node == -1 )
            return Edges[ i ];
        i = ++i % HASH_TABLE_SIZE;
    }
}

//
// This function is called to split an edge at the point defined by the Suffix argument.
//

int Edge::SplitEdge( Suffix &s )
{
    Remove();
    //act to split point
    Edge newEdge(  first_char_index,
                 first_char_index + s.last_char_index - s.first_char_index,
                 s.origin_node );
    Edge *new_edge = &newEdge;
    new_edge->Insert();
    
    
    
    

    //
    
    
    
    
    Nodes[ new_edge->end_node ].suffix_node = s.origin_node;
    first_char_index += s.last_char_index - s.first_char_index + 1;
    start_node = new_edge->end_node;
    Insert();
    Nodes[end_node].father=start_node;
    
        
 
    //update leaf counts from current node to 0
    
    /*
      Nodes[new_edge->end_node].leaf_count_beneath++;
    int fatherNode=s.origin_node;
    
    while (fatherNode!=0) {
        Nodes[fatherNode].leaf_count_beneath++;
        fatherNode = Nodes[fatherNode].father;
    }
    */
    
    return start_node;
    //
}


/*

void dump_edges( int current_n )
{
    cout << " Start  End  Suf  First Last  String\n";
    for ( int j = 0 ; j < HASH_TABLE_SIZE ; j++ ) {
        Edge *s = Edges + j;
        if ( s->start_node == -1 )
            continue;
        cout << setw( 5 ) << s->start_node << " "
        << setw( 5 ) << s->end_node << " "
        << setw( 3 ) << Nodes[ s->end_node ].suffix_node << " "
        << setw( 5 ) << s->first_char_index << " "
        << setw( 6 ) << s->last_char_index << "  ";
        int top;
        if ( current_n > s->last_char_index )
            top = s->last_char_index;
        else
            top = current_n;
        for ( int l = s->first_char_index ;
             l <= top;
             l++ )
            cout << T[ l ];
        cout << "\n";
    }
}
*/

//
// The canonical representation of a suffix for this algorithm
// requires that the origin_node by the closest node to the end
// of the tree. 
//

void Suffix::Canonize()
{
    if ( !Explicit() ) {
        Edge edge = Edge::Find( origin_node, T[ first_char_index ] );
        int edge_span = edge.last_char_index - edge.first_char_index;
        while ( edge_span <= ( last_char_index - first_char_index ) ) {
            first_char_index = first_char_index + edge_span + 1;
            origin_node = edge.end_node;
            if ( first_char_index <= last_char_index ) {
                edge = Edge::Find( edge.end_node, T[ first_char_index ] );
                edge_span = edge.last_char_index - edge.first_char_index;
            };
        }
    }
}

//
// heart of the algorithm.
//

void Suffix::AddPrefix(int current_index )
{
    int parent_node;
    int last_parent_node = -1;
    
    for ( ; ; ) {
        Edge edge;
        parent_node = origin_node;
        //
        // Step 1 is to try and find a matching edge for the given node.
        // If a matching edge exists, we are done adding edges, so we break
        // out of this big loop.
        //
        if ( Explicit() ) {
            edge = Edge::Find(origin_node, T[ current_index ] );
            if ( edge.start_node != -1 )
                break;
        } 
        else { //implicit node, a little more complicated
            edge = Edge::Find( origin_node, T[ first_char_index ] );
            int span =last_char_index - first_char_index;
            //if still match
            if ( T[ edge.first_char_index + span + 1 ] == T[ current_index ] )
                break;
            parent_node = edge.SplitEdge( *this );
        }
        //
        // We didn't find a matching edge, so we create a new one, add
        // it to the tree at the parent node position, and insert it
        // into the hash table.  When we create a new node, it also
        // means we need to create a suffix link to the new node from
        // the last node we visited.
        //
        
        
        Edge newEdge( current_index, N, parent_node );
        Edge *new_edge = &newEdge;
        new_edge->Insert();
        Nodes[new_edge->end_node].leaf_index=Node::Leaf++;
        
        //update leaf_count_beneath for upper node(icarus)
        /*
        int fatherNode = parent_node;
        while (fatherNode!=0) {
            Nodes[fatherNode].leaf_count_beneath++;
            fatherNode = Nodes[fatherNode].father;
        }
        */
        
        if ( last_parent_node > 0 )
            Nodes[ last_parent_node ].suffix_node = parent_node;
        last_parent_node = parent_node;
        //
        // This final step is where we move to the next smaller suffix  abc bc c
        //
        if ( origin_node == 0 )
            first_char_index++;
        else
            origin_node = Nodes[ origin_node ].suffix_node;
       Canonize();
    }
    if ( last_parent_node > 0 )
        Nodes[ last_parent_node ].suffix_node = parent_node;
    last_char_index++;  //Now the endpoint is the next active point
    Canonize();
}


//
// The validation code 
//
/*
char CurrentString[ MAX_LENGTH ];
char GoodSuffixes[ MAX_LENGTH ];
char BranchCount[ MAX_LENGTH * 2 ] = { 0 };

void validate()
{
    for ( int i = 0 ; i < N ; i++ )
        GoodSuffixes[ i ] = 0;
    walk_tree( 0, 0 );
    int error = 0;
    for ( int i = 0 ; i < N ; i++ )
        if ( GoodSuffixes[ i ] != 1 ) {
            cout << "Suffix " << i << " count wrong!\n";
            error++;
        }
    if ( error == 0 )
        cout << "All Suffixes present!\n";
    int leaf_count = 0;
    int branch_count = 0;
    for (int i = 0 ; i < Node::Count ; i++ ) {
        if ( BranchCount[ i ] == 0 )
            cout << "Logic error on node "
            << i
            << ", not a leaf or internal node!\n";
        else if ( BranchCount[ i ] == -1 )
            leaf_count++;
        else
            branch_count += BranchCount[ i ];
    }
    cout << "Leaf count : "
    << leaf_count
    << ( leaf_count == ( N + 1 ) ? " OK" : " Error!" )
    << "\n";
    cout << "Branch count : "
    << branch_count
    << ( branch_count == (Node::Count - 1) ? " OK" : " Error!" )
    << endl;
}

int walk_tree( int start_node, int last_char_so_far )
{
    int edges = 0;
    for ( int i = 0 ; i < 256 ; i++ ) {
        Edge edge = Edge::Find( start_node, i );
        if ( edge.start_node != -1 ) {
            if ( BranchCount[ edge.start_node ] < 0 )
                cerr << "Logic error on node "
                << edge.start_node
                << '\n';
            BranchCount[ edge.start_node ]++;
            edges++;
            int l = last_char_so_far;
            for ( int j = edge.first_char_index ; j <= edge.last_char_index ; j++ )
                CurrentString[ l++ ] = T[ j ];
            CurrentString[ l ] = '\0';
            if ( walk_tree( edge.end_node, l ) ) {
                if ( BranchCount[ edge.end_node ] > 0 )
                    cerr << "Logic error on node "
                    << edge.end_node
                    << "\n";
                BranchCount[ edge.end_node ]--;
            }
        }
    }
    //
    // If this node didn't have any child edges, it means we
    // are at a leaf node, and can check on this suffix.  We
    // check to see if it matches the input string, then tick
    // off it's entry in the GoodSuffixes list.
    //
    if ( edges == 0 ) {
        cout << "Suffix : ";
        for ( int m = 0 ; m < last_char_so_far ; m++ )
            cout << CurrentString[ m ];
        cout << "\n";
        GoodSuffixes[ strlen( CurrentString ) - 1 ]++;
        cout << "comparing: " << ( T + N - strlen( CurrentString ) + 1 )
        << " to " << CurrentString << endl;
        if ( strcmp(T + N - strlen(CurrentString) + 1, CurrentString ) != 0 )
            cout << "Comparison failure!\n";
        return 1;
    } else
        return 0;
}
*/

/******   other custom function   *******/

void Node::showNodeString(){
    //        std::string nodeString(&T[first_char_index],last_char_index-first_char_index);
    //    cout<<"\t"<<nodeString<<"\t"<<endl;
    
}

/*
int Suffix::countString(const string &query ){
    
    int currentNode = 0;
    int current_query_index = 0;
    int queryTemp;
    
    while (current_query_index<query.size()) {
        Edge tempEdge = Edge::Find(currentNode,query[current_query_index]);
        // cout<<tempEdge<<endl;
        //if find
        if (tempEdge.start_node!=-1) {
            string edgeString(&T[tempEdge.first_char_index],tempEdge.last_char_index-tempEdge.first_char_index+1);
//          cout<<T[tempEdge.first_char_index]<<edgeString<<endl;
            queryTemp = query.size()-current_query_index;
            if (edgeString.size()>=queryTemp){
                if (query.substr(current_query_index,queryTemp)!=edgeString.substr(0,queryTemp)) {
                    return 0;
                }
                else {
                    return Nodes[tempEdge.end_node].leaf_count_beneath;
                }
            }
            else {
                if (query.substr(current_query_index,edgeString.size())!=edgeString) {
                    return 0;
                }
                currentNode = tempEdge.end_node;
                current_query_index = current_query_index + edgeString.size();
            }
        }
        else {
            return 0;
        }
    }
    return -1;
}
*/

vector<int> Suffix::locateMotif(Motif& currentMotif){
    //implement2 find from edges walk tree
    //currentMotif.explainMotif();
    loci.clear();
    //cout<<currentMotif.expMotifs<<endl;
    int currentNode = 0;
    int current_query_index = 0;
    int queryTemp;
    string query = currentMotif.query;
    while (current_query_index<query.size()) {
        Edge tempEdge = Edge::Find(currentNode,query[current_query_index]);
        // cout<<tempEdge<<endl;
        //if find
        if (tempEdge.start_node!=-1) {
            string edgeString(&T[tempEdge.first_char_index],tempEdge.last_char_index-tempEdge.first_char_index+1);
            //          cout<<T[tempEdge.first_char_index]<<edgeString<<endl;
            queryTemp = query.size()-current_query_index;
            if (edgeString.size()>=queryTemp){
                if (query.substr(current_query_index,queryTemp)!=edgeString.substr(0,queryTemp)) {
                    break;
                }
                else {
                    // return edgeString.size()==queryTemp?Nodes[tempEdge.end_node].leaf_count_beneath:Nodes[currentNode].leaf_count_beneath;
                    if (Nodes[tempEdge.end_node].leaf_index!=-1){
                        loci.push_back(GenomeSize - (current_query_index+tempEdge.last_char_index-tempEdge.first_char_index+1));
                        break;
                    }
                    //traverse beneath nodes
                    traverseLoci(current_query_index+edgeString.size(),tempEdge.end_node);
                    break;
                }                    
            }
            else {
                if (query.substr(current_query_index,edgeString.size())!=edgeString) {
                    break;
                }
                currentNode = tempEdge.end_node;
                current_query_index = current_query_index + edgeString.size();
            }
        }
        else {
            break;
        }
    }
    return loci;
        
}


void Suffix::traverseLoci(int offset, int nodeIndex){
    char token[6]={'A','T','C','G','N','#'};
    
    for (int i=0; i<6; i++) {
        Edge edge = Edge::Find(nodeIndex, token[i]);
        if (edge.start_node != -1){
            if (Nodes[edge.end_node].leaf_index!=-1){
                loci.push_back(GenomeSize - (offset+edge.last_char_index-edge.first_char_index+1));
                continue;
            }
            else {
                traverseLoci(offset+edge.last_char_index-edge.first_char_index+1,edge.end_node);
                
            }  
        }
    }    
    return;
}

    

bool Suffix::isExistString(const string &query){
    int currentNode = 0;
    int current_query_index = 0;
    
    
    while (current_query_index<query.size()) {
        Edge tempEdge(Edge::Find(currentNode,query[current_query_index]));
        // cout<<tempEdge<<endl;
        //if find
        if (tempEdge.start_node!=-1) {
            string edgeString(&T[tempEdge.first_char_index],tempEdge.last_char_index-tempEdge.first_char_index+1);
        //cout<<T[tempEdge.first_char_index]<<edgeString<<endl;
            int minTemp = min(edgeString.size(),query.size()-current_query_index);
            if (query.substr(current_query_index,minTemp)!=edgeString.substr(0,minTemp)) {
                return false;
            }
            currentNode = tempEdge.end_node;
            current_query_index = current_query_index + minTemp;
            
            
        }
        else {
            return false;
        }
    }
    return true;

}


bool Suffix::initialize(){
    first_char_index=0;
    last_char_index=-1;
    origin_node=0;
    stringcount++;
    return true;
}




//
// STREED.CPP - Suffix tree creation - debug version
//
// Mark Nelson, updated December, 2006
//
// This code has been tested with Borland C++ and
// Microsoft Visual C++.
//
// This program gets a line of input, either from the
// command line or from user input.  It then creates
// the suffix tree corresponding to the given text.
//
// This program is intended to be a supplement to the
// code found in STREE.CPP.  It contains a extensive
// debugging information, which might make it harder
// to read.
//
// This version of the program also gets around the
// problem of requiring the last character of the
// input text to be unique.  It does this by overloading
// operator[] for the input buffer object.  When you select
// T[ N ], you will get a value of 256, which is obviously
// going to be a unique member of the character string.
// This overloading adds some complexity, which just might
// make the program a little harder to read!
//
// In addition, there is some overloading trickery that lets
// you send T[i] to the output stream, and send the 256 value
// as the string "<EOF>".  Another convenience that adds
// code and complexity.
//



//
//
// The maximum input string length this program
// will handle is defined here.  A suffix tree
// can have as many as 2N edges/nodes.  The edges
// are stored in a hash table, whose size is also
// defined here.  When I want to exercise the hash
// table a little bit, I set MAX_LENGTH to 6 and
// HASH_TABLE_SIZE to 13.
//
//






//
// A suffix in the tree is denoted by a Suffix structure
// that denotes its last character.  The canonical
// representation of a suffix for this algorithm requires
// that the origin_node by the closest node to the end
// of the tree.  To force this to be true, we have to
// slide down every edge in our current path until we
// reach the final node.


//
// This debug routine prints out the value of a
// Suffix object.  In order to print out the
// entire suffix string, I have to walk up the
// tree to each of the parent nodes.  This is
// handled by the print_parents() routine, which
// does this recursively.
//



ostream &operator<<( ostream &s, const Suffix &str )
{
    s << "("
    << str.origin_node
    << ",("
    << str.first_char_index
    << ","
    << str.last_char_index
    << ") ";
    s << "\"";
    print_parents( s, str.origin_node);
    for ( int i = str.first_char_index ;
         i <= str.last_char_index ;
         i++ )
        s << T[ i ];
    s << "\"";
    s << ")";
    return s;
}


void print_parents( ostream &s, int node)
{
    if ( node != 0 )
        for ( int i = 0 ; i < HASH_TABLE_SIZE ; i++ ) {
            if ( Edges[ i ].end_node == node ) {
                print_parents( s, Edges[ i ].start_node);
                for ( int j = Edges[ i ].first_char_index ;
                     j <= Edges[ i ].last_char_index
                     ; j++ )
                    s << T[ j ];
                return;
            }
        }
}


istream &operator>>( istream &s, Buffer &b )
{
    s >> b.data;
    assert( strlen( b.data ) < MAX_LENGTH );
    b.N = strlen( b.data );
    return s;
}

ostream &operator<<( ostream &s, Aux &a )
{
    if ( a.i == 256 )
        s << "<EOF>";
    else
        s << (char) a.i;
    return s;
}

ostream &operator<<( ostream &s, const Edge &edge )
{
    s << "Start, end nodes= "
    << edge.start_node
    << ", "
    << edge.end_node
    << " first, last = "
    << edge.first_char_index
    << ", "
    << edge.last_char_index
    << " \"";
    for ( int i = edge.first_char_index ; i <= edge.last_char_index ; i++ )
        s << T[ i ];
    s << "\"";
    return s;
}



//
// Adding a suffix line in AddPrefix() is really
// a simple operation.  All that needs to be done
// is to write out the correct value to the Nodes[]
// table in the correct place.  Since I've
// added some debug code here, it made sense to
// move it to a separate routine, even though it
// isn't being done that way in STREE.CPP
//

void Suffix::AddSuffixLink( int &last_parent, int parent)
{
    if ( last_parent > 0 ) {
        cout << "Creating suffix link from node "
        << last_parent
        << " to node "
        << parent
        << ".\n";
        Nodes[ last_parent ].suffix_node = parent;
    }
    last_parent = parent;
}

