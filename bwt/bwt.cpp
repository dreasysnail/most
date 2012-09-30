//
//  bwt.cpp
//  MOST
//
//  Created by zhang yizhe on 12-5-20.
//  Copyright (c) 2012年 SJTU. All rights reserved.
//
//
// Suffix tree creation
//
// We thank Mark Nelson for his open source
// implementation of Ukekknon algorithm!
//
//

#include "bwt.h"
extern int RegionSize;

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
string T;
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
        long int j = i;
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
    
    return start_node;
    //
}



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




/******   other custom function   *******/

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
                        loci.push_back(RegionSize - (current_query_index+tempEdge.last_char_index-tempEdge.first_char_index+1));
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
                loci.push_back(RegionSize - (offset+edge.last_char_index-edge.first_char_index+1));
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



/******** I/O methods ********/


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




