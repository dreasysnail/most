//
//  compress.cpp
//  bwt
//
//  Created by icarus on 11-10-19.
//  Copyright 2011年 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include <fstream>     //for fopen
using namespace std;

int main(int argc,char* argv[])
{
    FILE *fp;
    if (argc>2) {
        cout<<"too many parameter"<<endl;
    }
    else{
        fp = fopen(argv[1],"r");
    }
    char ch;
    fscanf(fp, "%c", ch);     //fget?   etc?
    cout<<ch<<endl;

}



