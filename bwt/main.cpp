//
//  main.cpp
//  bwt
//
//  Created by icarus on 11-10-4.
//  Copyright 2011年 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
using namespace std;

int main()
{	string sraw,sout;
	int b;
	while(cin>>sraw){
        sout=encode(sraw);
        cout<<"encode:"<<sout<<endl;
    }
	cin>>b;
    
}




