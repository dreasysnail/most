//
//  Clustering.h
//  most
//
//  Created by icarus on 13-1-2.
//  Copyright (c) 2013年 SJTU. All rights reserved.
//

#ifndef __most__Clustering__
#define __most__Clustering__

#include <iostream>
#include "common.h"
using std::cout;
using std::endl;

struct TreeNode {
    int left;
    int right;
    double distance;
};

//polar coordinates
struct PolarPos {
    float angle;
    float length;
};

class compareNode {
public:
    bool operator()(TreeNode tn1, TreeNode tn2) {
        return tn1.distance<tn2.distance;
    }
};

vector<TreeNode> SingleLinkClstr (const matrix_t &OriginDistMat);
vector<TreeNode> AverageLinkClstr (const matrix_t &OriginDistMat);
vector<TreeNode> CompleteLinkClstr (const matrix_t &OriginDistMat);
static double find_closest_pair(int n, matrix_t &distmatrix, int* ip, int* jp);
vector<int> cutTree (vector<TreeNode> Tree, int nclusters);
        //method: 0 single-linkage 1 average-linkage 2 complete-linkage
int GapStatistics (matrix_t &distmatrix, vector<TreeNode> Tree,char method=1);
float CalLogWk(const vector<TreeNode> &tn, int k, const matrix_t &distmatrix);
int GetRepElement(const matrix_t &distmatrix, vector<int> clusterLabel, int clstrID);
vector<float> GetRowSum(const matrix_t &distmatrix);
bool IsMatrix(const matrix_t &distmatrix);


vector<matrix_t > SubMatrice(matrix_t matrix, vector<int> label);
float CalWithinGrpDist(matrix_t &distmatrix);
matrix_t DistMatrixOfUniformDispersion(int nPoints, float var, int seed);

#endif /* defined(__most__Clustering__) */
