//
//  Clustering.cpp
//  most
//
//  Created by icarus on 13-1-2.
//  Copyright (c) 2013年 SJTU. All rights reserved.
//  Reference:
//  1. Sibson, R. (1973). SLINK: An optimally efficient algorithm for the single-link cluster method. The Computer Journal, 16(1): 30-34.
//  2. M. J. L. de Hoon, S. Imoto, J. Nolan, and S. Miyano: Open Source Clustering Software. Bioinformatics, 20 (9): 1453--1454 (2004).
#include "Clustering.h"
#include "time.h"
#include <cmath>
#include <set>
using std::max;


/*
 
 Purpose
 =======
 
 This routine performs single-linkage hierarchical clustering, using
 either the distance matrix directly, if available, or by calculating the
 distances from the data array. 
 The output of this algorithm is identical to conventional single-linkage
 hierarchical clustering, but is much more memory-efficient and faster. Hence,
 it can be applied to large data sets, for which the conventional single-
 linkage algorithm fails due to lack of memory.
 
 
 Arguments
 =========
 
 distmatrix (input)
 The distance matrix. If the distance matrix is passed by the calling routine
 treecluster, it is used by pslcluster to speed up the clustering calculation.
 The pslcluster routine does not modify the contents of distmatrix, and does
 not deallocate it. If distmatrix is NULL, the pairwise distances are calculated
 
 Return value
 ============
 
 An array of TreeNode structs, describing the
 hierarchical clustering solution consisting of nelements-1 nodes.
 
 negative left/right means node
 
 ========================================================================
*/
vector<TreeNode> SingleLinkClstr (const matrix_t &OriginDistMat) {
    matrix_t distmatrix = OriginDistMat;
    //element: leaf   node: non-leaf nodes
    const int nelements = int(distmatrix.size());
    const int nnodes = nelements - 1;

    
    
    TreeNode tempT = {0,0,0};
    vector<TreeNode> result;
    result.assign(nelements,tempT);
    
    vector<float> temp;
    vector<int> index, activeVector;
    temp.assign(nnodes,0);
    index.assign(nelements,0);
        
    //vector: 0,1,2...nNode
    for (int i = 0; i < nnodes; i++)
        activeVector.push_back(i);
    
    //   if(distmatrix)
    for (int i = 0; i < nelements; i++) {
        //current A: i
        result[i].distance = BIGNUM;
        for (int j = 0; j < i; j++) temp[j] = distmatrix[i][j];
        //
        for (int j = 0; j < i; j++) {
            //current B: K
            int k = activeVector[j];
            if (result[j].distance >= temp[j]) {
                if (result[j].distance < temp[k]) temp[k] = result[j].distance;
                result[j].distance = temp[j];
                activeVector[j] = i;
            }
            else if (temp[j] < temp[k])
                temp[k] = temp[j];
            
        }
        for (int j = 0; j < i; j++) {
            if (result[j].distance >= result[activeVector[j]].distance) activeVector[j] = i;
        }
    }
    
    for (int i = 0; i < nnodes; i++) result[i].left = i;
    sort(result.begin(),result.end(),compareNode());
    
    for (int i = 0; i < nelements; i++) index[i] = i;
    for (int i = 0; i < nnodes; i++) {
        int j = result[i].left;
        int k = activeVector[j];
        result[i].left = index[j];
        result[i].right = index[k];
        index[k] = -i-1;
    }
    //negative left/right means node rather than leaf
    return result;
}


vector<TreeNode> AverageLinkClstr (const matrix_t &OriginDistMat) {
        matrix_t distmatrix = OriginDistMat;
        const int nelements = int(distmatrix.size());
        const int nnodes = nelements - 1;
        int j;
        int n;
    
        TreeNode tempT = {0,0,0};
        vector<TreeNode> result;
        result.assign(nelements,tempT);
    
        vector<int> clusterid,number;
        clusterid.assign(nelements,0);
        number.assign(nelements,0);
 

        
        /* Setup a list specifying to which cluster a gene belongs, and keep track
         * of the number of elements in each cluster (needed to calculate the
         * average). */
        for (j = 0; j < nelements; j++) {
            number[j] = 1;
            clusterid[j] = j;
        }
        
        for (n = nelements; n > 1; n--) {
            int sum;
            int is = 1;
            int js = 0;
            result[nelements-n].distance = find_closest_pair(n, distmatrix, &is, &js);
            
            /* Save result */
            result[nelements-n].left = clusterid[is];
            result[nelements-n].right = clusterid[js];
            
            /* Fix the distances */
            sum = number[is] + number[js];
            for (j = 0; j < js; j++)
            { distmatrix[js][j] = distmatrix[is][j]*number[is]
                + distmatrix[js][j]*number[js];
                distmatrix[js][j] /= sum;
            }
            for (j = js+1; j < is; j++)
            { distmatrix[j][js] = distmatrix[is][j]*number[is]
                + distmatrix[j][js]*number[js];
                distmatrix[j][js] /= sum;
            }
            for (j = is+1; j < n; j++)
            { distmatrix[j][js] = distmatrix[j][is]*number[is]
                + distmatrix[j][js]*number[js];
                distmatrix[j][js] /= sum;
            }
            
            for (j = 0; j < is; j++) distmatrix[is][j] = distmatrix[n-1][j];
            for (j = is+1; j < n-1; j++) distmatrix[j][is] = distmatrix[n-1][j];
            
            /* Update number of elements in the clusters */
            number[js] = sum;
            number[is] = number[n-1];
            
            /* Update clusterids */
            clusterid[js] = n-nelements-1;
            clusterid[is] = clusterid[n-1];
        }

        return result;
}

vector<TreeNode> CompleteLinkClstr (const matrix_t &OriginDistMat) {
    matrix_t distmatrix = OriginDistMat;
    const int nelements = int(distmatrix.size());
    const int nnodes = nelements - 1;
    int j;
    int n;
    
    TreeNode tempT = {0,0,0};
    vector<TreeNode> result;
    result.assign(nelements,tempT);
    
    vector<int> clusterid;
    clusterid.assign(nelements,0);

    
    
    
    /* Setup a list specifying to which cluster a gene belongs */
    for (j = 0; j < nelements; j++) {
        clusterid[j] = j;
    }
    
    for (n = nelements; n > 1; n--) {
        int is = 1;
        int js = 0;
        result[nelements-n].distance = find_closest_pair(n, distmatrix, &is, &js);
        
        /* Fix the distances */
        for (j = 0; j < js; j++)
            distmatrix[js][j] = max(distmatrix[is][j],distmatrix[js][j]);
        for (j = js+1; j < is; j++)
            distmatrix[j][js] = max(distmatrix[is][j],distmatrix[j][js]);
        for (j = is+1; j < n; j++)
            distmatrix[j][js] = max(distmatrix[j][is],distmatrix[j][js]);
        
        for (j = 0; j < is; j++) distmatrix[is][j] = distmatrix[n-1][j];
        for (j = is+1; j < n-1; j++) distmatrix[j][is] = distmatrix[n-1][j];
        
        /* Update clusterids */
        result[nelements-n].left = clusterid[is];
        result[nelements-n].right = clusterid[js];
        clusterid[js] = n-nelements-1;
        clusterid[is] = clusterid[n-1];
    }
    return result;
}


static
double find_closest_pair(int n, matrix_t &distmatrix, int* ip, int* jp)
/*
 This function searches the distance matrix to find the pair with the shortest
 distance between them. The indices of the pair are returned in ip and jp; the
 distance itself is returned by the function.
 
 n          (input) int
 The number of elements in the distance matrix.
 
 distmatrix (input) double**
 A ragged array containing the distance matrix. The number of columns in each
 row is one less than the row index.
 
 ip         (output) int*
 A pointer to the integer that is to receive the first index of the pair with
 the shortest distance.
 
 jp         (output) int*
 A pointer to the integer that is to receive the second index of the pair with
 the shortest distance.
 */
{ int i, j;
    double temp;
    double distance = distmatrix[1][0];
    *ip = 1;
    *jp = 0;
    for (i = 1; i < n; i++)
    { for (j = 0; j < i; j++)
    { temp = distmatrix[i][j];
        if (temp<distance)
        { distance = temp;
            *ip = i;
            *jp = j;
        }
    }
    }
    return distance;
}






vector<int> cutTree (vector<TreeNode> Tree, int nclusters)

/*
 Purpose
 =======
 
 The cuttree routine takes the output of a hierarchical clustering routine, and
 divides the elements in the tree structure into clusters based on the
 hierarchical clustering result. The number of clusters is specified by the user.
 
 Arguments
 =========
 
 tree           (input) Node[nelements-1]
 The clustering solution. Each node in the array describes one linking event,
 with tree[i].left and tree[i].right representig the elements that were joined.
 The original elements are numbered 0..nelements-1, nodes are numbered
 -1..-(nelements-1).
 
 nclusters      (input) int
 The number of clusters to be formed.
 
 Return
 ======
 clusterid      (output) vector<int>
 The number of the cluster to which each element was assigned. Space for this
 array should be allocated before calling the cuttree routine. If a memory
 error occured, all elements in clusterid are set to -1.
 
 ========================================================================
 */

{
    int i, j, k;
    int nelements = Tree.size();
    //cluster index
    int icluster = 0;
    const int njoin = nelements-nclusters; // number of nodes to join
    vector<int> clusterid(nelements,-1);
    vector<int> nodeid(njoin,-1);
    for (i = nelements-2; i >= njoin; i--) {
        k = Tree[i].left;
        if (k>=0)
        { clusterid[k] = icluster;
            icluster++;
        }
        k = Tree[i].right;
        if (k>=0)
        { clusterid[k] = icluster;
            icluster++;
        }
    }
    for (i = njoin-1; i >= 0; i--) {
        if(nodeid[i]<0) {
            j = icluster;
            nodeid[i] = j;
            icluster++;
        }
        else j = nodeid[i];
        k = Tree[i].left;
        if (k<0)
            nodeid[-k-1] = j;
        else
            clusterid[k] = j;
        k = Tree[i].right;
        if (k<0)
            nodeid[-k-1] = j;
        else
            clusterid[k] = j;
    }
    return clusterid;
}

//Best clustering number: Gap Statistics.
int GapStatistics (matrix_t &distmatrix, vector<TreeNode> Tree, char method)

/*
 Purpose
 =======
 
 The GapStatistics routine takes the output of a hierarchical clustering cutting result, and evaluate the gap statistics based on reference uniform distribution.
 
 Arguments
 =========
 
 clusterid      (input) vector<int>
 The number of elements that were clustered.
 
 tree           (input) Node[nelements-1]
 The clustering solution. Each node in the array describes one linking event,
 with tree[i].left and tree[i].right representig the elements that were joined.
 The original elements are numbered 0..nelements-1, nodes are numbered
 -1..-(nelements-1).
 
 nclusters      (input) int
 The number of clusters to be formed.
 
 clusterid      (output) int[nelements]
 The number of the cluster to which each element was assigned. Space for this
 array should be allocated before calling the cuttree routine. If a memory
 error occured, all elements in clusterid are set to -1.
 
 Result
 ======
 
 
 
 ========================================================================
 */

{
    int bestNclstr=-1;
    float totalVariance = CalWithinGrpDist(distmatrix);
    vector<TreeNode> tn;
    
    switch (method) {
        case 0:
            tn = SingleLinkClstr(distmatrix);
            break;
        case 1:
            tn = AverageLinkClstr(distmatrix);
            break;
        case 2:
            tn = CompleteLinkClstr(distmatrix);
            break;
        default:
            tn = AverageLinkClstr(distmatrix);
            break;
    }
    
    vector<vector<TreeNode> > refTns;
    vector<matrix_t > refDistMatrix;
    for (int i=0; i<GAP_PERMUTATION; i++) {
        refDistMatrix.push_back(DistMatrixOfUniformDispersion(distmatrix.size(), totalVariance, i));
        
        switch (method) {
            case 0:
                refTns.push_back(SingleLinkClstr(refDistMatrix[i]));
                break;
            case 1:
                refTns.push_back(AverageLinkClstr(refDistMatrix[i]));
                break;
            case 2:
                refTns.push_back(CompleteLinkClstr(refDistMatrix[i]));
                break;
            default:
                refTns.push_back(AverageLinkClstr(refDistMatrix[i]));
                break;
        }

    }
    //cout<<totalVariance<<" "<<CalWithinGrpDist(refDistMatrix[0])<<endl;
    //for each cluster number
    vector<float> gap, Sk;
    int min_num = atoi(option["minclstr"].c_str());
    int max_num = std::min(atoi(option["maxclstr"].c_str()), int(distmatrix.size()));
    for (int k=min_num; k<max_num; k++){
        float logWk = CalLogWk(tn, k, distmatrix);
        float SumLogWkB = 0;
        float SumLogWkBSquare = 0;
        for (int i=0; i<GAP_PERMUTATION; i++) {
            float logWkB = CalLogWk(refTns[i], k, refDistMatrix[i]);
            SumLogWkB += logWkB;
            SumLogWkBSquare += logWkB*logWkB;
        }
        gap.push_back(SumLogWkB/GAP_PERMUTATION-logWk);
        //total standard error
        Sk.push_back( sqrtf((SumLogWkBSquare/GAP_PERMUTATION-SumLogWkB*SumLogWkB/GAP_PERMUTATION/GAP_PERMUTATION)*(1+1/float(GAP_PERMUTATION))));
        if (gap.size()>=2 && (gap[gap.size()-2]>=gap[gap.size()-1]-GAP_SD_MULTIPLIER*Sk[gap.size()-1])){
            bestNclstr = gap.size()-2;
            break;
        }
    }
    return bestNclstr==-1?1:bestNclstr+min_num;
}

float CalLogWk(const vector<TreeNode> &tn, int k, const matrix_t &distmatrix){
    vector<int> clusterLabel;
    clusterLabel = cutTree(tn, k);
    //0=>clstr0's submatrix....
    vector<matrix_t > clusterMatrice = SubMatrice(distmatrix, clusterLabel);
    float Wk = 0;
    for (int nClstr=0; nClstr<k; nClstr++) {
        Wk += CalWithinGrpDist(clusterMatrice[nClstr]);
    }
    return logf(Wk);
}


// retrieve index(in terms of sliced matrix) of representative element(with smallest distance)of cluster No.K from distance matrix
int GetRepElement(const matrix_t &distmatrix, vector<int> clusterLabel, int clstrID) {
    if (!IsMatrix(distmatrix)) printAndExit("FATAL: Abnormal distance matrix!");
    vector<matrix_t > clusterMatrice = SubMatrice(distmatrix, clusterLabel);
    matrix_t thisCluster = clusterMatrice[clstrID];
    vector<float> rowSum = GetRowSum(thisCluster);
    vector<float>::iterator it = min_element(rowSum.begin(),rowSum.end());
    int currentKey = it-rowSum.begin();
    // get primary key
    /*
    int counter=0;
    for (int primaryKey=0; primaryKey<clusterLabel.size(); primaryKey++) {
        if (clusterLabel[primaryKey]==clstrID) {
            counter++;
            if (counter==currentKey+1) {
                cout<<currentKey<<" "<<primaryKey<<endl;
                return primaryKey;
            }
        }
    }
    return -1;
    */
    return currentKey;
}

vector<float> GetRowSum(const matrix_t &distmatrix) {
    if (!IsMatrix(distmatrix)) printAndExit("FATAL: Abnormal distance matrix!");
    vector<float> rowSum;
    for (int i=0; i<distmatrix.size(); i++) {
        float tempSum=0;
        for (int j=0; j<distmatrix[0].size(); j++) {
            tempSum += distmatrix[i][j];
        }
        rowSum.push_back(tempSum);
    }
    return rowSum;
}

bool IsMatrix(const matrix_t &distmatrix) {
    if (distmatrix.size()==0) {
        return false;
    }
    if (distmatrix.size()==1) {
        return true;
    }
    for (int i=1; i<distmatrix.size(); i++) {
        if (distmatrix[i].size()!=distmatrix[0].size()) {
            return false;
        }
    }
    return true;
}

//cut matrix by group label

vector<matrix_t > SubMatrice(matrix_t matrix, vector<int> label) {
    //total cluster number
    std::set<int> tempSet(label.begin(),label.end());
    int nCluster = tempSet.size();
    vector<vector<vector<float> > > clusterNo2matrix;
    //may have problem
    clusterNo2matrix.reserve(nCluster);
    for (int i=0; i<nCluster; i++) {
        vector<int> index;
        for (int j=0; j<matrix.size(); j++) {
            if (label[j]==i) {
                index.push_back(j);
            }
        }
        matrix_t submatrix;
        InitMatrix(submatrix, index.size(), index.size(), 0);
        for (int j=0; j<index.size(); j++) {
            for (int k=0; k<index.size(); k++) {
                submatrix[j][k] = matrix[index[j]][index[k]];
            }
        }
        clusterNo2matrix.push_back(submatrix);
    }
    return clusterNo2matrix;
}



// calculate within-group-distance for given distance matrix
// within-group-distance = (1/2n)*within-group-sum-of-squares
float CalWithinGrpDist(matrix_t &distmatrix) {
    float distAll = 0;
    for (int i=0; i<distmatrix.size(); i++) {
        for (int j=0; j<distmatrix[0].size(); j++) {
            distAll += distmatrix[i][j];
        }
    }
    return distAll/2/distmatrix.size();
}



//MC simulation 
matrix_t DistMatrixOfUniformDispersion(int nPoints, float var, int seed) {
// dispersion bound in circle, it can be proved that the expectation of distance between two random points is 4/9R. r=var*9/2/(nPoints-1) ???
    vector<PolarPos> Points;
    float rBound = var*9/4/(nPoints-1);
    //if no seed, all matrix are same
    srand((unsigned)(time(NULL)+seed));
    for (int i=0; i<nPoints; i++) {
        PolarPos tempPoint = {rand()/(float)(RAND_MAX/(2*PI)),rand()/(float)(RAND_MAX/rBound)};
        Points.push_back(tempPoint);
    }
    matrix_t distmatrix;
    InitMatrix(distmatrix, nPoints, nPoints, 0);
    for (int i=0; i<nPoints; i++) {
        for (int j=0; j<=i; j++) {
            float distSquare = Points[i].length*Points[i].length + Points[j].length*Points[j].length - 2*Points[i].length*Points[j].length*cosf(Points[i].angle-Points[j].angle);
            //if (distSquare<0) cout<<Points[i].length<<Points[j].length<<endl;
            distmatrix[i][j] = sqrtf(distSquare);
            distmatrix[j][i] = distmatrix[i][j];
        }
    }
    return distmatrix;
}




