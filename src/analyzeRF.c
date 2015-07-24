/******************************************************************************
    Copyright (C) 2015 Arthur Vigil

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
******************************************************************************/

#include <R.h>
#include <Rmath.h>

/* Node status */
#define NODE_TERMINAL -1
#define NODE_TOSPLIT  -2
#define NODE_INTERIOR -3

void predictTree(double *x, int n, int mdim, int *treemap,
                      int *nodestatus, double *xbestsplit,
                      int *bestvar, int *nodeclass,
                      int treeSize, int *cat, int nclass,
                      int *jts, int *nodex,
                      int *frequency, int *abundant, int *deficient);

void _traceTree(int *treemap, double *x, int *bestvar, double *bestsplit,
               int *nodepred, int *endnode, int *prediction,
               int *frequency, int *abundance, int *deficiency){
    int k=0, v=0;
    double t= 0.;
    while(treemap[(k*3)+2] != -1){
        //Rprintf("k = %d %d %d %d %d\n", k+1, treemap[(k*3)],treemap[(k*3)+1],treemap[(k*3)+2], bestvar[k]);
        v = bestvar[k]-1;  // variable used by kth node
        t = bestsplit[k];   // threshold used at kth node
        frequency[v]++;
        if(x[v]<=t){
            k = treemap[(k*3)]-1;
            deficiency[v]++;
            //Rprintf("\t%f <= %f, -> %d\n", x[v], t, k+1);
        } else{
            k = treemap[(k*3)+1]-1;
            abundance[v]++;
            //Rprintf("\t%f > %f -> %d\n", x[v], t, k+1);
        }
    }
}

/*
parameters:
 mdim - number of features
 ntest - number of observations
 nclass - number of classes
 maxcat - max categories for a feature
 nrnodes - max number of nodes in a given tree
 ntree - number of trees in the forest
 x - data
 xbestsplit - variable used to split at a node
 pid -
 cutoff - voting threshold for choosing a class
 countts - accumulated votes for each class
 treemap - stores structure of trees in forest
 nodestatus - status of each node (eg terminal?) in tree
 cat - number of categories for each feature
 nodeclass - prediction of each (terminal) node in a tree
 jts - predictions for each tree
 jet - predictions for each observation
 bestvar - variable used to split at node in tree
 node -
 treesize - number of nodes in a tree
 keeppred - keep predictions for each tree?
 prox - keep proximities?
 proxMat - proximity matrix
 nodes - should terminal node indicators be returned
 */
void _predict(int *mdim, int *ntest, int *nclass, int *maxcat,
                 int *nrnodes, int *ntree, double *x, int *response, double *xbestsplit,
                 double *pid, double *cutoff, double *countts, int *treemap,
                 int *nodestatus, int *cat, int *nodeclass, int *jts,
                 int *jet, int *bestvar, int *node, int *treeSize,
                 int *keepPred, int *prox, double *proxMat, int *nodes,
                 int *frequency, int *abundant, int *deficient) {
    int i, j, k, n, idxNodes, offset1, offset2, *junk, ntie, decision;
    double crit, cmax;

    memset(countts, 0, (*nclass * *ntest) * sizeof(double));
    idxNodes = 0;
    offset1 = 0;
    offset2 = 0;
    junk = NULL;

    // track variable frequency for each tree
    // a[tree][observation][feature]
    int *local_frequency = malloc((*ntest) * (*mdim) * sizeof(int));
    int *local_abundant = malloc((*ntest) * (*mdim) * sizeof(int));
    int *local_deficient = malloc((*ntest) * (*mdim) * sizeof(int));

    for (j = 0; j < *ntree; ++j) {
        memset(local_frequency, 0, (*ntest) * (*mdim) * sizeof(int));
        memset(local_abundant, 0, (*ntest) * (*mdim) * sizeof(int));
        memset(local_deficient, 0, (*ntest) * (*mdim) * sizeof(int));

        /* predict by the j-th tree */
        predictTree(x, *ntest, *mdim, treemap + 2*idxNodes,
                         nodestatus + idxNodes, xbestsplit + idxNodes,
                         bestvar + idxNodes, nodeclass + idxNodes,
                         treeSize[j], cat, *nclass,
                         jts + offset1, node + offset2,
                         local_frequency, // f[j]
                         local_abundant,  // a[j]
                         local_deficient);   // d[j]

        /* accumulate votes: */
        for (n = 0; n < *ntest; ++n) {
            countts[jts[n + offset1] - 1 + n * *nclass] += 1.0;
        }

        idxNodes += *nrnodes;
        if (*keepPred) offset1 += *ntest;
        if (*nodes)    offset2 += *ntest;

        for(i = 0; i<*ntest;i++){
            // decision of jth tree on ith observation (using base 0 indexing)
            decision = jts[*ntest*j + i] - 1;
            // tree j incorrectly classifies case i
            if(decision != (response[i]-1)){
                continue;
            }
            else{
                for (k=0; k<*mdim;k++){
                    frequency[decision * *mdim + k] += local_frequency[(i * *mdim) + k];
                    abundant[decision * *mdim + k] += local_abundant[(i * *mdim) + k];
                    deficient[decision * *mdim + k] += local_deficient[(i * *mdim) + k];
                }
            }
        }
    }

    /* Aggregated prediction is the class with the maximum votes/cutoff */
    for (n = 0; n < *ntest; ++n) {
        cmax = 0.0;
        ntie = 1;
        for (j = 0; j < *nclass; ++j) {
            crit = (countts[j + n * *nclass] / *ntree) / cutoff[j];
            if (crit > cmax) {
                jet[n] = j + 1;
                cmax = crit;
                ntie = 1;
            }
            /* Break ties at random: */
            if (crit == cmax) {
                if (unif_rand() < 1.0 / ntie) jet[n] = j + 1;
                ntie++;
            }
        }
    }

    free(local_frequency);
    free(local_deficient);
    free(local_abundant);
}

/*
parameters:
 x - array of data
 n - number of data observations
 mdim - number features for each observation
 treemap - map of current tree
 nodestatus - status of each node in tree
 xbestsplit - split threshold for each node in tree
 bestvar - variable used to split for each node in tree
 nodeclass - class vote of node (only makes sense for terminal nodes)
 treeSize - number of nodes in tree
 cat - ?
 nclass - number of classes
 jts - tree's prediction for each observation
 nodex - terminating node
 */
void predictTree(double *x, int n, int mdim, int *treemap,
                      int *nodestatus, double *xbestsplit,
                      int *bestvar, int *nodeclass,
                      int treeSize, int *cat, int nclass,
                      int *jts, int *nodex,
                      int *frequency, int *abundant, int *deficient) {
    int m, i, j, k;

    for (i = 0; i < n; ++i) {
        // make prediction for ith case
        k = 0;
        while (nodestatus[k] != NODE_TERMINAL) {
            m = bestvar[k] - 1;
            frequency[ mdim*i + m ]++;     // frequency for nth case, mth feature
            /* Split by a numerical predictor */
            if (x[m + i * mdim] <= xbestsplit[k]) {
                k = treemap[k * 2] - 1;
                deficient[ mdim*i + m ]++; // deficiency for nth case, mth feature
            } else {
                k = treemap[1 + k * 2] - 1;
                abundant[ mdim*i + m ]++;  // abundance for nth case, mth feature
            }
        }
        /* Terminal node: assign class label */
        jts[i] = nodeclass[k];
        nodex[i] = k + 1;
    }
}
