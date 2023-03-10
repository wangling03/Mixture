//------------------------------------------------------------------------------
// The code was written by Vlad I. Morariu 
// and is copyrighted under the Lesser GPL: 
//
// Copyright (C) 2007 Vlad I. Morariu
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 or later.
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
// See the GNU Lesser General Public License for more details. 
// You should have received a copy of the GNU Lesser General Public
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, 
// MA 02111-1307, USA.  
//
// The author may be contacted via email at: morariu(at)cs(.)umd(.)edu 
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// File    : sample.cpp
// Purpose : Show example of how figtree() function can be used to evaluate
//           gauss transforms.  
// Author  : Vlad I. Morariu       morariu(at)cs(.)umd(.)edu
// Date    : 2007-06-25
// Modified: 2008-01-23 to add comments and make the sample clearer
//------------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h> // for memset
#include <math.h>   // for abs
#include "figtree.h"

#ifdef WIN32
  // Link to figtree.dll if WIN32.  If not, then we assume that 
  // libfigtree.so or libfigtree.a are linked to from the Makefile.
  // The locations of figtree.dll and ANN.dll must either be in the PATH
  // environment variable, or in the same location from which sample.exe is
  // is executed.
  #pragma comment(lib,"../lib/figtree.lib")
#endif

// This program runs all figtree methods with different parameters and tests
// whether the maximum absolute error is below epsilon.

#ifndef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#endif

int main()
{
  int ds[] = {1,2,3,4,5,6};
  double hs[] = {.01,.02,.04,.08,.16,.32,.64};
  int nds = sizeof(ds)/sizeof(ds[0]);
  int nhs = sizeof(hs)/sizeof(hs[0]);

  // The number of targets (vectors at which gauss transform is evaluated).
  int M = 1000;

  // The number of sources which will be used for the gauss transform. 
  int N = 1000;

  // Number of weights.  For each set of weights a different Gauss Transform is computed, 
  // but by giving multiple sets of weights at once some overhead can be shared.
  int W = 2;

  // Desired maximum absolute error after normalizing output by sum of weights.
  // If the weights, q_i (see below), add up to 1, then this is will be the 
  // maximum absolute error.
  // The smaller epsilon is, the more accurate the results will be, at the
  // expense of increased computational complexity.
  double epsilon = 1e-2;

  // try each parameter combination ntrial times
  int ntrials = 3;

  for( int dind = 0; dind < nds; dind++ )
  {
    for( int hind = 0; hind < nhs; hind++ )
    {
      for( int tind = 0; tind < ntrials; tind++ )
      {
        // The dimensionality of each sample vector.
        int d = ds[dind];

        // The bandwidth.  NOTE: this is not the same as standard deviation since 
        // the Gauss Transform sums terms exp( -||x_i - y_j||^2 / h^2 ) as opposed
        // to  exp( -||x_i - y_j||^2 / (2*sigma^2) ).  Thus, if sigma is known, 
        // bandwidth can be set to h = sqrt(2)*sigma.
        double h = hs[hind];
    
        // The source array.  It is a contiguous array, where
        // ( x[i*d], x[i*d+1], ..., x[i*d+d-1] ) is the ith d-dimensional sample.
        double * x = new double[d*N];
        // initialize with random values between 0 and 1
        for( int i = 0; i < d*N; i++ )
          x[i] = rand()/(double(RAND_MAX)+1);

        // The target array.  It is a contiguous array, where
        // ( y[j*d], y[j*d+1], ..., y[j*d+d-1] ) is the jth d-dimensional sample.
        double * y = new double[d*M];
        // initialize with random values between 0 and 1
        for( int i = 0; i < d*M; i++ )
          y[i] = rand()/(double(RAND_MAX)+1);

        // The weight array.  The ith weight is associated with the ith source sample.
        // To evaluate the Gauss Transform with the same sources and targets, but 
        // different sets of weights, add another row of weights and set W = 2.  
        double * q = new double[W*N];
        double * Q = new double[W]; // sum of each set of qs (error is relative to Q)
        // initialize with random values between 0 and 1
        for( int i = 0; i < W; i++ )
        {
          Q[i] = 0;
          for( int j = 0; j < N; j++ )
          {
            q[i*N+j] = rand()/(double(RAND_MAX)+1);
            Q[i] += q[i*N+j];
          }
        }

        // allocate array into which the result of the Gauss Transform will be stored for each
        // target sample.  The first M elements will correspond to the Gauss Transform computed
        // with the first set of weights, second M elements will correspond to the G.T. computed
        // with the second set of weights, etc.
        double * g_sf = new double[W*M];
        double * g_sf_tree = new double[W*M];
        double * g_ifgt_u = new double[W*M];
        double * g_ifgt_tree_u = new double[W*M];
        double * g_ifgt_nu = new double[W*M];
        double * g_ifgt_tree_nu = new double[W*M];
    
        // initialize all output arrays to zero
        memset( g_sf         , 0, sizeof(double)*W*M );
        memset( g_sf_tree     , 0, sizeof(double)*W*M );
        memset( g_ifgt_u     , 0, sizeof(double)*W*M );
        memset( g_ifgt_tree_u , 0, sizeof(double)*W*M );
        memset( g_ifgt_nu    , 0, sizeof(double)*W*M );
        memset( g_ifgt_tree_nu, 0, sizeof(double)*W*M );

        // evaluate gauss transform using direct (slow) method
        figtree( d, N, M, W, x, h, q, y, epsilon, g_sf, FIGTREE_EVAL_DIRECT );

        // evaluate gauss transform using direct method with approximate nearest neighbors
        figtree( d, N, M, W, x, h, q, y, epsilon, g_sf_tree, FIGTREE_EVAL_DIRECT_TREE );

        // evaluate gauss transform using FIGTREE (truncated series), estimating parameters with and without
        //   the assumption that sources are uniformly distributed
        figtree( d, N, M, W, x, h, q, y, epsilon, g_ifgt_u, FIGTREE_EVAL_IFGT, FIGTREE_PARAM_UNIFORM, 1 );
        figtree( d, N, M, W, x, h, q, y, epsilon, g_ifgt_nu, FIGTREE_EVAL_IFGT, FIGTREE_PARAM_NON_UNIFORM, 1 );

        // evaluate gauss transform using FIGTREE (truncated series), estimating parameters with and without
        //   the assumption that sources are uniformly distributed
        figtree( d, N, M, W, x, h, q, y, epsilon, g_ifgt_tree_u, FIGTREE_EVAL_IFGT_TREE, FIGTREE_PARAM_UNIFORM, 1 );
        figtree( d, N, M, W, x, h, q, y, epsilon, g_ifgt_tree_nu, FIGTREE_EVAL_IFGT_TREE, FIGTREE_PARAM_NON_UNIFORM, 1 );

        // compute absolute error of the Gauss Transform at each target and for all sets of weights.
        double err_sf_tree      = fabs(g_sf_tree[0]      - g_sf[0])/Q[0];
        double err_ifgt_u       = fabs(g_ifgt_u[0]       - g_sf[0])/Q[0];
        double err_ifgt_nu      = fabs(g_ifgt_nu[0]      - g_sf[0])/Q[0];
        double err_ifgt_tree_u  = fabs(g_ifgt_tree_u[0]  - g_sf[0])/Q[0];
        double err_ifgt_tree_nu = fabs(g_ifgt_tree_nu[0] - g_sf[0])/Q[0];

        for( int i = 0; i < W; i++)
          for( int j = 0; j < M; j++ )
          {
            err_sf_tree      = MAX(err_sf_tree,      fabs(g_sf_tree[i*M+j]     -g_sf[i*M+j])/Q[i]);
            err_ifgt_u       = MAX(err_ifgt_u,       fabs(g_ifgt_u[i*M+j]      -g_sf[i*M+j])/Q[i]);
            err_ifgt_nu      = MAX(err_ifgt_nu,      fabs(g_ifgt_nu[i*M+j]     -g_sf[i*M+j])/Q[i]);
            err_ifgt_tree_u  = MAX(err_ifgt_tree_u,  fabs(g_ifgt_tree_u[i*M+j] -g_sf[i*M+j])/Q[i]);
            err_ifgt_tree_nu = MAX(err_ifgt_tree_nu, fabs(g_ifgt_tree_nu[i*M+j]-g_sf[i*M+j])/Q[i]);
          }

        // deallocate memory
        delete [] g_sf;
        delete [] g_sf_tree;
        delete [] g_ifgt_u;
        delete [] g_ifgt_nu;
        delete [] g_ifgt_tree_u;
        delete [] g_ifgt_tree_nu;

        // print out results for all six ways to evaluate
        printf("\n\n\n");
        printf("Errors( d=%i, h=%e ):\n  sf-tree: %6.4e\n  ifgt-u: %6.4e\n  ifgt-nu: %6.4e\n  ifgt_tree_u: %6.4e\n  ifgt_tree_nu: %6.4e\n", 
                 d, h, err_sf_tree, err_ifgt_u, err_ifgt_nu, err_ifgt_tree_u, err_ifgt_tree_nu );
        if( (err_sf_tree < epsilon) &&
            (err_ifgt_u < epsilon) && 
            (err_ifgt_nu < epsilon) && 
            (err_ifgt_tree_u < epsilon) &&
            (err_ifgt_tree_nu < epsilon) )
        {
          printf("\n\n\nTest Passed.\n\n\n\n");
        }
        else
        {
          printf("\n\n\nTEST FAILED!\n\n\n\n");
          return -1;
        }
      }
    }
  }
   
  return 0;
}
