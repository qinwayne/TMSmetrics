/* 
   Matlab call: H = sweep_kl (X,m_min,m_max,TAU);
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include "mex.h"
#include "matrix.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int pp, m, m1, m2;
	double *D_2, *X, *b, *DV1, *DV2, *H, *TAU;
	int p, k, i, j, n, ind1, ind2, nn, tau, n_emb, n_tau, count;
	double aux1, aux2, h;

	/* ************************* */
	/* Check number of arguments */
	/* ************************* */
	if(nrhs != 4)
		mexErrMsgTxt("Wrong number of inputs!");
	else if(nlhs>1)
		mexErrMsgTxt("Too many output arguments!");


	/* ***************** */
	/* Initialise Inputs */
	/* ***************** */
	
	ind1 = -1;
	X = mxGetPr(prhs[++ind1]);
	pp = (int) mxGetN(prhs[ind1]);
	k = (int) mxGetM(prhs[ind1]);
	if (pp==1)
		pp = k;

	b = mxGetPr(prhs[++ind1]);
	m1 = (int) *b;

	b = mxGetPr(prhs[++ind1]);
	m2 = (int) *b;

	b = mxGetPr(prhs[++ind1]);
	TAU = b;
	n_tau = (int) mxGetN(prhs[ind1]);
	k = (int) mxGetM(prhs[ind1]);
	if (n_tau==1)
		n_tau = k;

	nn = pp-(m2-1)*TAU[n_tau-1];
	n_emb = m2-m1+1;

	/* ****************** */
	/* Initialise Outputs */
	/* ****************** */
	ind1 = -1;
	plhs[++ind1] = mxCreateDoubleMatrix(n_emb,n_tau, mxREAL);
	H = (double*) mxGetPr(plhs[ind1]);


	/* ******************** */
	/* Array Initialisation */
	/* ******************** */
	D_2 = (double*) calloc (nn*nn, sizeof(double));


	/* ********* */
	/* Here Goes */
	/* ********* */
	for (n=count=0; n<n_tau; n++)
	  {
	  tau = (int) TAU[n];

			/* Up to m1 */
	  /* Compute Distances */
	  for (i=0; i<nn; i++)
	    {
	    DV1 = &X[i];
	    D_2[i*nn+i] = 100000.0;
	    for (j=i+1; j<nn; j++)
		{
		DV2 = &X[j];
		for (k=aux1=ind1=0; k<m1; k++,ind1+=tau)
			aux1 += (aux2=DV1[ind1]-DV2[ind1])*aux2;
		D_2[i*nn+j] = aux1;
		D_2[j*nn+i] = aux1;
		}
	    }
	  
	  /* Compute Nearest Neighbour Distances */
	  for (i=ind1=H[count]=aux1=0; i<nn; i++)
	    {
	    for (j=0,ind2=ind1; j<nn; j++,ind1++)
		if ((D_2[ind1]<D_2[ind2]) && (D_2[ind1]!=0))
			ind2 = ind1;
	    h = (double) nn*sqrt(D_2[ind2]);
	    H[count] += log(h+(h==0));
	    aux1 += (double) (h==0);
	    }

	  /* Compute Kozachenko-Leonenko Estimate */
	  H[count] /= (double) nn-aux1;
	  H[count] += 1.270347180559945;

	  count ++;


			/* (m1+1) -> m2 */
	  for (m=m1+1; m<=m2; m++,count++)
	    {
	    /* Incrementally Compute Distances */
	    for (i=0; i<nn; i++)
	      {
	      DV1 = &X[i];
	      for (j=i+1; j<nn; j++)
		{
		DV2 = &X[j];

		ind1 = (m-1)*tau;
		D_2[i*nn+j] += (aux2=DV1[ind1]-DV2[ind1])*aux2;
		D_2[j*nn+i] += (aux2=DV1[ind1]-DV2[ind1])*aux2;
		}
	      }

	    /* Compute Nearest Neighbour Distances */
	    for (i=ind1=aux1=H[count]=0; i<nn; i++)
	      {
	      for (j=0,ind2=ind1; j<nn; j++,ind1++)
		if ((D_2[ind1]<D_2[ind2]) && (D_2[ind1]!=0))
			ind2 = ind1;
		h = (double) nn*sqrt(D_2[ind2]);
		H[count] += log(h+(h==0));
		aux1 += (double) (h==0);
	      }

	    /* Compute Kozachenko-Leonenko Estimate */
	    H[count] /= (double) nn-aux1;
	    H[count] += 1.270347180559945;
	    }
	  }
	

	free (D_2);

}



