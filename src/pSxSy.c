#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>
#include "llik.h"

/******************************************************************************/
/*                                                                            */
/*                          MULTIPLE RELATED PAIRS                            */
/*                                P(Sx, Sy)                                   */
/*                     mixed radix incrementing approach                      */
/*                                                                            */ 
/******************************************************************************/

/* ixy: 0-based sorted indices - which ux are in uy (or uxy); likewise for iyx */
/* mmax = min(|Sxy|, nm) [mmax <- min(sum(vxy), nm)], r[i] = 0 for i > nm     */
/* sprob is on log scale, so only non-zero ones, can't go beyond mmax */

///***!!! DO NOT CHANGE ORIGINAL VARIABLES, they are NOT copied!!!

double probSxSy(double *rm, double *sprob, int nm, int mmax) 
{
  int i, ikeep, m, nmall;
  double prob, sum;
  int ibd[nm], ibdmax[nm];
  double logr[nm], log1r[nm];

  /* remove r = 0 and r = 1, update nm and m accordingly */
  nmall = nm;
  m = 0;  // sum(ibd)
  ikeep = 0;
  prob  = 0; 
  for (i = 0; i < nmall; i++) {
    //    if (fabs(rm[i]) < EPS) {             // rm[i] = 0
    if (rm[i] == 0) {         // same representation (?)
      nm--;
      //    } else if (fabs(rm[i] - 1) < EPS) {  // rm[i] = 1
    } else if (rm[i] == 1) {  // same representation (?)
      m++;
      nm--;
    } else {
      logr[ ikeep] = log(rm[i]);
      log1r[ikeep] = log(1 - rm[i]);
      ibd[  ikeep] = 0;
      prob        += log1r[ikeep];
      ikeep++;
    }
  }
  if (mmax < m) {  // number of r = 1 is greater than number of shared alleles
    return 0;
  }

  /* in line with unified algorithm */
  for (i = 0; i < nm; i++) {
    ibdmax[i] = 1;
  }
  for (i = mmax - m; i < nm; i++) {  /* if mmax < nm */ //???
    ibdmax[i] = 0;
  }

  sum = exp(prob + sprob[m]);  // # of r = 1 is greater than # of shared alleles
  i = nm - 1; 
  while (!equalArr(ibd, ibdmax, nm)) {  //equalArr() returns 1 if nm = 0
    if (ibd[i] == 1 || m == mmax) {
      if (ibd[i] == 0) {
	i--;
	continue; 
      }
      m -= ibd[i]; // can be 0!!!      
      ibd[i] = 0;
      prob += log1r[i] - logr[i]; 
      i--;
    } else {
      ibd[i] = 1;
      m++;
      prob += logr[i] - log1r[i];
      sum += exp(prob + sprob[m]);  // that's where we'll add 0's and 1's
      i = nm - 1;
    }
  } 
  return sum;
}
 
double probSxSyEqr(double *rm, double *sprob, int nm, int mmax, double *cbinom) 
{
  int m;
  double sum, logr, log1r;

  /* edge cases */
  if (*rm == 0) {  // or rm[0]  // same representation (?) 
    //  if (fabs(*rm - 0) < EPS) {  // *rm = 0
    return exp(sprob[0]);
  }
  if (*rm == 1) {               // same representation (?)
    //  if (fabs(*rm - 1) < EPS) {  // *rm = 1
    if (mmax < nm) {
      return 0;
    } else {
      return exp(sprob[mmax]);  
    }
  }

  /* 0 < r < 1 */
  logr  = log(*rm);
  log1r = log(1 - *rm);
  sum = 0;
  for (m = 0; m <= mmax; m++) {
    sum += exp(cbinom[m] + m*logr + (nm - m)*log1r + sprob[m]);
  }
  return sum;
}

void probSxSyCond(int *vx, int *vy, double *logpxy, double *logj, double *factj,
		  int nx, int ny, int nux, int nuy, int nuxy, int *ixy, int *iyx,
		  double combx, double comby, double *sprob, int *mmax, int nm)
{
  int vxy[nuxy], s[nuxy], vxc[nuxy], vyc[nuxy], vmax[nuxy];
  int i, numx, numy, nums, nsi, mm;
  double prob; 

  prob = combx + comby;
  numx = nx; // do we need nx, ny (can start with numx, numy) // numerator
  numy = ny;
  nums = 0;

  /* find vxy (for Sxy), initiate s with 0's, calculate mmax, vmax */
  *mmax = 0;
  for (i = 0; i < nuxy; i++) {
    s[i]   = 0;                       /* first subset - empty */
    vxc[i] = vx[ixy[i]];
    vyc[i] = vy[iyx[i]];
    vxy[i] = MIN(vxc[i], vyc[i]);
    *mmax += vxy[i];
  }
  *mmax = MIN(*mmax, nm);  
  vMax(vxy, nuxy, *mmax, vmax);  /* vmax equal vxy if *mmax = nm */
  mm = *mmax;  // almost no gain  

  /* initialize sprob with 0's, fill in the independence case */
  sprob[0] = exp(prob);    /* independence, Sxy is an empty set */  
  for (i = 1; i <= mm; i++) sprob[i] = 0; 

  /* go through the subsets */
  i = nuxy - 1;

  // !!! note factj starts with 0 (but not logj) !!!
  while (!equalArr(s, vmax, nuxy)) {
    if (s[i] == vxy[i] || nums == mm) {
      nsi = s[i]; 
      if (nsi == 0) {  // small gain in efficiency
	i--;
	continue;
      }
      prob -= factj[nums] + factj[numx]   + factj[numy]  // subtract old
            - factj[nsi]  - factj[vxc[i]] - factj[vyc[i]];
      s[i] = 0;
      vxc[i] += nsi;
      vyc[i] += nsi;
      nums -= nsi;
      numx += nsi;
      numy += nsi;
      prob += factj[nums] + factj[numx]   + factj[numy]  // add new
                          - factj[vxc[i]] - factj[vyc[i]]
            + logpxy[i]*nsi;
      i--;
    } else {
      s[i]++;
      vxc[i]--;
      vyc[i]--;
      nums++;
      numx--;
      numy--; 
      prob += logj[nums - 1] - logj[numx]   - logj[numy]
	    - logj[s[i] - 1] + logj[vxc[i]] + logj[vyc[i]]
	    - logpxy[i]; 
      sprob[nums] += exp(prob);    
      i = nuxy - 1;
    }
  }
}

void vMax(int *base, int nbase, int sumlim, int *vmax)
{
  int i, rem;
  for (i = 0; i < nbase; i++) {
    vmax[i] = 0;
  }
  rem = sumlim;
  i = 0;
  while (rem > 0) {
    vmax[i] = MIN(base[i], rem);
    rem -= base[i];
    i++;
  }
}

int equalArr(int *arr1, int *arr2, int narr)
{
  int i;
  for (i = 0; i < narr; i++) {
    if (arr1[i] != arr2[i]) {
      return 0;
    }
  }
  return 1;
}


  
