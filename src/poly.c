/*
*				poly.c
*
* Manage polynomials.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	AstrOmatic WCS library
*
*	Copyright:		(C) 1998-2011 Emmanuel Bertin -- IAP/CNRS/UPMC
*
*	License:		GNU General Public License
*
*	AstrOmatic software is free software: you can redistribute it and/or
*	modify it under the terms of the GNU General Public License as
*	published by the Free Software Foundation, either version 3 of the
*	License, or (at your option) any later version.
*	AstrOmatic software is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*	You should have received a copy of the GNU General Public License
*	along with AstrOmatic software.
*	If not, see <http://www.gnu.org/licenses/>.
*
*	Last modified:		20/12/2011
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include	<math.h>
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>

#include	"poly.h"

#define	QCALLOC(ptr, typ, nel) \
		{if (!(ptr = (typ *)calloc((size_t)(nel),sizeof(typ)))) \
		  qerror("Not enough memory for ", \
			#ptr " (" #nel " elements) !");;}

#define	QMALLOC(ptr, typ, nel) \
		{if (!(ptr = (typ *)malloc((size_t)(nel)*sizeof(typ)))) \
		  qerror("Not enough memory for ", \
			#ptr " (" #nel " elements) !");;}

/********************************* qerror ************************************/
/*
I hope it will never be used!
*/
void	qerror(char *msg1, char *msg2)
  {
  fprintf(stderr, "\n> %s%s\n\n",msg1,msg2);
  exit(-1);
  }


/****** poly_init ************************************************************
PROTO   polystruct *poly_init(long *group, long ndim, long *degree, long ngroup)
PURPOSE Allocate and initialize a polynom structure.
INPUT   1D array containing the group for each parameter,
        number of dimensions (parameters),
        1D array with the polynomial degree for each group,
        number of groups.
OUTPUT  polystruct pointer.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 30/08/2011
 ***/
polystruct	*poly_init(long *group, long ndim, long *degree, long ngroup)
  {
   void	qerror(char *msg1, char *msg2);
   polystruct	*poly;
   char		str[512];
   long		nd[POLY_MAXDIM];
   long		*groupt,
		d,g,n, num,den, dmax;

  QCALLOC(poly, polystruct, 1);
  if ((poly->ndim=ndim) > POLY_MAXDIM)
    {
    sprintf(str, "The dimensionality of the polynom (%ld) exceeds the maximum\n"
		"allowed one (%d)", ndim, POLY_MAXDIM);
    qerror("*Error*: ", str);
    }

  if (ndim)
    QMALLOC(poly->group, long, poly->ndim);
    for (groupt=poly->group, d=ndim; d--;)
      *(groupt++) = *(group++)-1;

  poly->ngroup = ngroup;
  if (ngroup)
    {
    group = poly->group;	/* Forget the original *group */

    QMALLOC(poly->degree, long, poly->ngroup);

/*-- Compute the number of context parameters for each group */
    memset(nd, 0, ngroup*sizeof(long));
    for (d=0; d<ndim; d++)
      {
      if ((g=group[d])>=ngroup)
        qerror("*Error*: polynomial GROUP out of range", "");
      nd[g]++;
      }
    }

/* Compute the total number of coefficients */
  poly->ncoeff = 1;
  for (g=0; g<ngroup; g++)
    {
    if ((dmax=poly->degree[g]=*(degree++))>POLY_MAXDEGREE)
      {
      sprintf(str, "The degree of the polynom (%ld) exceeds the maximum\n"
		"allowed one (%d)", poly->degree[g], POLY_MAXDEGREE);
      qerror("*Error*: ", str);
      }

/*-- There are (n+d)!/(n!d!) coeffs per group = Prod_(i<=d)(n+i)/Prod_(i<=d)i */
    n = nd[g];
    d = dmax>n? n: dmax;
    for (num=den=1; d; num*=(n+dmax--), den*=d--);
    poly->ncoeff *= num/den;
    }

  QMALLOC(poly->basis, double, poly->ncoeff);
  QCALLOC(poly->coeff, double, poly->ncoeff);

  return poly;
  }


/****** poly_end *************************************************************
PROTO   void poly_end(polystruct *poly)
PURPOSE Free a polynom structure and everything it contains.
INPUT   polystruct pointer.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP, Leiden observatory & ESO)
VERSION 09/04/2000
 ***/
void	poly_end(polystruct *poly)
  {
  if (poly)
    {
    free(poly->coeff);
    free(poly->basis);
    free(poly->degree);
    free(poly->group);
    free(poly);
    }
  }


/****** poly_func ************************************************************
PROTO   double poly_func(polystruct *poly, double *pos)
PURPOSE Evaluate a multidimensional polynom.
INPUT   polystruct pointer,
        pointer to the 1D array of input vector data.
OUTPUT  Polynom value.
NOTES   Values of the basis functions are updated in poly->basis.
AUTHOR  E. Bertin (IAP)
VERSION 03/03/2004
 ***/
double	poly_func(polystruct *poly, double *pos)
  {
   double	xpol[POLY_MAXDIM+1];
   double      	*post, *xpolt, *basis, *coeff, xval;
   long double	val;
   long		expo[POLY_MAXDIM+1], gexpo[POLY_MAXDIM+1];
   long	       	*expot, *degree,*degreet, *group,*groupt, *gexpot,
			d,g,t, ndim;

/* Prepare the vectors and counters */
  ndim = poly->ndim;
  basis = poly->basis;
  coeff = poly->coeff;
  group = poly->group;
  degree = poly->degree;
  if (ndim)
    {
    for (xpolt=xpol, expot=expo, post=pos, d=ndim; --d;)
      {
      *(++xpolt) = 1.0;
      *(++expot) = 0;
      }
    for (gexpot=gexpo, degreet=degree, g=poly->ngroup; g--;)
      *(gexpot++) = *(degreet++);
    if (gexpo[*group])
      gexpo[*group]--;
    }

/* The constant term is handled separately */
  val = *(coeff++);
  *(basis++) = 1.0;
  *expo = 1;
  *xpol = *pos;

/* Compute the rest of the polynom */
  for (t=poly->ncoeff; --t; )
    {
/*-- xpol[0] contains the current product of the x^n's */
    val += (*(basis++)=*xpol)**(coeff++);
/*-- A complex recursion between terms of the polynom speeds up computations */
/*-- Not too good for roundoff errors (prefer Horner's), but much easier for */
/*-- multivariate polynomials: this is why we use a long double accumulator */
    post = pos;
    groupt = group;
    expot = expo;
    xpolt = xpol;
    for (d=0; d<ndim; d++, groupt++)
      if (gexpo[*groupt]--)
        {
        ++*(expot++);
        xval = (*(xpolt--) *= *post);
        while (d--)
          *(xpolt--) = xval;
        break;
        }
      else
        {
        gexpo[*groupt] = *expot;
        *(expot++) = 0;
        *(xpolt++) = 1.0;
        post++;
        }
    }

  return (double)val;
  }

