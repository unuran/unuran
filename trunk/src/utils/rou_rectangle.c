/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      rou_rectangle.c                                              *
 *                                                                           *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      The bounding rectangle for the RoU-methods is computed numerically   *
 *                                                                           *
 *****************************************************************************
 *   (c) 2000 Wolfgang Hoermann and Josef Leydold                            *
 *   Dept. for Statistics, University of Economics, Vienna, Austria          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                         *
 *   59 Temple Place, Suite 330, Boston, MA 02111-1307, USA                  *
 *                                                                           *
 *****************************************************************************
 *****************************************************************************
 *                                                                           *
 *   REFERENCES:                                                             *
 *   [1] Wakefield J.C., Gelfand A.E., Smith A.F.M.                          *
 *       Efficient generation of random variates via the ratio-of-uniforms   *
 *       method.                                                             *
 *       Statistics and Computing (1991) 1, pp (129-133)                     *
 *                                                                           *
 *   [2] Hoermann, W., Leydold J., and Derflinger, G. (2004):                *
 *       Automatic non-uniform random variate generation, Springer, Berlin.  *
 *       Section 2.4, Algorithm 2.9 (RoU), p.35                              *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <unur_source.h>
#include <distr/distr_source.h>
#include <distr/cvec.h>
#include <utils/fmax_source.h>
#include <utils/hooke_source.h>
#include <utils/matrix_source.h>
#include <utils/unur_fp_source.h>
#include <utils/rou_rectangle_source.h>
#include <uniform/urng.h>

/*---------------------------------------------------------------------------*/
/* Constants                                                                 */

/* Convergence parameters for the hooke optimization algorithm */

#define ROU_HOOKE_RHO     (0.5)
#define ROU_HOOKE_EPSILON (1.e-7)
#define ROU_HOOKE_MAXITER (10000)

/* Scaling factor for the computed minimum bounding rectangle.               */
/* The computed rectangle  (0, vmax)x(umin[d], umax[d]) is scaled by this    */
/* factor, i.e. :                                                            */
/* vmax = vmax * ( 1+ ROU_RECT_SCALING)                                     */
/* umin[i] = umin[i] - (umax[i]-umin[i])*ROU_RECT_SCALING/2.                */
/* umax[i] = umax[i] + (umax[i]-umin[i])*ROU_RECT_SCALING/2.                */
#define ROU_RECT_SCALING (1.e-4)


static double _unur_rou_rectangle_aux_vmax(double *x, void *p );
static double _unur_rou_rectangle_aux_umin(double *x, void *p );
static double _unur_rou_rectangle_aux_umax(double *x, void *p );
/*---------------------------------------------------------------------------*/
/* Auxiliary functions used in the computation of the bounding rectangle     */
/*---------------------------------------------------------------------------*/


#define PDF(x)    _unur_cvec_PDF((x),(distr))    /* call to PDF              */

/*---------------------------------------------------------------------------*/


double
_unur_rou_rectangle_aux_vmax(double *x, void *p )
     /*----------------------------------------------------------------------*/
     /* Auxiliary function used in the computation of the bounding rectangle */
     /*----------------------------------------------------------------------*/
{
  struct ROU_RECTANGLE *rr;
  rr = p; /* typecast from void* to unur_rou_rectangle* */

  return -pow( _unur_cvec_PDF((x),(rr->distr)) ,
	       1./(1.+ rr->r * rr->dim) );
}

/*---------------------------------------------------------------------------*/

double
_unur_rou_rectangle_aux_umin(double *x, void *p)
     /*----------------------------------------------------------------------*/
     /* Auxiliary function used in the computation of the bounding rectangle */
     /*----------------------------------------------------------------------*/
{
  struct ROU_RECTANGLE *rr;
  rr = p; /* typecast from void* to unur_rou_rectangle* */

  return ( (x[rr->aux_dim] - rr->center[rr->aux_dim])
	   * pow( _unur_cvec_PDF((x),(rr->distr)),
		  rr->r / (1.+ rr->r * rr->dim) ) );
}

/*---------------------------------------------------------------------------*/

double
_unur_rou_rectangle_aux_umax(double *x, void *p)
     /*----------------------------------------------------------------------*/
     /* Auxiliary function used in the computation of the bounding rectangle */
     /*----------------------------------------------------------------------*/
{
  return (- _unur_rou_rectangle_aux_umin(x,p)) ;
}

/*---------------------------------------------------------------------------*/

int
_unur_rou_rectangle( struct ROU_RECTANGLE *rr )
     /*----------------------------------------------------------------------*/
     /* compute universal bounding rectangle                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{

  struct unur_funct_vgeneric faux; /* function to be minimized/maximized    */
  double *xstart, *xend, *xumin, *xumax; /* coordinate arrays used in maximum/minimum calculations */
  int d, dim; /* index used in dimension loops (0 <= d < dim) */
  int hooke_iters_vmax;  /* actual number of min/max iterations = return value of hooke()*/
  int hooke_iters_umin;  /* actual number of min/max iterations = return value of hooke()*/
  int hooke_iters_umax;  /* actual number of min/max iterations = return value of hooke()*/
  double scaled_epsilon; /* to be used in the hooke algorithm */


  /* dimension of the distribution */
  dim = rr->dim;

  /* allocate memory for the coordinate vectors */
  xstart = _unur_xmalloc(dim * sizeof(double));
  xend   = _unur_xmalloc(dim * sizeof(double));
  xumin  = _unur_xmalloc(dim * sizeof(double));
  xumax  = _unur_xmalloc(dim * sizeof(double));

  /* calculation of vmax */
      faux.f = (UNUR_FUNCT_VGENERIC*) _unur_rou_rectangle_aux_vmax;
      faux.params = rr;

      /* starting point */
      memcpy(xstart, rr->center, dim * sizeof(double));

      hooke_iters_vmax = _unur_hooke( faux, dim, xstart, xend,
				      ROU_HOOKE_RHO, ROU_HOOKE_EPSILON, ROU_HOOKE_MAXITER);

      rr->vmax = -faux.f(xend, faux.params);

      if (hooke_iters_vmax >= ROU_HOOKE_MAXITER) {
	 scaled_epsilon = ROU_HOOKE_EPSILON * rr->vmax;
	 if (scaled_epsilon>ROU_HOOKE_EPSILON) scaled_epsilon=ROU_HOOKE_EPSILON;

         /* recalculating extremum with scaled_epsilon and new starting point */
         memcpy(xstart, xend, dim * sizeof(double));
         hooke_iters_vmax = _unur_hooke( faux, dim, xstart, xend,
                              ROU_HOOKE_RHO, scaled_epsilon , ROU_HOOKE_MAXITER);
         rr->vmax = -faux.f(xend, faux.params);
         if (hooke_iters_vmax >= ROU_HOOKE_MAXITER) {
           _unur_warning("RECTANGLE" , UNUR_ERR_GENERIC, "Bounding rect uncertain (vmax)");
         }
      }

  /* calculation of umin and umax */

    for (d=0; d<dim; d++) {

      /* setting coordinate dimension to be used by the auxiliary functions */
      rr->aux_dim  = d;

      /* starting point at center */
      memcpy(xstart, rr->center, dim * sizeof(double));

      /*-----------------------------------------------------------------------------*/
      /* calculation for umin */

      faux.f = (UNUR_FUNCT_VGENERIC*) _unur_rou_rectangle_aux_umin;
      faux.params = rr;

      hooke_iters_umin = _unur_hooke( faux, dim, xstart, xend,
                           ROU_HOOKE_RHO, ROU_HOOKE_EPSILON, ROU_HOOKE_MAXITER);
      rr->umin[d] = faux.f(xend, faux.params);

      /* storing actual endpoint in case we need a recalculation */
      memcpy(xumin, xend, dim * sizeof(double));

      /*-----------------------------------------------------------------------------*/
      /* and now, an analogue calculation for umax */

      faux.f = (UNUR_FUNCT_VGENERIC*) _unur_rou_rectangle_aux_umax;
      faux.params = rr;

      hooke_iters_umax = _unur_hooke( faux, dim, xstart, xend,
                           ROU_HOOKE_RHO, ROU_HOOKE_EPSILON, ROU_HOOKE_MAXITER);
      rr->umax[d] = -faux.f(xend, faux.params);

      /* storing actual endpoint in case we need a recalculation */
      memcpy(xumax, xend, dim * sizeof(double));

      /*-----------------------------------------------------------------------------*/
      /* checking if we need to recalculate umin */
      if (hooke_iters_umin >= ROU_HOOKE_MAXITER) {
	 scaled_epsilon = ROU_HOOKE_EPSILON * (rr->umax[d]-rr->umin[d]);
	 if (scaled_epsilon>ROU_HOOKE_EPSILON) scaled_epsilon=ROU_HOOKE_EPSILON;

         /* recalculating extremum with scaled_epsilon and new starting point */
         faux.f = (UNUR_FUNCT_VGENERIC*) _unur_rou_rectangle_aux_umin;
         faux.params = rr;

         memcpy(xstart, xumin, dim * sizeof(double));
         hooke_iters_umin = _unur_hooke( faux, dim, xstart, xend,
                              ROU_HOOKE_RHO, scaled_epsilon , ROU_HOOKE_MAXITER);
         rr->umin[d] = faux.f(xend, faux.params);
         if (hooke_iters_umin >= ROU_HOOKE_MAXITER) {
           _unur_warning("RECTANGLE" , UNUR_ERR_GENERIC, "Bounding rect uncertain (umin)");
         }
      }

      /* checking if we need to recalculate umax */
      if (hooke_iters_umax >= ROU_HOOKE_MAXITER) {
	 scaled_epsilon = ROU_HOOKE_EPSILON * (rr->umax[d]-rr->umin[d]);
	 if (scaled_epsilon>ROU_HOOKE_EPSILON) scaled_epsilon=ROU_HOOKE_EPSILON;

         /* recalculating extremum with scaled_epsilon and new starting point */
         faux.f = (UNUR_FUNCT_VGENERIC*) _unur_rou_rectangle_aux_umax;
         faux.params = rr;

         memcpy(xstart, xumax, dim * sizeof(double));
         hooke_iters_umax = _unur_hooke( faux, dim, xstart, xend,
                              ROU_HOOKE_RHO, scaled_epsilon , ROU_HOOKE_MAXITER);
         rr->umin[d] = faux.f(xend, faux.params);
         if (hooke_iters_umax >= ROU_HOOKE_MAXITER) {
           _unur_warning("RECTANGLE" , UNUR_ERR_GENERIC, "Bounding rect uncertain (umax)");
         }
      }

      /*-----------------------------------------------------------------------------*/
      /* additional scaling of boundary rectangle */
      rr->vmax = rr->vmax * ( 1+ ROU_RECT_SCALING);
      rr->umin[d] = rr->umin[d] - (rr->umax[d]-rr->umin[d])*ROU_RECT_SCALING/2.;
      rr->umax[d] = rr->umax[d] + (rr->umax[d]-rr->umin[d])*ROU_RECT_SCALING/2.;

    }

  free(xstart); free(xend); free(xumin); free(xumax);

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_vnrou_rectangle() */


#undef PDF
#undef ROU_HOOKE_RHO
#undef ROU_HOOKE_EPSILON
#undef ROU_HOOKE_MAXITER
#undef ROU_RECT_SCALING

/*---------------------------------------------------------------------------*/

