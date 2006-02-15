/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: mvtdr_struct.h                                                    *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares structures for method MVTDR                              *
 *         (Multi-Variate Transformed Density Rejection)                     *
 *                                                                           *
 *****************************************************************************
     $Id$
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2000 Wolfgang Hoermann and Josef Leydold                  *
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
 *****************************************************************************/


/**
  structures
**/

typedef struct s_vertex           /* data for one vertex */
{
  struct s_vertex *next;          /* pointer to next vertex in list */
  int index;                      /* index of vertex */
  double *coord;                /* coordinates of spanning vector(norm = 1), "vertex" */
  double norm;                    /* norm of vertex */
} VERTEX;


typedef struct s_cone             /* a cone */
{
  struct s_cone *next;            /* pointer to next cone in list */
/* #if DEBUG > 0 */
/*   int index;                      /\* index of cone *\/ */
/* #endif */
  int level;                      /* level of triangulation */
  VERTEX **v;                   /* list of vertices of the cone */
  double *center;               /* barycenter of cone */
  double logdetf;                    /* log determinant - log((dim-1)!) for cone */
  double alpha;                   /* parameters for hat function */
  double beta;
  double *gv;                   /* <g,v> for all vertices v */
  double logai;                      /* (log of) coefficient for marginal density */
  double tp;                      /* coordinate of touching point */
  double Hi;                      /* volume under hat in cone */
  double Hsum;                    /* accumulated sum of volumes */
  double fp;                      /* value of density at touching point */
/* #if RECTANGLE == 1 */
/*   double height;                  /\* height of pyramid *\/ */
/*   double tdrg_vol;                /\* volume below hat of gamma density *\/ */
/* #endif */
/* #if DEBUG > 1000 */
/*   double fp;                      /\* value of density at touching point *\/ */
/*   double g[N];                    /\* direction of sweep-plane *\/ */
/*   double Vi;                      /\* volume of triangle of distance 1 *\/ */
/* #endif */
} CONE;


typedef struct s_edge_table       /* hash table for edges */
{
  int  index[2];                  /* index of incident vertices */
  VERTEX *vertex;                 /* index of corresponding vertex (=barycenter) */
  struct s_edge_table *next;
} E_TABLE;

typedef struct s_tp_arg           /* argument for tp function */
{
  double t;                      /* touching point */
  double logH;                  /* log of volume below hat */
  CONE *c;                        /* parameters */
  UNUR_GEN *gen;
  int valid;                   /* bolean to store if x in domain */  /* remove ! */
} TP_ARG;


/*---------------------------------------------------------------------------*/
/* Information for constructing the generator                                */

struct unur_mvtdr_par { 

/* #if RECTANGLE == 1 */
/*   double rl[N];                   /\* lower bounds for rectangle *\/ */
/*   double ru[N];                   /\* upper bounds for rectangle *\/ */
/* #endif */

  int max_cones;                  /* maximum number of cones (at least 2^(N+T_STEPS_MIN) */
  int steps_min;                  /* minimum number of triangulation steps */

#if MODE == 1
  double mode_to_boundary;        /* move mode to boundary if |mode - boundary| / length < MODE_TO_BOUNDARY */
#endif

};


/*---------------------------------------------------------------------------*/
/* The generator object                                                      */

struct unur_mvtdr_gen { 
  int    dim;               /* dimension of distribution                     */

  const double *center;          /* center of distribution */

  CONE *cone;                        /* root of list of cones */
  CONE *last_cone;                   /* pointer to last cone in list */
  int n_cone;                      /* number of cones */

  VERTEX *vertex;                      /* root of list of vertices */
  VERTEX *last_vertex;                 /* pointer to last vertex in list */
  int n_vertex;                      /* number of vertices */

  E_TABLE **etable;                   /* pointer to edge table */
  int etable_size;                    /* size of edge table */

  CONE **guide;                   /* pointer to guide table */
  int guide_size;                 /* size of guide table */

  double *S;                      /* working array for storing point on simples */
  double *g;                      /* working array for vector g (direction of sweeping plane) */
  double *tp_coord;               /* working array for storing coordinates of touching point of hat */
  double *tp_mcoord;              /* working array for storing coordinates of touching point of hat moved into center */
  double *tp_Tgrad;               /* working array for storing gradient of transformed density at tp */

  double Htot;                    /* total volume below hat */
  int steps_min;                  /* minimum number of triangulation steps */
  int n_steps;                  /* (highest) number of triangulation steps */

/* #if RECTANGLE == 1                /\* rectangle moved to mode == origin *\/ */
/*   double rl[N];                   /\* lower bounds for rectangle *\/ */
/*   double ru[N];                   /\* upper bounds for rectangle *\/ */
/* #endif */

/* #if RECTANGLE == 1 */
/*   double max_gamma;               /\* maximum value for gamma variaties *\/ */
/* #endif */
/* #if DEBUG */
/*   double volume_density;          /\* volume below density (for debugging only) *\/ */
/*   int db_n_rpoint_accept;         /\* number of accepted random points (for debugging only) *\/ */
/*   int db_n_rpoint_reject;         /\* number of rejected random points (for debugging only) *\/ */
/* #if RECTANGLE == 1 */
/*   int db_n_rpoint_reject_outside; /\* number of rejected random points outside domain (for debugging only) *\/ */
/* #endif */

};

/*---------------------------------------------------------------------------*/
