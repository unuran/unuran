/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      cstd.c                                                       *
 *                                                                           *
 *   TYPE:      continuous univariate random variate                         *
 *   METHOD:    generators for standard distribution (from CRAND)            *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * ..... beschreibung ....                                                   *
 *                                                                           *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   author: Wolfgang.Hoermann @ statistik.wu-wien.ac.at                     *
 *           Josef.Leydold @ statistik.wu-wien.ac.at                         *
 *                                                                           *
 *   last modification: Thu Dec  2 14:52:19 CET 1999                         *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 1999 Wolfgang Hoermann and Josef Leydold                  *
 *   Dept. for Statistics, University of Economics, Vienna, Austria          *
 *                                                                           *
 *                                                                           *
 *   This library is free software; you can redistribute it and/or           *
 *   modify it under the terms of the GNU Library General Public             *
 *   License as published by the Free Software Foundation; either            *
 *   version 2 of the License, or (at your option) any later version.        *
 *                                                                           *
 *   This library is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU       *
 *   Library General Public License for more details.                        *
 *                                                                           *
 *   You should have received a copy of the GNU Library General Public       *
 *   License along with this library; if not, write to the Free              *
 *   Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.      *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <ctype.h>
#include <string.h>

#include <unur_methods.h>
#include <unur_stdgen.h>

#include <unur_cookies.h>
#include <unur_errno.h>
#include <unur_math.h>
#include <unur_utils.h>

/*---------------------------------------------------------------------------*/

#define GENTYPE "CSTD"         /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *cstd_create( struct unur_par *par );

static int cstd_parse_definition( struct unur_gen *gen, char *definition );
static int cstd_check_param( struct unur_gen *gen );

#if UNUR_DEBUG & UNUR_DB_INFO
static void cstd_info_init( struct unur_par *par, struct unur_gen *gen, int succeeded );
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define PAR     par->data.cstd
#define GEN     gen->data.cstd
#define SAMPLE  gen->sample.cont

/*---------------------------------------------------------------------------*/

/* distributions and methods                                                 */
/* #define UNUR_MASK_DISTR    0x000ffff0UL  --> unuran_defs.h                */

/* Exponential distribution                                                  */
#define DIST_EXP       0x00000010UL

/* Gamma distribution                                                        */
#define DIST_GAMMA     0x00000020UL

/* Normal distribution                                                       */
#define DIST_NORMAL    0x00000030UL

/*---------------------------------------------------------------------------*/

/* store id numbers in an array                                              */

#define UNUR_MAX_DIST_LEN     64   /* maximal length of string for distribution name */

struct cstd_id {
  char name[UNUR_MAX_DIST_LEN];    /* name of distribution                   */
  unsigned long id;                /* id number                              */
};

/* the following table should be sorted by the frequency of the use of the 
   corresponding distribution.                                               */
static struct cstd_id dist_table[] = {
  {"normal",      DIST_NORMAL},
  {"exponential", DIST_EXP},
  {"gamma",       DIST_GAMMA},
  
  {"",            0UL}             /* dummy: end of table .. */
};

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  User Interface                                                         **/
/*****************************************************************************/

struct unur_par *
 unur_cstd_new( char *definition )
/*---------------------------------------------------------------------------*/
/* get default parameters                                                    */
/*                                                                           */
/* parameters:                                                               */
/*   definition ... string containing description of distribution            */
/*                                                                           */
/* return:                                                                   */
/*   default parameters (pointer to structure)                               */
/*                                                                           */
/* error:                                                                    */
/*   return NULL                                                             */
/*---------------------------------------------------------------------------*/
{ 
  struct unur_par *par;

  /* allocate structure */
  par = _unur_malloc(sizeof(struct unur_par));
  COOKIE_SET(par,CK_CSTD_PAR);

  /* copy input */
  PAR.definition = definition;

  par->method   = UNUR_METH_CSTD;   /* method and default variant            */
  par->set      = 0UL;              /* inidicate default parameters          */    
  par->urng     = unur_get_default_urng(); /* use default urng               */

  _unur_set_debugflag_default(par); /* set default debugging flags           */
  _unur_set_genid(par,GENTYPE);     /* set generator identifier              */
  
  /* routine for starting generator */
  par->init = unur_cstd_init;

  return par;

} /* end of unur_cstd_new() */

/*****************************************************************************/

struct unur_gen *
 unur_cstd_init( struct unur_par *par )
/*---------------------------------------------------------------------------*/
/* initialize new generator                                                  */
/*                                                                           */
/* parameters:                                                               */
/*   params  pointer to paramters for building generator object              */
/*                                                                           */
/* return:                                                                   */
/*   pointer to generator object                                             */
/*                                                                           */
/* error:                                                                    */
/*   return NULL                                                             */
/*---------------------------------------------------------------------------*/
{ 
  struct unur_gen *gen;
  int succeeded;

  /* check arguments */
  CHECK_NULL(par,NULL);
  COOKIE_CHECK(par,CK_CSTD_PAR,NULL);

  /* create a new empty generator object */
  gen = cstd_create(par);
  if (!gen) { free(par); return NULL; }

  succeeded = cstd_parse_definition( gen, PAR.definition );

  if (succeeded)
    succeeded = cstd_check_param(gen);

#if UNUR_DEBUG & UNUR_DB_INFO
  /* write info into log file */
  if (gen->debug) cstd_info_init(par,gen,succeeded);
#endif

  /* free parameters */
  free(par);
  if (!succeeded) {
    unur_cstd_free(gen);
    return NULL;
  }

  return gen;

} /* end of unur_cstd_init() */

/*****************************************************************************/

double
 unur_cstd_sample( struct unur_gen *gen )
/*---------------------------------------------------------------------------*/
/* sample from generator                                                     */
/*                                                                           */
/* parameters:                                                               */
/*   gen ... pointer to generator object                                     */
/*                                                                           */
/* return:                                                                   */
/*   double (sample from random variate)                                     */
/*                                                                           */
/* error:                                                                    */
/*   return 0.                                                               */
/*---------------------------------------------------------------------------*/
{ 

  /* check arguments */
  CHECK_NULL(gen,0.);
  COOKIE_CHECK(gen,CK_CSTD_GEN,0.);

  switch (gen->method & UNUR_MASK_DISTR) {
  case DIST_NORMAL:
    return (GEN.dist_param[0] + GEN.dist_param[1] * nbm(gen->urng));
  case DIST_EXP:
    return ( -GEN.dist_param[0] * log(1. - _unur_call_urng(gen)) );
  case DIST_GAMMA:
    return (gammarand( GEN.dist_param[0],gen->urng ) * GEN.dist_param[1] + GEN.dist_param[2]);
  default:
    _unur_error(gen->genid,UNUR_ERR_INIT,"unknown distribution.");
    return 0;
  }

} /* end of unur_cstd_sample() */

/*****************************************************************************/

void
 unur_cstd_free( struct unur_gen *gen )
/*---------------------------------------------------------------------------*/
/* deallocate generator object                                               */
/*                                                                           */
/* parameters:                                                               */
/*   gen ... pointer to generator object                                     */
/*---------------------------------------------------------------------------*/
{ 

  /* check arguments */
  if( !gen ) /* nothing to do */
    return;

  /* magic cookies */
  COOKIE_CHECK(gen,CK_CSTD_GEN,/*void*/);
  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free memory */
  _unur_free_genid(gen);
  free(gen);

} /* end of unur_cstd_free() */

/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/

static struct unur_gen *
 cstd_create( struct unur_par *par )
/*---------------------------------------------------------------------------*/
/* allocate memory for generator                                             */
/*                                                                           */
/* parameters:                                                               */
/*   par ... pointer to parameter for building generator object              */
/*                                                                           */
/* return:                                                                   */
/*   pointer to (empty) generator object with default settings               */
/*                                                                           */
/* error:                                                                    */
/*   return NULL                                                             */
/*---------------------------------------------------------------------------*/
{
  struct unur_gen *gen;
  int i;

  /* check arguments */
  CHECK_NULL(par,NULL);
  COOKIE_CHECK(par,CK_CSTD_PAR,NULL);

  /* allocate memory for generator object */
  gen = _unur_malloc( sizeof(struct unur_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_CSTD_GEN);

  /* routines for sampling and destroying generator */
  SAMPLE = unur_cstd_sample;
  gen->destroy = unur_cstd_free;

  /* copy some parameters into generator object */
  _unur_copy_urng_pointer(par,gen);  /* pointer to urng into generator object*/
  _unur_copy_debugflag(par,gen);     /* copy debugging flags into generator object */
  _unur_copy_genid(par,gen);         /* copy generator identifier            */

  /* indicates method and variant 
     the bits for indicating distribution are cleared */
  gen->method = par->method & ~UNUR_MASK_DISTR;

  /* init parameters of distribution */
  GEN.n_dist_param = 0;
  for (i=0; i<UNUR_MAX_DIST_PARAMS; i++)
      GEN.dist_param[i] = 0.;

  /* return pointer to (almost empty) generator object */
  return(gen);
  
} /* end of cstd_create() */

/*****************************************************************************/
/**  Definition parsing code                                                **/
/*****************************************************************************/

static int
 cstd_parse_definition( struct unur_gen *gen, char *def )
/*---------------------------------------------------------------------------*/
/* parse definition string, split into its components                        */
/*                                                                           */
/* parameters:                                                               */
/*   gen  ... pointer to generator object                                    */
/*   def   ... definition string like "normal(0.,1.)"                        */
/*                                                                           */
/* return:                                                                   */
/*   1 ... if successful                                                     */
/*   0 ... otherwise                                                         */
/*                                                                           */
/* error:                                                                    */
/*   return 0                                                                */
/*                                                                           */
/* comment:                                                                  */
/*   have the following components:                                          */
/*      string ( number [, number [, number [...]]] )                        */
/*   where the string is interpreted as the name of the distribution and     */
/*   number as a parameter.                                                  */
/*                                                                           */
/* warning:                                                                  */
/*   the input is not not checked against typos.                             */
/*   case sensitive: always use lower case letters.                          */
/*---------------------------------------------------------------------------*/
{
  static char *string = NULL;
  static int len_string = 0;

  char *dist, *param;
  int len_def;
  int n_param;
  int i;

  /* check arguments */
  CHECK_NULL(gen,0.);
  COOKIE_CHECK(gen,CK_CSTD_GEN,0.);

  /* check for empty input */
  len_def = strlen(def);
  if (!len_def) {
    _unur_error(gen->genid,UNUR_ERR_INIT,"definition string empty.");
    return 0;
  }

  /* we have to copy the definition string into a temporary array */
  if (len_string < len_def) {
    /* temporary array too short */
    len_string = 2 * len_def;    /* to avoid too many malloc calls */
    free(string);                   /* would realloc be faster? */
    string = _unur_malloc( (len_string + 1) * sizeof(char) );
  }
  strcpy(string,def);

  /* tokenize definition string */

  /* name of distribution */
  dist = strtok(string, " (");

  /* parameters of distribution */
  for( n_param = 0; n_param <= UNUR_MAX_DIST_PARAMS; n_param++ ) {
    param = strtok(NULL, " (,)");
    if (!param) break;

    /* check number of parameters */
    if (n_param == UNUR_MAX_DIST_PARAMS) {
      _unur_error(gen->genid,UNUR_ERR_INIT,"too many parameters.");
      return 0;
    }

    /* store parameter */
    GEN.dist_param[n_param] = strtod(param,NULL);
  }

  /* store number of parameters in structure */
  GEN.n_dist_param = n_param;

  /* get id number of distribution */
  i = 0;
  gen->method &= ~UNUR_MASK_DISTR;     /* clear bits for indicating distribution */
  while(dist_table[i].name[0] != '\0') {  /* empty string ? */
    if (strcmp(dist,dist_table[i].name) == 0) {
      /* non ISO C alternative: strcasecmp for non case-sensitive comparison */
      gen->method |= dist_table[i].id;
      break;
    }
    i++; /* try next distribution */
  }

  /* distribution known ? */
  if (!(gen->method & UNUR_MASK_DISTR)) {
    _unur_error(gen->genid,UNUR_ERR_INIT,"unkown distribution.");
    return 0;
  }

  /* o.k. */
  return 1;
} /* end of cstd_parse_definition() */

/*****************************************************************************/
/**  Check input                                                            **/
/*****************************************************************************/

static int
 cstd_check_param( struct unur_gen *gen )
/*---------------------------------------------------------------------------*/
/* check parameters for distribution                                         */
/* if no optional parameters are given, use defaults                         */
/* if essential parameter is missing, make error message and abort           */
/*                                                                           */
/* parameters:                                                               */
/*   gen ... pointer to generator object                                     */
/*                                                                           */
/* return:                                                                   */
/*   1 ... if successful                                                     */
/*   0 ... otherwise                                                         */
/*                                                                           */
/* error:                                                                    */
/*   return 0                                                                */
/*---------------------------------------------------------------------------*/
{
  switch (gen->method & UNUR_MASK_DISTR) {

  case DIST_NORMAL:
    /* set defaults */
    switch (GEN.n_dist_param) {
    case 0: GEN.dist_param[0] = 0.;    /* mean (location)   */
    case 1: GEN.dist_param[1] = 1.;    /* std. dev. (scale) */
    case 2:
      break;
    default:
      _unur_error(gen->genid,UNUR_ERR_INIT,"normal: invalid number of parameters.");
      return 0;
    }
    GEN.n_dist_param = 2;
    /* check parameters */
    if (GEN.dist_param[1] <= 0.) {
      _unur_error(gen->genid,UNUR_ERR_INIT,"normal: std deviation [1] <= 0.");
      return 0;
    }
    return 1;

  case DIST_EXP:
    /* set defaults */
    switch (GEN.n_dist_param) {
    case 0: GEN.dist_param[0] = 1.;    /* scale    */
    case 1: GEN.dist_param[1] = 0.;    /* location */
    case 2:
      break;
    default:
      _unur_error(gen->genid,UNUR_ERR_INIT,"exponential: invalid number of parameters.");
      return 0;
    }
    GEN.n_dist_param = 2;
    /* check parameters */
    if (GEN.dist_param[0] <= 0.) {
      _unur_error(gen->genid,UNUR_ERR_INIT,"exponential: scale parameter [0] <= 0.");
      return 0;
    }
    return 1;

  case DIST_GAMMA:
    /* set defaults */
    switch (GEN.n_dist_param) {
    case 0: GEN.dist_param[0] = 1.;    /* shape    */
    case 1: GEN.dist_param[1] = 1.;    /* scale    */
    case 2: GEN.dist_param[2] = 0.;    /* location */
    case 3:
      break;
    default:
      _unur_error(gen->genid,UNUR_ERR_INIT,"gamma: invalid number of parameters.");
      return 0;
    }
    GEN.n_dist_param = 3;
    /* check parameters */
    if (GEN.dist_param[0] <= 0.) {
      _unur_error(gen->genid,UNUR_ERR_INIT,"gamma: shape parameter [0] <= 0.");
      return 0;
    }
    if (GEN.dist_param[1] <= 0.) {
      _unur_error(gen->genid,UNUR_ERR_INIT,"gamma: scale parameter [1] <= 0.");
      return 0;
    }
    return 1;

  default:
    _unur_error(gen->genid,UNUR_ERR_INIT,"unknown distribution.");
    return 0;
  }

} /* end of cstd_check_param() */

/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

#if UNUR_DEBUG & UNUR_DB_INFO

static void
 cstd_info_init( struct unur_par *par, struct unur_gen *gen, int succeeded )
/*---------------------------------------------------------------------------*/
/* write info about generator into logfile                                   */
/*                                                                           */
/* parameters:                                                               */
/*   par ... pointer to parameter for building generator object              */
/*   gen ... pointer to generator object                                     */
/*---------------------------------------------------------------------------*/
{
  FILE *log;
  int i;
  char *dist_name;

  log = unur_get_log();

  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(log,"%s: method  = generator for standard distribution\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

  fprintf(log,"%s: sampling routine = unur_cstd_sample()\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

  /* distribution */
  fprintf(log,"%s: distribution:\n",gen->genid);
  fprintf(log,"%s:    definition (input): \"%s\"\n",gen->genid,PAR.definition);

  /* get name of distribution */
  i = 0;
  dist_name = NULL;
  while(dist_table[i].id != 0UL) {
    if ( dist_table[i].id == (gen->method & UNUR_MASK_DISTR) ) {
      dist_name = dist_table[i].name;
      break;
    }
    i++;
  }
  
  if (dist_name) {
    fprintf(log,"%s:    name = %s\n",gen->genid,dist_name);
    fprintf(log,"%s:    parameters: %d",gen->genid,GEN.n_dist_param);
    if( !succeeded ) 
      fprintf(log,"\t INVALID!\n");
    else
      fprintf(log,"\n");
    for (i=0; i<GEN.n_dist_param; i++)
      fprintf(log,"%s:\tparam[%d] = %g\n",gen->genid,i,GEN.dist_param[i]);
  }
  else
    fprintf(log,"%s:    UNKNOWN distribution\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

  if( !succeeded )
    fprintf(log,"%s: INIT failed!\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

} /* end of cstd_info_init() */

#endif

/*****************************************************************************/

