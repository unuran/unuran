/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  Tests                                                                    *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <stdio.h>
#include <malloc.h>
#include <config.h>
#include <time.h>

#include <prng.h>
#include <unuran.h>
#include <unuran_tests.h>

/*---------------------------------------------------------------------------*/

#if UNUR_URNG_TYPE != UNUR_URNG_PRNG
#error UNUR_URNG_TYPE must be set to UNUR_URNG_PRNG in unuran_config.h
#endif

/*---------------------------------------------------------------------------*/

#define CHI_TEST_INTERVALS 100

/*---------------------------------------------------------------------------*/
/* enable/disable tests                                                      */
/* (comment out the methods or distributions that are not tested)            */

/* methods                                                                   */
#define T_DAU
#define T_DGT
#define T_SROU
#define T_STDR
#define T_UNIF

/* continuous univariate distributions                                       */
#define D_BETA
#define D_CAUCHY
#define D_GAMMA
#define D_NORMAL
#define D_UNIFORM

/* discrete univariate distributions                                         */
#define D_PV_RANDOM
#define D_PV_GEOMETRIC

/*---------------------------------------------------------------------------*/
/* global variables                                                          */

/* tresholds */
#define PVAL_LIMIT 1e-3             /* treshold for p-value for stat. test   */

/* distributions */
struct list_distr {
  UNUR_DISTR *distr;                /* pointer to distribution object        */
  unsigned    type;                 /* type of distribution                  */
  double      c_max;                /* maximal value for c to be T_c concave */
};

/* distribution types  */
#define T_Tconcave   0x00000001     /* univariate continuous T_c-concave     */
#define T_fpv        0x00000100     /* univariate discrete finite prob. vector */

extern struct list_distr *list_of_distr;  /* pointer to list of distributions */
extern int n_distr;                 /* number of distributions               */

/*---------------------------------------------------------------------------*/
/* True and false                                                            */

#ifndef TRUE
#define TRUE   (1)
#endif

#ifndef FALSE
#define FALSE  (0)
#endif

#ifndef M_LN10
#define M_LN10  2.302585092994046
#endif

/*---------------------------------------------------------------------------*/
/* make list of distributions                                                */
int make_list_of_distributions( struct list_distr **list_of_distr );

/*---------------------------------------------------------------------------*/
