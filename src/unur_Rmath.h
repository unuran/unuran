
/*---------------------------------------------------------------------------*/
/*                                                                           */
/*  Use Rmath library from R project.                                        */
/*                                                                           */
/*---------------------------------------------------------------------------*/

#ifdef HAVE_LIBRMATH

/* we have to distinguish between two cases:                                 */
#  ifdef R_UNURAN
/*   Rmath for 'Runuran': nothin special to do. */
#  else
/*   Rmath standalone library. */
#    define MATHLIB_STANDALONE
#  endif

/* include header file */
#  include <Rmath.h>

#endif

/*---------------------------------------------------------------------------*/
