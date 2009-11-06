
#if defined(HAVE_LIBRMATH)
/* Rmath standalone library from R project. */
#  define MATHLIB_STANDALONE
#  include <Rmath.h>
#  define HAVE_R_FUNCTIONS
#elif defined(R_UNURAN)
/* Rmath for 'Runuran' */
#  include <Rmath.h>
#  define HAVE_R_FUNCTIONS
#endif

