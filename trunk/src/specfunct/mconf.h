/*
This file was considerably changed for "mycephes"

We were only interested in the files:
gamma.c
igamma.c
incbet.c
ndtr.c
ndtri.c

to use these files we needed the auxiliary files:

isnan.c
mtherr.c
polevl.c

The aim was to enhance portability for ANSI-C without expert knowledge
of the floating point unit.


The main change is that we changed the numeric constants (in cephes
included in const.c) from variables to defines and thus moved them
into this header file (mconf.h).
We changed these constants using the constants provided by ANSI-c in the
math.h file.
We changed also the definitions of MAXGAM (in incbet.c) and of MAXSTIR
(in gamma.c).

Everything concerning NANS and INFINITY was moved into isnan.c.
(Also the defines to turn it off or on.)
We have only tested this version with NAN and INFINITY turned off, but we have
not changed anything concerning that question. So turning this support on
should (could) work.

We have added the define in the next line to make it possible to suppress the error 
messages. 
*/
/*#define DEBUG */
/*If DEBUG is not defined the error messages are not printed, but nothing else
  was changed in mtherr.c  */

/*							mconf.h
 *
 *	Common include file for math routines
 *
 *
 *
 * SYNOPSIS:
 *
 * #include "mconf.h"
 *
 *
 *
 * DESCRIPTION:
 *
 * This file contains definitions for error codes that are
 * passed to the common error handling routine mtherr()
 * (which see).
 *
 * The file also includes a conditional assembly definition
 * for the type of computer arithmetic (IEEE, DEC, Motorola
 * IEEE, or UNKnown).
 * 
 * For Digital Equipment PDP-11 and VAX computers, certain
 * IBM systems, and others that use numbers with a 56-bit
 * significand, the symbol DEC should be defined.  In this
 * mode, most floating point constants are given as arrays
 * of octal integers to eliminate decimal to binary conversion
 * errors that might be introduced by the compiler.
 *
 * For little-endian computers, such as IBM PC, that follow the
 * IEEE Standard for Binary Floating Point Arithmetic (ANSI/IEEE
 * Std 754-1985), the symbol IBMPC should be defined.  These
 * numbers have 53-bit significands.  In this mode, constants
 * are provided as arrays of hexadecimal 16 bit integers.
 *
 * Big-endian IEEE format is denoted MIEEE.  On some RISC
 * systems such as Sun SPARC, double precision constants
 * must be stored on 8-byte address boundaries.  Since integer
 * arrays may be aligned differently, the MIEEE configuration
 * may fail on such machines.
 *
 * To accommodate other types of computer arithmetic, all
 * constants are also provided in a normal decimal radix
 * which one can hope are correctly converted to a suitable
 * format by the available C language compiler.  To invoke
 * this mode, define the symbol UNK.
 *
 * An important difference among these modes is a predefined
 * set of machine arithmetic constants for each.  The numbers
 * MACHEP (the machine roundoff error), MAXNUM (largest number
 * represented), and several other parameters are preset by
 * the configuration symbol.  Check the file const.c to
 * ensure that these values are correct for your computer.
 *
 * Configurations NANS, INFINITIES, MINUSZERO, and DENORMAL
 * may fail on many systems.  Verify that they are supposed
 * to work on your computer.
 */
/*
Cephes Math Library Release 2.3:  June, 1995
Copyright 1984, 1987, 1989, 1995 by Stephen L. Moshier
*/

/*#define DEBUG*/ 
/*If DEBUG is not defined no error messages are issued */

/* Constant definitions for math error conditions
 */

#define DOMAIN		1	/* argument domain error */
#define SING		2	/* argument singularity */
#define OVERFLOW	3	/* overflow range error */
#define UNDERFLOW	4	/* underflow range error */


/* UNKnown arithmetic, invokes coefficients given in
 * normal decimal format.  Beware of range boundary
 * problems (MACHEP, MAXLOG, etc. in const.c) and
 * roundoff problems in pow.c:
 * (Sun SPARCstation)
 */
#define UNK 1



/* Get ANSI function prototypes, if you want them. */
#if 1
#define ANSIPROT 1
int mtherr ( char *, int );
#else
int mtherr();
#endif

/* Variable for error reporting.  See mtherr.c.  */
extern int merror;


/*these defines were given as extern variables in const.c */
/*these were the old definitions for UNK */
/* double MACHEP =  1.11022302462515654042E-16;    2**-53 */
/* double MAXLOG =  7.08396418532264106224E2;      log 2**1022 */
/* double MINLOG = -7.08396418532264106224E2;      log 2**-1022 */
/* double MAXNUM =  1.79769313486231570815E308;     2**1024*(1-MACHEP) */
/* double PI     =  3.14159265358979323846;        pi */
/* double SQRTH  =  7.07106781186547524401E-1;     sqrt(2)/2 */


#include <math.h>
#define MACHEP  DBL_EPSILON   
#define MAXLOG  log(DBL_MAX)
#define MINLOG  log(DBL_MIN)
#define MAXNUM  DBL_MAX
#define PI       3.14159265358979323846   
#define SQRTH    7.07106781186547524401E-1



/*the below lines were transferred to isnan.c*/
/* so look there if you want to activate INFINITIES or NAN */

/* Define to ask for infinity support, else undefine. 
#undef INFINITIES 

 Define to ask for support of numbers that are Not-a-Number,
   else undefine.  This may automatically define INFINITIES in some files. 
#undef NANS 

#ifdef NANS
#define BIGENDIAN 0   0 for IBMPC, 1 for bigendian floating point unit 
 for UNURAN bigendian is only used for the function isnan() in isnan.c
   and only if NANS is dfined
#endif
#ifdef INFINITIES
double INFINITY = 1.0/0.0;  99e999; 
#else
double INFINITY =  1.79769313486231570815E308;     2**1024*(1-MACHEP) 
#endif
#ifdef NANS
double NAN = 1.0/0.0 - 1.0/0.0;
#else
double NAN = 0.0;
#endif


*/






















