/**
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

/********************************************************************/
/**                                                                **/
/**   cranduni.c                                                   **/
/**                                                                **/
/********************************************************************/              

/*------------------------------------------------------------------*/

/* the header files */
#include <stdio.h>

/*------------------------------------------------------------------*/

/* Type of generator (0-3, see below                                */
#define UNIFORM 2                     

/* value for the start seed (for details see below)                 */
#define START_SEED 1

/*------------------------------------------------------------------*/

static unsigned long int x = START_SEED;          /* standard seed  */

/*------------------------------------------------------------------*/

/*------------------------------------------------------------------*/
#if UNIFORM == 0
/*------------------------------------------------------------------*/

/*************************************************************
 * Ahrens J.H., Dieter U., Grube A. (1970)                   *
 * m = 2^32, a = 663608941, c = 0                            *
 * Pseudo-random numbers a new proposal for the choice of    *
 * multiplicators, Computing 6, 1970, P. 121-138             *
 *************************************************************/

double uniform()
{
  x*=663608941;
  return(x*2.328306436538696e-10);
}

/* This function allows to start the uniform generator with  *
 * a certain seed. The seed must be of the form 4k+1,        *
 * otherwise it is changed                                   */

void uinit(n)
unsigned long int n;
{
  x=4*(n/4)+1;
  printf("UINIT 0\n");
}

/* This function returns the current state of the seed       */

unsigned long int uget()
{ printf("uget0\n");
  return(x);
}
/*------------------------------------------------------------------*/
#endif
/*------------------------------------------------------------------*/
#if UNIFORM == 1
/*------------------------------------------------------------------*/

/*************************************************************
 * Marsaglia G. (1972), m = 2^32, a = 69069, c = 1           *
 * The structure of linear congruential sequences, in:       *
 * Applications of Number Theory to Numerical Analysis, S.K. *
 * Zaremba, ed., Academic Press, New York 1972               *
 *************************************************************/

double uniform()
{
  x=(69069*(x)+1);
  return(x*2.328306436538696e-10+1.164153218269348e-10);
}

/* This function allows to start the uniform generator with  */
/* a certain seed                                            */

void uinit(n)
unsigned long int n;
{
  x=n;
  printf("UINIT 1\n");
}

/* This function returns the current state of the seed       */

unsigned long int uget()
{ printf("uget1\n");
  return(x);
}
/*------------------------------------------------------------------*/
#endif
/*------------------------------------------------------------------*/
#if UNIFORM == 2
/*------------------------------------------------------------------*/

/*************************************************************
 * Fishman G.S., Moore L.R. (1986)                           *
 * m = 2^31-1, a = 742938285, c = 0                          *
 * An exhaustive analysis of multiplicative congruential     *
 * number generators with modulus 2^31-1, SIAM Journal Sci.  *
 * Stat. Comput., 1986, P. 24-45                             *
 *************************************************************/

#define A 742938285
#define AHI (A>>15)
#define ALO (A&0x7FFF)

double uniform()
{
  unsigned long int xhi,xlo,mid;

  xhi= x>>16;
  xlo= x&0xFFFF;
  mid=AHI*xlo+(ALO<<1)*xhi;
  x=AHI*xhi+(mid>>16)+ALO*xlo;
  if (x&0x80000000) x-=0x7FFFFFFF;
  x+=((mid&0xFFFF)<<15);
  if (x&0x80000000) x-=0x7FFFFFFF;
  return (x*4.656612875245797e-10);
}

/* This function allows to start the uniform generator with  */
/* a certain seed. The seed is ok if it is greater equal 1   */
/* and lower equal  2^{31}-2 !!!                             */

void uinit(n)
unsigned long int n;
{
  n%=0x7FFFFFFF; /* n wird mod 2^{31}-1 gerechnet */
  if (n==0) n=1;
  x=n;
  printf("UINIT 2\n");
}

/*This function returns the current state of the seed        */

unsigned long int uget()
{ printf("uget2\n");
  return(x);
}
/*------------------------------------------------------------------*/
#endif
/*------------------------------------------------------------------*/
#if UNIFORM == 3
/*------------------------------------------------------------------*/

/*************************************************************
 * L'Ecuyer P., Blouin F., (1989)                            *
 * m = 2^31-1                                                *
 * multiple recursive linear congruential generator          *
 * Multiple Recursive and Matrix Linear Congruential         *
 * Generators, Université Laval, Ste-Foy, Quebec,Canada 1989 *
 *************************************************************/

static unsigned long int xurn1=1, xurn2=69070, xurn3=475628535,
                         xurn4=1129920461, xurn5=772999773;
#define A1 107374182
#define A1HI (A1>>15)
#define A1LO (A1&0x7FFF)
#define A5 104480
#define A5HI (A5>>15)
#define A5LO (A5&0x7FFF)

double uniform()
{
  unsigned long int h,xurnhi,xurnlo,mid;

  STEP_COUNTER;

  xurnhi=xurn5>>16;
  xurnlo=xurn5&0xFFFF;
  mid=A5HI*xurnlo+(A5LO<<1)*xurnhi;
  xurn5=A5HI*xurnhi+(mid>>16)+A5LO*xurnlo;
  if (xurn5&0x80000000) xurn5-=0x7FFFFFFF;
  xurn5+=((mid&0xFFFF)<<15);
  if (xurn5&0x80000000) xurn5-=0x7FFFFFFF;
  h=xurn5;
  xurn5=xurn4;xurn4=xurn3;xurn3=xurn2;xurn2=xurn1;
  xurnhi= xurn1>>16;
  xurnlo= xurn1&0xFFFF;
  mid=A1HI*xurnlo+(A1LO<<1)*xurnhi;
  xurn1=A1HI*xurnhi+(mid>>16)+A1LO*xurnlo;
  if (xurn1&0x80000000) xurn1-=0x7FFFFFFF;
  xurn1+=((mid&0xFFFF)<<15);
  if (xurn1&0x80000000) xurn1-=0x7FFFFFFF;
  xurn1+=h;
  if (xurn1&0x80000000) xurn1-=0x7FFFFFFF;

  return ((double)xurn1*4.656612875245797e-10+1.32e-10);
}

/* This function allows to start the uniform generator with  */
/* a single value to make it compatible with uinit of        */
/* uniform 0, 1 and 2.                                       */
/* from the one value passed to uinit consecutive values are */
/* computed to initialize xurn1 to xurn5.                    */

void uinit(n)
unsigned long int n;
{
  xurn1=n%0x7FFFFFFF;n=(n*69069+1)&0xFFFFFFFF;
  xurn2=n%0x7FFFFFFF;n=(n*69069+1)&0xFFFFFFFF;
  xurn3=n%0x7FFFFFFF;n=(n*69069+1)&0xFFFFFFFF;
  xurn4=n%0x7FFFFFFF;n=(n*69069+1)&0xFFFFFFFF;
  xurn5=n%0x7FFFFFFF;
}

/* This function allows to start the uniform generator with  */
/* a complete seed; it passes an array of                    */
/* 5 unsigned long ints; the elements of the array are ok if */
/* they are greater equal 1 or lower equal  2^{31}-2 !!!     */

void uinitl(n)
unsigned long int n[5];
{
  int i;
  for(i=0;i<=4;++i)
     {n[i]%=0x7FFFFFFF;/*n wird mod 2^{31}-1 gerechnet*/
      if (n[i]==0)
          n[i]=1;
     }
  xurn1=n[0];xurn2=n[1];xurn3=n[2];xurn4=n[3];xurn5=n[4];
}

/* This function puts the current state into an array of 5   */
/* unsigned long ints                                        */

void ugetl(unsigned long int res[])

{ printf("uget3\n");
  printf("xurn1 = %ld\n",xurn1);res[0]=xurn1;
  printf("xurn2 = %ld\n",xurn2);res[1]=xurn2;
  printf("xurn3 = %ld\n",xurn3);res[2]=xurn3;
  printf("xurn4 = %ld\n",xurn4);res[3]=xurn4;
  printf("xurn5 = %ld\n",xurn5);res[4]=xurn5;
}

/* This function returns the current state of the first seed       */
/* for compatibility with generator 0 to 2*/
unsigned long int uget()
{ printf("uget3\n");
  return(xurn1);
}

/*------------------------------------------------------------------*/
#endif
/*------------------------------------------------------------------*/




