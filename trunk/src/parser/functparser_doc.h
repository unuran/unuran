
/*---------------------------------------------------------------------------*/
/*

=NODE  StringFunct    Function String
=UP StringAPI [25]

=DESCRIPTION

In unuran it is also possible to define functions (e.g. CDF or PDF) as
strings. As you can see in Example 2 (@ref{Example_2_str}) it is very
easy to define the PDF of a distribution object by means of a string. 
The possibilities using this string interface are more restricted than
using a pointer to a routine coded in C (@ref{Example_2}). 
But the differences in evaluation time is small.
When a distribution object is defined using this string interface then
of course the same conditions on the given density or CDF must be
satisfied for a chosen method as for the standard API.
This string interface can be used for both within the UNURAN string
API using the unur_str2gen() call, and for calls that define the
density or CDF for a particular distribution object as done with
(e.g.) the call unur_distr_cont_set_pdfstr().
Here is an example for the latter case:

@smallexample
  unur_distr_cont_set_pdfstr(distr,"1-x*x");
@end smallexample


@subheading Syntax

The syntax for the function string is case insensitive, white space is
ingnored. The expressions are similar to most programming languages
and mathamatical programs (see also the examples below). Especially it
is influenced by C. The usual preceedence rules are used (from highest
to lowest preceedence: parenthesis, functions, power, multiplication,
addition, relation operators).
Relation operators can be used as indicator functions, i.e. the term
@code{(x>1)} is evaluted as @code{1} if this relation is satisfied,
and as @code{0} otherwise.

The first unknown symbol (letter or word) is interpreted as the
variable of the function. It is recommended to use @code{x}.
Only one variable can be used.


@subheading List of symbols

@cartouche
@noindent
@b{Numbers and constants}
@sp 1
@multitable   {xxxxxxxxxxx} {xxxxxxxxxxxxxxxxxxxxxxxxxxxxx} {xxxxxxxxxxxxxxxxxxxxxxx} 
@item @code{[numbers]} @tab numbers @tab 123   
@item @code{.}     @tab comma       @tab 165.567  
@item @code{e}     @tab exponet     @tab  123e-7
@end multitable
@end cartouche

@sp 1
@cartouche
@noindent
@b{Relation operators (Indicator functions)}
@sp 1
@multitable {xxxxxxxxxxx} {xxxxxxxxxxxxxxxxxxxxxxxxxxxxx} {xxxxxxxxxxxxxxxxxxxxxxx}
@item Symbol        @tab Explanation              @tab Examples 
@item @code{<}      @tab @i{less than}            @tab @code{(x<1)}
@item @code{=}      @tab @i{equal}                @tab @code{(2=x)}
@item @code{==}     @tab @i{same as} @code{=}     @tab @code{(x==3)}
@item @code{>}      @tab @i{greater than}         @tab @code{(x>0)}
@item @code{<=}     @tab @i{less than or equal}   @tab @code{(x<=1)}      
@item @code{!=}     @tab @i{not equal}            @tab @code{(x!0)}
@item @code{<>}     @tab @i{same as} @code{!=}    @tab @code{(x<>pi)}
@item @code{>=}     @tab @i{greater or equal}     @tab @code{(x>=1)}
@end multitable
@end cartouche

@sp 1
@cartouche
@noindent
@b{Special symbols}
@sp 1
@multitable {xxxxxxxxxxx} {xxxxxxxxxxxxxxxxxxxxxxxxxxxxx} {xxxxxxxxxxxxxxxxxxxxxxx}
@item @code{(}   @tab open bracket @tab 2*(3+x)      
@item @code{)}   @tab close bracket @tab 2*(3+x)
@item @code{,}   @tab 
@end multitable
@end cartouche

@sp 1
@cartouche
@noindent
@b{Arithmetic operators}
@sp 1
@multitable  {xxxxxxxxxxx} {xxxxxxxxxxxxxxxxxxxxxxxxxxxxx} {xxxxxxxxxxxxxxxxxxxxxxx}
@item @code{+}      @tab @i{addition}             @tab @code{2+x}
@item @code{-}      @tab @i{subtraction}          @tab @code{2-x}
@item @code{*}      @tab @i{multiplication}       @tab @code{2*x}
@item @code{/}      @tab @i{division}             @tab @code{x/2}
@item @code{^}      @tab @i{power}                @tab @code{x^2}
@end multitable
@end cartouche

@sp 1
@cartouche
@noindent
@b{system constants}
@sp 1
@multitable  {xxxxxxxxxxx} {xxxxxxxxxxxxxxxxxxxxxxxxxxxxx} {xxxxxxxxxxxxxxxxxxxxxxx}
@item @code{pi}     @tab @i{pi = 3.1415@dots{}}   @tab @code{3*pi+2}
@item @code{e}      @tab @i{Euler's constant = 2.7182@dots{}} 
                                                  @tab @code{3*e+2}
  						       (do not cofuse with @code{3e2})
@end multitable
@end cartouche

@sp 1
@cartouche
@noindent
@b{Functions}
@sp 1
@multitable  {xxxxxxxxxxx} {xxxxxxxxxxxxxxxxxxxxxxxxxxxxx} {xxxxxxxxxxxxxxxxxxxxxxx}
@item @code{mod}    @tab @code{mod(m,n)} @i{ gives the remainder on devision m by n}   @tab mod(x,2)   
@item @code{exp}    @tab @i{exponential function 
                         (same as @code{e^x})}    @tab @code{exp(-x^2)} (same as @code{e^(-x^2)})
@item @code{log}    @tab @i{natural logarithm}    @tab @code{log(x)}
@item @code{sin}    @tab @i{sine}                 @tab @code{sin(x)}            
@item @code{cos}    @tab @i{cosine}               @tab @code{cos(x)}
@item @code{tan}    @tab @i{tangent}              @tab @code{tan(x)}
@item @code{sec}    @tab @i{secant}               @tab @code{sec(x*2)}
@item @code{sqrt}   @tab @i{square root}          @tab @code{sqrt(2*x)}
@item @code{abs}    @tab @i{absolute value}       @tab @code{abs(x)}
@item @code{sgn}    @tab @i{sign function}        @tab @code{sign(x)*3} 
@end multitable
@end cartouche

@sp 1
@cartouche
@noindent
@b{Variable}
@sp 1
@multitable  {xxxxxxxxxxx} {xxxxxxxxxxxxxxxxxxxxxxxxxxxxx} {xxxxxxxxxxxxxxxxxxxxxxx}
@item @code{x}         @tab variable         @tab 3*x^2
@end multitable
@end cartouche


Especially with relation operators it is very easy to define functions in pieces :
So you can define the piecewise linear function defined throug the points (-1,0),(0,1),(1,0) with 
@example
(x>-1)*(x<0)*(1+x) + (x>=0)*(x<1)*(1-x) .
@end example  



@sp 2
Several complex examples for string-defined 'distributions':

@smallexample 
-4.7285e-7*x

2.894736*10^2

(sin( log(3*x*(cos( 3*x^3-4.6789/(x+4)]))))-1

exp(x^2)*(sin(x*cos(x^2-1))+1)*((x-3*pi*x)<1)

(3*(2<>x)and(x>2))+x

(x>-1)*(x<0)*(1+x) + (x>=0)*(x<1)*(1-x
@end smallexample

=EON
*/
/*---------------------------------------------------------------------------*/





/* junk ....

@subheading symbols and examples

The string isn't case-sensitive and use the usual precendence (the operators in the table below are sorted from highest to lowest precedance).
The first unknown character string will be interpreted as the variable. Parameters must written as real values.

                    
In the following table you can see all symbols you can use sorted by symbol groups:

@sp 1
@multitable   {xxxxxxxxxxx} {xxxxxxxxxxxxxxxxxxxxxxxxxxxxx} {xxxxxxxxxxxxxxxxxxxxxxx}
@item         symbols   @tab explanation @tab examples 
@end multitable


@sp 1
@noindent 
@multitable   {xxxxxxxxxxx} {xxxxxxxxxxxxxxxxxxxxxxxxxxxxx} {xxxxxxxxxxxxxxxxxxxxxxx} 
@item         @b{[figures]} @tab numbers @tab 123   
@item         @b{. -}     @tab comma, negative numbers       @tab -165.567  
@item         @b{e}     @tab exponet     @tab  123e-7
@item         @b{pi}       @tab  pi: 3,1415...   @tab 3*pi+2
@item         @b{e}        @tab  exponential constant: 2,7182... @tab 3*e+2   
@end multitable

@sp 1
@noindent
@b{special functions}
@multitable  {xxxxxxxxxxx} {xxxxxxxxxxxxxxxxxxxxxxxxxxxxx} {xxxxxxxxxxxxxxxxxxxxxxx}
@item         @b{exp}      @tab exponential function @tab exp(-x^2)
@item         @b{log}      @tab natural logarithm    @tab log(x)
@item         @b{sin}      @tab sine                 @tab sin(x)             
@item         @b{cos}      @tab cosine               @tab cos(x)
@item         @b{tan}      @tab tangent              @tab tan(x)
@item         @b{sec}      @tab secant               @tab sec(x*2)
@item         @b{sqrt}     @tab square root          @tab sqrt(2*x)
@item         @b{abs}      @tab absolute value       @tab abs(x)
@item         @b{sgn}      @tab sign function        @tab sign(x)*3 
@end multitable 

@sp 1
@noindent
@multitable {xxxxxxxxxxx} {xxxxxxxxxxxxxxxxxxxxxxxxxxxxx} {xxxxxxxxxxxxxxxxxxxxxxx} 
@item        @b{(}   @tab left parenthesis @tab 2*(3+x)      
@item        @b{)}   @tab right parenthesis @tab 2*(3+x)
@end multitable

@sp 1
@noindent
@b{relation operators}
@multitable {xxxxxxxxxxx} {xxxxxxxxxxxxxxxxxxxxxxxxxxxxx} {xxxxxxxxxxxxxxxxxxxxxxx} 
@item       @b{<}         @tab less             @tab 3*(x<1) 
@item       @b{=}, @b{==} @tab equal            @tab 3*(2=x), 3*(2==x)
@item       @b{>}         @tab greater          @tab 3*(x>1)
@item       @b{<=}        @tab less or equal    @tab 3*(x<=1)      
@item       @b{<> , !=}   @tab not equal        @tab 3*(x<>1), (x!=1)
@item       @b{>=}        @tab greater or equal @tab 3*(x>=1)
@end multitable


@sp 1
@noindent
@b{binary operators}
@multitable  {xxxxxxxxxxx} {xxxxxxxxxxxxxxxxxxxxxxxxxxxxx} {xxxxxxxxxxxxxxxxxxxxxxx}
@item         @b{^}        @tab  power           @tab x^2
@item         @b{*}        @tab  multiplication  @tab 2*x
@item         @b{/}        @tab  division        @tab x/2
@item         @b{+}        @tab  addition        @tab 2+x
@item         @b{-}        @tab  subtraction     @tab 2-x
@end multitable

@sp 2
Several examples for string-defined 'distributions':

@smallexample 

exp(-x^2/2)/(sqrt(2*pi))

exp(-sqrt(1+x^2)+x/2)

(x<1)*(x>0)*sin(x)
 

@end smallexample

@sp 2
@noindent
@b{Remarks:}
@itemize @bullet
@item
In case of doubt with the predendence please use parentheses.

@item
The figure 'e' is used both as a system constant written after an operator and as a seperator between mantisse and exponent for the basis 10. 

@item
Relation operators are very useful to define functions on a limited support and to define piecewise functions. 


@item
Of course it ist better to use for standard distributions build in functions.

@item
Especially with relation operators it is very easy to define functions in pieces :
So you can define the piecewise linear function defined throug the points (-1,0),(0,1),(1,0) with 
@example
(x>-1)*(x<0)*(1+x) + (x>=0)*(x<1)*(1-x) .
@end example 

@item 
It is better to truncate the support as to define a density function which is 0 outside of a special area.

@sp 2

*/








