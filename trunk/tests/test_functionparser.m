(*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
     $Id$ 
 *****************************************************************************
 *                                                                           *
 *  Create file with results of CDF, PDF and derivatives of PDF              *
 *  for various distributions and make C test file.                          *
 *                                                                           *
 *****************************************************************************)

(* === List of tests ========================================================*)

(* Eingeben der zu testenden Funktion im Format : 
      {"Funktion_als_string",{kleinster_Funktionswert,maximaler_Funktionswert,\
Anzahl_Funktionswerte}       *)
Testsample = {
         {"3.45+x", {-2, 2, 5}},                (* arithmetic operator *)
         {"x-4.673", {-2, 2, 5}},
         {"3*x", {-2, 2, 5}},
         {"x/5.8", {-2, 2, 5}},
	 {"x^12-3.345+5", {-2, 2, 5}},
      	 {"3*x", {-2, 2, 5}},
         {"2.894736*10^2", {-2, 2, 5}},

	 {"3*(x^5-x^4/(1.5-x))", {-2, 2, 5}},    (*  brackets  *)

        {"3*(x>1)" ,{-2,2,5}},
        {"3*(x<1)" ,{-2,2,5}},
        {"3*(x>=1)",{-2,2,5}},
        {"3*(x<=1)",{-2,2,5}},
        {"3*(x==1)",{-2,2,5}},
        {"3*(2==x)",{-2,2,5}},
        {"3*(x<>1)",{-2,2,5}},
        {"3*(2<>x)",{-2,2,5}},
 
        {"3*not[2<>x]",{-2,2,5}},
        {"3*(2<>x)and(x>2)",{-2,2,5}},
        {"3*(2<x)or(x<-1)",{-2,2,5}},
    
	{"exp[-4*X]", {-2, 2, 5}},               (*   system functions  *)
        
        {"ln[x]",{1,6,6}},
        {"log[3,x]",{1,6,6}},  
        {"sqrt[x]",{1,6,6}},
        {"sin[x]",{1,6,6}},
        {"mod[x,3]",{1,6,6}},
        {"sgn[x]",{-2,2,5}},
        {"sec[x]",{-2,2,5}},
      	{"exp[-x^2]+Log[2,4]-Pi*Sin[x+x*2]", {-5, 2, 4}},
        {"Sin[x]*3*Ln[x]",   {2, 4, 2}  },
      	{"exp[x^2]*(cos[x]<1)", {-3, 8, 5}},
	{"abs[x]-3*x", {-2, 2, 5}},

        {"exp[x^2]*(sin[x*cos[x^2-1]]+1)*((x-3*pi*x)<1)", {-3148,789,5}} };
  


(* === Set Constants ========================================================*)

(* name of datafile file for running tests *)
DataFile = "t_functionparser.data"; 


(* === Define cosntants and functions for function parser ===================*)

(* --- Constants ----------------------------------------------------------- *)

pi = Pi;
e  = E;

(* --- Functions ----------------------------------------------------------- *)

mod[x_,y_]  := Mod[x,y];
exp[x_]     := Exp[x];
ln[x_]      := Log[x];
log[x_, y_] := Log[x, y];
sin[x_]     := Sin[x];
cos[x_]     := Cos[x];
tan[x_]     := Tan[x];
sec[x_]     := 1/Cos[x];
sqrt[x_]    := Sqrt[x];
abs[x_]     := Abs[x];
sgn[x_]     := Sign[x];    

not[x_]     := Not[x];

(* === Define derivatives for these functions (according to function parser) *)

Unprotect[Derivative];

(* --- Relation Operators -------------------------------------------------- *)
(* Derivative[1][And][x_] := 0;
 Derivative[1][Or][x_] := 0; *)

Derivative[1][Not] [x_]  := 0;
Derivative[1][Sign][x_] := 0;

Derivative[1, 0][Unequal][x_, y_] := 0;
Derivative[0, 1][Unequal][x_, y_] := 0;

Derivative[1, 0][Greater][x_, y_] := 0;
Derivative[0, 1][Greater][x_, y_] := 0;

Derivative[1, 0][GreaterEqual][x_, y_] := 0;
Derivative[0, 1][GreaterEqual][x_, y_] := 0;

Derivative[1, 0][Less][x_, y_] := 0;
Derivative[0, 1][Less][x_, y_] := 0;

Derivative[1, 0][LessEqual][x_, y_] := 0;
Derivative[0, 1][LessEqual][x_, y_] := 0;

Derivative[1, 0][Mod][x_, y_] := 0;
Derivative[0, 1][Mod][x_, y_] := 0;


(* --- Functions ----------------------------------------------------------- *)

Derivative[1][Abs][x_] := sgn[x];

(* ------------------------------------------------------------------------- *)

Protect[Derivative];


(* === Write results for an expression into data file ====================== *)

UnurWriteData[expression_,points__] := Module [
	(*	expression ... function term                                 *)
	(*      points     ... list {x_min, x_max, number of points}         *)

	(* local variables                                                   *)
	{
	  funct,   (* string use to compute expression with Mathematica      *)
	  x,       (* argument for function                                  *)
	  xmin,    (* minimal value for x                                    *)
	  xmax,    (* maximal value for x                                    *)
	  xstep    (* step width                                             *)
	},

	(* convert to lower case letters                                     *)
	funct = ToLowerCase[expression];

	(* print function string into data file                              *)
	WriteString[DataFile,
		    "function=",
		    UnurTransformExpression[funct],
		    "\n"]; 

	(* get values for x                                                  *)
	xmin = points[[1]];
	xmax = points[[2]];
	xstep = (xmax - xmin) / (points[[3]]-1);
	xmax += xstep/2;

	(* print function and its derivative at all given points *)
	Do[ UnurWriteLine[x,funct],{x,xmin,xmax,xstep}];

	(* add blank line *)
	WriteString[DataFile,"\n"]; 

]; (* end of UnurWriteData[] *)


(* --- Transform expression into string for function parser ---------------- *)

UnurTransformExpression[expression_] := Module [
	(*	expression ... function term                                 *)

	(* local variables                                                   *)
	{fstr},

	(* Replace square brackets by parenthesis                            *)
	fstr = StringReplace[ expression, {"[" -> "(", "]" -> ")"}];

	(* return result *)
	Return[fstr];

]; (* end of UnurTransformExpression[] *)


(* --- Compute function and write line into date file       ---------------- *)

UnurWriteLine[xarg_,funct_] := Module [
	(*      xarg  ... argument                                           *)
	(*	funct ... function term                                      *)

	(* local variables                                                   *)
	{
	  xval, (* numerical value of argument x                             *)
	  fx,   (* value of function at x                                    *)
	  dfx,  (* value of derivative of function at x                      *)
	  fstr  (* functtion term                                            *)
        },

	(* argument *)
	xval = N[xarg];
	WriteString[DataFile, CForm[xval],"\t"];

	(*  replace '<>' -> '!='   *)
        fstr = StringReplace[funct, {"<>" -> "!="}];
 
        (*  replace 'and' -> '&&', 'or' -> '||'   *)
        fstr = StringReplace[fstr,{"and" -> "&&","or" -> "||"}]; 

					(* WriteString["stdout",fstr,"\n"]; *)
	(* function *)
	fx = N[ToExpression[fstr]] /. x -> xval /. {True -> 1, False -> 0};
	WriteString[DataFile, CForm[fx],"\t"];

	(* derivative *)
	dfx = N[ D[ ToExpression[fstr], x ]] /. x -> xval /. {True -> 1, False -> 0};
	WriteString[DataFile, CForm[dfx],"\n"];

]; (* end of UnurWriteLine[] *)

(* === Main ================================================================ *)

Do [
   fstr   = Testsample[[i]][[1]];
   points = Testsample[[i]][[2]];
   UnurWriteData[fstr, points],
   {i, 1, Length[Testsample]}
];

(* === Exit ================================================================ *)

Quit[]


