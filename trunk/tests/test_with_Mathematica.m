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

(* === Set Constants ========================================================*)

(* name of C file for running tests *)
RunFileName = "run_test_with_Mathematica.c";

(* sample size for tests *)
RunSampleSize = 50;   
RunSampleSizeSlow = 10;   (* for slow routines *)

(* constants *)
isCONT = 1;
isDISCR = 2;

(* === Load statistics packages =============================================*)
<<Statistics`ContinuousDistributions`
<<Statistics`DiscreteDistributions`

(* === Declare additional distributions ==================================== *)

(* --- Defaults (to avoid error with C program) ---------------------------- *)

CDF[___] := 0.;
PDF[___] := 0.;

(* --- Powerexponential Distribution --------------------------------------- *)

PowerexponentialDistribution/: Domain[PowerexponentialDistribution[___]] := Interval[{-Infinity, Infinity}]

PowerexponentialDistribution/: PDF[PowerexponentialDistribution[r_], x_] :=
	If [ x>0, 
		Exp[-x^r]/(2 Gamma[1+1/r]),
	(* Else *)
		Exp[-(-x)^r]/(2 Gamma[1+1/r]) ];

PowerexponentialDistribution/: CDF[PowerexponentialDistribution[r_], x_] := 
	If [ x<0,
		Integrate[ Exp[-(-t)^r], {t,-Infinity,x}] / (2 Gamma[1+1/r]),
	(* Else *)
		1/2 + Integrate[Exp[-(t)^r], {t,0,x}] / (2 Gamma[1+1/r]) ];

PowerexponentialDistribution/: Random[PowerexponentialDistribution[r_]] := 
	(* Warning!! This not a generator for the Powerexponential Distribution!!       *)
	(* FindRoot[ PDF[PowerexponentialDistribution[r],x] == Random[],{x,0.3}][[1,2]] *)
	(*  is correct but very slow. So we use a Laplace random variate instead !!     *)
	If[ Random[] < 0.5, 1, -1] * If[ Random[] < Max[0.3,(1-1/r)], Random[], 1-Log[Random[]]/r ];

(* --- Laplace Distribution ------------------------------------------------ *)

(* Mathematica's implementation of the Laplace Distribution is wrong!!       *)
(* We have to fix it!                                                        *)

Unprotect[LaplaceDistribution];
Clear[LaplaceDistribution];

LaplaceDistribution/: Domain[LaplaceDistribution[___]] := Interval[{-Infinity, Infinity}]

LaplaceDistribution/: PDF[LaplaceDistribution[m_,s_], x_] :=
	If [ x-m < 0,
		Exp[(x-m)/s] / (2 s),
	(* Else *)
		Exp[(-x+m)/s] / (2 s) ];

LaplaceDistribution/: CDF[LaplaceDistribution[m_,s_], x_] :=
	If [ x-m < 0,
		Integrate[ Exp[(t-m)/s] / (2 s), {t,-Infinity,x}],
	(* Else *)
		1/2 + Integrate[ Exp[(-t+m)/s] / (2 s), {t,m,x}] ];

LaplaceDistribution/: Random[LaplaceDistribution[m_,s_]] := 
	m + If[ Random[]>0.5, 1,-1 ] * s * Log[Random[]];

(* === End of Distributions ================================================ *)

(* === Format and write output ============================================= *)

(* --- Write result for parameter/argument combination into file ----------- *)

UnurTestDistrResultLine[stream_, distr_, dtype_, fparams__, x_] := Module [
	(*	stream  ... output stream                                    *)
	(*   	distr   ... distribution                                     *)
	(*	dtype   ... type of distribution (CONT|DISCR|...)            *)
	(*	fparams ... parameters for distribution                      *)
	(*	x       ... argument where `distr' should be evaluated       *)

	(* local variables *)
	{i,Fx},

	(* number of parameters *)
	WriteString[stream, Length[fparams]," "];

	(* list of parameters *)
	Table[
		WriteString[stream, CForm[ N[ fparams[[i]]+0., $MachinePrecision+1]]," "],
		{i,1,Length[fparams]} ];

	(* argument *)
	WriteString[stream, CForm[ N[x+0.,$MachinePrecision+1] ]," "];

	(* CDF *)
	WriteString[stream, CForm[ N[ CDF[ Apply[distr,fparams], x ]+0., $MachinePrecision+1]]," "];
	
	(* PDF *)
	WriteString[stream, CForm[ N[ PDF[ Apply[distr,fparams], x+0. ]+0., $MachinePrecision+1]]," "];

	(* derivative of PDF *)
	Switch [ dtype,
		isCONT,
			Fx = N[ 0. + D[ PDF[ Apply[distr,fparams], t ], t] /. t->x, $MachinePrecision+1],
		isDISCR, (* discrete distribution: no derivative *)
			Fx = 0.,
		_, (* default: should not happen *)
			Fx = 0.
	];
	WriteString[stream, CForm[ Fx ]];

	(* end *)
	WriteString[stream, "\n"];

]; (* end of UnurTestDistrResultLine[] *)

(* --- make file with results for given distribution ----------------------- *)

UnurTestDistrResultFile[dname_, dtype_, runfile_, fparbd__, size_, distribution___] := Module [ 
	(*	dname   ... UNURAN name of distribution                      *)
	(*	dtype   ... type of distribution (CONT|DISCR|...)            *)
	(*	runfile ... output stream for C file running tests           *)
	(*	fparbd  ... bounds for parameters for distribution           *)
	(*	size    ... size of sample                                   *)
	(*	distribution ... (optional)				     *)
	
	(* local variables *)
	{stream, datafilename, distrstring, i, j, nfparams, fparams, x},

	(* compose a file name for output *)
	datafilename = "t_distr_" <> dname <> ".data";
	Print[ datafilename ];

	(* open output stream *)
	stream = OpenWrite[datafilename];

	(* distribution *)
	If [ Length[{distribution}] > 0,
		(* distribution given as third argument *)
		distr = distribution,
	(* else *)
		(* get distribution from dname *)
		distrstring = ToUpperCase[ StringTake[dname,1] ]
				<> StringDrop[dname,1] 
				<> "Distribution";
		distr = ToExpression[distrstring];
	];

	(* number of parameters for distribution *)
	nfparams = Length[fparbd];

	Do[
		(* make parameters for distribution at random *)
		fparams = 
			Table[ Block [ {prec}, 
				(* read precision (3rd optional argument) *)
				prec = If [ Length[fparbd[[j]]]==2, 1000, fparbd[[j,3]] ];
				(* choose parameter at random *)
				fparbd[[j,1]] + Round[(fparbd[[j,2]]-fparbd[[j,1]])*prec*Random[]]/prec ],
			{j,nfparams} ];	

		(* make an argument at random (use importance sampling) *)
		x = Round[1000 * Random[ Apply[distr,fparams] ]] / 1000;

		(* print line with this parameters *)
		UnurTestDistrResultLine[stream,distr,dtype,fparams,x],

	{i,size} ];

	(* close output stream *)
	Close[stream];

	(* print test line into C file running tests *)

	(* list of parameters. we use upper bound *)
	For [ i=0, i<nfparams, i++,
		WriteString[runfile,"   fparams[",i,"] = ",N[fparbd[[i+1,2]]],";\n"];
	];
	WriteString[runfile,"   distr = unur_distr_",dname,"(fparams,",nfparams,");\n"];
	WriteString[runfile,"   test_cdf_pdf( distr,\"",datafilename,"\", 1. );\n"];
        WriteString[runfile,"   unur_distr_free(distr);\n\n"];

]; (* end of UnurTestDistrResultFile[] *)

(* --- Open C file for output----------------------------------------------- *)

UnurTestRunOpen[filename_] := Module [
	(*	filename ... name of output file                             *)

	(* local variables *)
	{stream},

	(* open file for writing *)
	stream = OpenWrite[filename];

	(* write file preamble *)
	WriteString[stream,"\
/***************************************************************************** \n\
 *                                                                           * \n\
 *          UNURAN -- Universal Non-Uniform Random number generator          * \n\
 *                                                                           * \n\
 *****************************************************************************/\n\
                                                                               \n\
/*---------------------------------------------------------------------------*/\n\
/*  #define DEBUG 1 */                                                         \n\
/*---------------------------------------------------------------------------*/\n\
#include \"test_with_Mathematica.h\"                                           \n\
/*---------------------------------------------------------------------------*/\n\
                                                                               \n\
int main()                                                                     \n\
{                                                                              \n\
   double fparams[UNUR_DISTR_MAXPARAMS];  /* parameters of distribution      */\n\
   UNUR_DISTR *distr;                     /* distribution object             */\n\
   FILE *UNURANLOG;                       /* unuran log file                 */\n\
                                                                               \n\
   /* open log file for unuran and set output stream for unuran messages */    \n\
   UNURANLOG = fopen( \"t_Mathematica_unuran.log\",\"w\" );                    \n\
   if (UNURANLOG == NULL) exit (-1);                                           \n\
   unur_set_stream( UNURANLOG );                                               \n\
                                                                               \n"
	];

	(* return stream handler *)
	Return[ stream ];

]; (* end of UnurTestRunOpen[] *)

(* --- Close output file -------------------------------------------------- *)

UnurTestRunClose[stream_] := Module [
	(*	stream ... output stream                                    *)

	(* local variables: none *)
	{},

	(* write file postamble *)
	WriteString[stream,"\
                                                                               \n\
   /* close log file */                                                        \n\
   fclose(UNURANLOG);                                                          \n\
                                                                               \n\
   exit (0);                                                                   \n\
}                                                                              \n"
	];

	(* close stream *)
	Close[stream];

]; (* end of UnurTestRunClose[] *)

(* === Start of Main ======================================================= *)

(* Open C file *)
runfile = UnurTestRunOpen[ RunFileName ];

(* --- List of Continuous Distributions ------------------------------------ *)
(*
(* Beta *)
fparams = {{1,10}, {1,100}};
UnurTestDistrResultFile["beta", isCONT, runfile, fparams, RunSampleSize];

(* Cauchy *)
fparams = {{-100,100},{1/100,100}};
UnurTestDistrResultFile["cauchy", isCONT, runfile, fparams, RunSampleSize];
*)
(* Chi *)
fparams = {{1,100}};
UnurTestDistrResultFile["chi", isCONT, runfile, fparams, RunSampleSize];
(*
(* Chisquare *)
fparams = {{1,100}};
UnurTestDistrResultFile["chisquare", isCONT, runfile, fparams, RunSampleSize, ChiSquareDistribution];

(* Exponential -- parameters differ *)
fparams = {{1/100,100}};
ed[mu_] = ExponentialDistribution[1/mu];
UnurTestDistrResultFile["exponential", isCONT, runfile, fparams, RunSampleSize, ed];

(* ExtremeValue I *)
fparams = {{-100,100},{1/100,100}};
UnurTestDistrResultFile["extremeI",isCONT, runfile, fparams, RunSampleSize, ExtremeValueDistribution];

(* Gamma *)
fparams = {{1/2,10},{1/100,100}};
UnurTestDistrResultFile["gamma", isCONT, runfile, fparams, RunSampleSize];

(* Laplace *)
fparams = {{-100,100}, {1/100,100}};
UnurTestDistrResultFile["laplace", isCONT, runfile, fparams, RunSampleSize];

(* Logistic -- reverse order of parameters *)
fparams = {{1/100,100},{-100,100}};
ld[beta_,mu_] = LogisticDistribution[mu,beta];
UnurTestDistrResultFile["logistic", isCONT, runfile, fparams, RunSampleSize, ld];

(* Normal *)
fparams = {{-100,100}, {1/100,100}};
UnurTestDistrResultFile["normal", isCONT, runfile, fparams, RunSampleSize];

(* Pareto *)
fparams = {{1/100,100},{1/100,100}};
UnurTestDistrResultFile["pareto", isCONT, runfile, fparams, RunSampleSize];

(* Powerexponential *)
fparams = {{1/1000,3}};
UnurTestDistrResultFile["powerexponential", isCONT, runfile, fparams, RunSampleSizeSlow];

(* Rayleigh *)
fparams = {{1/100,100}};
UnurTestDistrResultFile["rayleigh", isCONT, runfile, fparams, RunSampleSize];

(* Student *)
fparams = {{1/100,100}};
UnurTestDistrResultFile["student", isCONT, runfile, fparams, RunSampleSize, StudentTDistribution];

(* Weibull *)
fparams = {{1/2,10},{1/100,100}};
UnurTestDistrResultFile["weibull", isCONT, runfile, fparams, RunSampleSize];
*)

(* --- List of Discrete Distributions -------------------------------------- *)
(*
(* Binomial *)
fparams = {{2,1000,1(*as parameter is discrete*)},{1/1000,999/1000}};
UnurTestDistrResultFile["binomial", is_DISCR, runfile, fparams, RunSampleSize];

(* Geometric *)
fparams = {{1/1000, 999/1000}};
UnurTestDistrResultFile["geometric", is_DISCR, runfile, fparams, RunSampleSize];

(* Hypergeometric; order of parameters is different *)
hgd[N_,M_,n_] = HypergeometricDistribution[n,M,N];
fparams = {{100,1000,1},{1,99,1},{1,99,1}};
UnurTestDistrResultFile["hypergeometric", is_DISCR, runfile, fparams, RunSampleSize, hgd];

(* Logarithmic *)
fparams = {{1/1000, 999/1000}};
UnurTestDistrResultFile["logarithmic", is_DISCR, runfile, fparams, RunSampleSize, LogSeriesDistribution];

(* NegativeBinomial; order of parameters is different *)
nb[p_,n_] = NegativeBinomialDistribution[n,p];
fparams = {{1/1000,999/1000},{1,100,1(*as parameter is discrete*)}};
UnurTestDistrResultFile["negativebinomial", is_DISCR, runfile, fparams, RunSampleSize, nb];

(* Poisson *)
fparams = {{1/100,100}};
UnurTestDistrResultFile["poisson", is_DISCR, runfile, fparams, RunSampleSize];

*)

(* --- Done ---------------------------------------------------------------- *)

(* close C file *)
UnurTestRunClose[runfile];

(* === Exit ================================================================ *)

Quit[]







