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

If [ Environment["srcdir"] != $Failed, 
	SrcDir = Environment["srcdir"],
(* Else *)
        SrcDir =  "./" ];

(* name of C file for running tests *)
RunFileName = "t_StdDistr.c";

(* sample size for tests *)
RunSampleSize = 100;

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

(* --- Extreme Value II Distribution --------------------------------------- *)

ExtremeValueIIDistribution/: Domain[ExtremeValueIIDistribution[k_,zeta_,theta_]] := 
	Interval[{zeta, Infinity}];

ExtremeValueIIDistribution/: PDF[ExtremeValueIIDistribution[k_,zeta_,theta_], x_] :=
	If [ x > zeta,
		Exp[ -((x-zeta)/theta)^(-k)] * ((x-zeta)/theta)^(-k-1) * k/theta,
	(* else *)
		0
	];

ExtremeValueIIDistribution/: CDF[ExtremeValueIIDistribution[k_,zeta_,theta_], x_] :=
	If [ x > zeta,
		Exp[ -((x-zeta)/theta)^(-k)],
	(* else *)
		0
	];

ExtremeValueIIDistribution/: Random[ExtremeValueIIDistribution[k_,zeta_,theta_]] :=
	zeta + theta * Exp[ -Log[ -Log[Random[]] ]/k ];

(* --- Laplace Distribution ------------------------------------------------ *)

(* Mathematica's implementation of the Laplace Distribution is wrong!!       *)
(* We have to fix it!                                                        *)

Unprotect[LaplaceDistribution];
Clear[LaplaceDistribution];

LaplaceDistribution/: Domain[LaplaceDistribution[___]] := 
	Interval[{-Infinity, Infinity}];

LaplaceDistribution/: PDF[LaplaceDistribution[m_,s_], x_] :=
	If [ x-m < 0,
		Exp[(x-m)/s] / (2 s),
	(* Else *)
		Exp[(-x+m)/s] / (2 s) ];

LaplaceDistribution/: CDF[LaplaceDistribution[m_,s_], x_] :=
	If [ x-m < 0,
		Exp[(-m + x)/s]/2,
	(* Else *)
		1 - Exp[(m - x)/s]/2 ];

LaplaceDistribution/: Random[LaplaceDistribution[m_,s_]] := 
	m + If[ Random[]>0.5, 1,-1 ] * s * Log[Random[]];

(* --- Lomax Distribution -------------------------------------------------- *)

LomaxDistribution/: Domain[LomaxDistribution[___]] := 
	Interval[{0,1}];

LomaxDistribution/: PDF[LomaxDistribution[a_,C_], x_] :=
	(x+C)^(-(a+1)) * a * C^a;

LomaxDistribution/: CDF[LomaxDistribution[a_,C_], x_] :=
	1 - C^a/(C + x)^a;

LomaxDistribution/: Random[LomaxDistribution[a_,C_]] :=
	-C + (-(C^a/(-1 + Random[])))^(1/a);

(* --- Powerexponential Distribution --------------------------------------- *)

PowerexponentialDistribution/: Domain[PowerexponentialDistribution[___]] := 
	Interval[{-Infinity, Infinity}];

PowerexponentialDistribution/: PDF[PowerexponentialDistribution[r_], x_] :=
	If [ x>0, 
		Exp[-x^r]/(2 Gamma[1+1/r]),
	(* Else *)
		Exp[-(-x)^r]/(2 Gamma[1+1/r]) ];

PowerexponentialDistribution/: CDF[PowerexponentialDistribution[r_], x_] := 
	If [ x<0,
		Gamma[1/r, (-x)^r] / (2*r*Gamma[1 + 1/r]),
	(* Else *)
		1/2 + (Gamma[1 + 1/r] - Gamma[1/r, x^r]/r) / (2 Gamma[1+1/r]) ];

PowerexponentialDistribution/: Random[PowerexponentialDistribution[r_]] := 
	(* Warning!! This not a generator for the Powerexponential Distribution!!       *)
	(* FindRoot[ PDF[PowerexponentialDistribution[r],x] == Random[],{x,0.3}][[1,2]] *)
	(*  is correct but very slow. So we use a Laplace random variate instead !!     *)
	If[ Random[] < 0.5, 1, -1] * If[ Random[] < Max[0.3,(1-1/r)], Random[], 1-Log[Random[]]/r ];

(* --- Triangular Distribution --------------------------------------------- *)

TriangularDistribution/: Domain[TriangularDistribution[___]] := 
	Interval[{0,1}];

TriangularDistribution/: PDF[TriangularDistribution[H_], x_] :=
	Which [ x > 0 && x <= H,
			2*x / H,
		x > H && x < 1,
			2*(1-x) / (1-H),
		True,
			0
	];

TriangularDistribution/: CDF[TriangularDistribution[H_], x_] := 
	Which [ x <= 0,
			0,
		x > 0 && x <= H,
			x*x / H,
		x > H && x < 1,
			(H + x * (x-2))/(H-1),
		True,
			1
	];

TriangularDistribution/: Random[TriangularDistribution[H_]] := 
	(* Warning!! This not a generator for the Triangular Distribution!!       *)
	(* We simply use uniform distribution!!                                   *)
	Random[];

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
		WriteString[stream, CForm[ N[ fparams[[i]] ]]," "],
		{i,1,Length[fparams]} ];

	(* argument *)
	WriteString[stream, CForm[x]," "];

	(* CDF *)
	WriteString[stream, CForm[ N[ CDF[ Apply[distr,fparams], x ]]]," "];
	
	(* PDF *)
	WriteString[stream, CForm[ N[ PDF[ Apply[distr,fparams], x ]]]," "];

	(* derivative of PDF *)
	Switch [ dtype,
		isCONT,
			Fx = N[ 0. + D[ PDF[ Apply[distr,fparams], t ], t] /. t->x],
		isDISCR, (* discrete distribution: no derivative *)
			Fx = 0,
		_, (* default: should not happen *)
			Fx = 0
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
	(*		    (use isDISCR as 3rd entry in list for            *)
	(*		    discrete parameters)                             *)
	(*	size    ... size of sample                                   *)
	(*	distribution ... (optional)				     *)
	
	(* local variables *)
	{stream, datafilename, distrstring, i, j, nfparams, fparams, x},

	(* compose a file name for output *)
	datafilename = "../" <> SrcDir <> "tests/t_distr_" <> dname <> ".data";
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
			Table[ 
				If [ Length[fparbd[[j]]]>2 && fparbd[[j,3]] == isDISCR,
					Round[ fparbd[[j,1]] + Random[] * (fparbd[[j,2]]-fparbd[[j,1]]) ],
				(* else: isCONT *)
					N[ fparbd[[j,1]] + Random[] * (fparbd[[j,2]]-fparbd[[j,1]]) ]
				],
			{j,nfparams} ];	

		(* make an argument at random (use importance sampling) *)
		x = Random[ Apply[distr,fparams] ];
		If [ dtype == isCONT, x = N[x] ];

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
	WriteString[runfile,"   if (test_cdf_pdf( TESTLOG, distr,\"",datafilename,"\" ) == 0)\n"];
        WriteString[runfile,"      n_failed++;\n"];
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
#include \"test_StdDistr.h\"                                                   \n\
/*---------------------------------------------------------------------------*/\n\
                                                                               \n\
int main()                                                                     \n\
{                                                                              \n\
   double fparams[UNUR_DISTR_MAXPARAMS];  /* parameters of distribution      */\n\
   UNUR_DISTR *distr;                     /* distribution object             */\n\
   FILE *UNURANLOG;                       /* unuran log file                 */\n\
   FILE *TESTLOG;                         /* test log file                   */\n\
   int n_failed = 0;                      /* number of failed tests          */\n\
                                                                               \n\
   /* open log file for unuran and set output stream for unuran messages */    \n\
   UNURANLOG = fopen( \"t_Mathematica_unuran.log\",\"w\" );                    \n\
   if (UNURANLOG == NULL) exit (-1);                                           \n\
   unur_set_stream( UNURANLOG );                                               \n\
                                                                               \n\
   /* open log file for testing */                                             \n\
   TESTLOG = fopen( \"t_Mathematica_test.log\",\"w\" );                        \n\
   if (TESTLOG == NULL) exit (-1);                                             \n\
                                                                               \n\
   /* write header into log file */                                            \n\
   {                                                                           \n\
      time_t started;                                                          \n\
      fprintf(TESTLOG,\"\\nUNURAN - Universal Non-Uniform RANdom number generator\\n\\n\");\n\
      if (time( &started ) != -1)                                              \n\
         fprintf(TESTLOG,\"%s\",ctime(&started));                              \n\
      fprintf(TESTLOG,\"\\n======================================================\\n\\n\");\n\
      fprintf(TESTLOG,\"(Search for string \\\"data\\\" to find new section.)\\n\\n\");\n\
   }                                                                           \n\
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
   fclose(TESTLOG);                                                            \n\
                                                                               \n\
   if (n_failed > 0)                                                           \n\
      exit (-1);                                                               \n\
   else                                                                        \n\
      exit (0);                                                                \n\
}                                                                              \n"
	];

	(* close stream *)
	Close[stream];

]; (* end of UnurTestRunClose[] *)

(* === Start of Main ======================================================= *)

(* Open C file *)
runfile = UnurTestRunOpen[ RunFileName ];

(* --- List of Continuous Distributions ------------------------------------ *)

(* Beta *)
fparams = {{1,10}, {1,100}};
UnurTestDistrResultFile["beta", isCONT, runfile, fparams, RunSampleSize];

(* Cauchy *)
fparams = {{-100,100},{1/100,100}};
UnurTestDistrResultFile["cauchy", isCONT, runfile, fparams, RunSampleSize];

(* Chi *)
fparams = {{1,100}};
UnurTestDistrResultFile["chi", isCONT, runfile, fparams, RunSampleSize];

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

(* ExtremeValue II *)
fparams = {{1/100,100}, {-100,100}, {1/100,100}};
UnurTestDistrResultFile["extremeII",isCONT, runfile, fparams, RunSampleSize, ExtremeValueIIDistribution];

(* Gamma *)
fparams = {{1/2,10}, {1/100,100}};
UnurTestDistrResultFile["gamma", isCONT, runfile, fparams, RunSampleSize];

(* Laplace *)
fparams = {{-100,100}, {1/100,100}};
UnurTestDistrResultFile["laplace", isCONT, runfile, fparams, RunSampleSize];

(* Lomax *)
fparams = {{1/100,100}, {1/100,100}};
UnurTestDistrResultFile["lomax", isCONT, runfile, fparams, RunSampleSize];

(* Logistic *)
fparams = {{-100,100},{1/100,100}};
UnurTestDistrResultFile["logistic", isCONT, runfile, fparams, RunSampleSize];

(* Normal *)
fparams = {{-100,100}, {1/100,100}};
UnurTestDistrResultFile["normal", isCONT, runfile, fparams, RunSampleSize];

(* Pareto *)
fparams = {{1/100,100},{1/100,100}};
UnurTestDistrResultFile["pareto", isCONT, runfile, fparams, RunSampleSize];

(* Powerexponential *)
fparams = {{1/1000,3}};
UnurTestDistrResultFile["powerexponential", isCONT, runfile, fparams, RunSampleSize];

(* Rayleigh *)
fparams = {{1/100,100}};
UnurTestDistrResultFile["rayleigh", isCONT, runfile, fparams, RunSampleSize];

(* Student *)
fparams = {{1/100,100}};
UnurTestDistrResultFile["student", isCONT, runfile, fparams, RunSampleSize, StudentTDistribution];

(* Triangular *)
fparams = {{0,1}};
UnurTestDistrResultFile["triangular", isCONT, runfile, fparams, RunSampleSize];

(* Uniform *)
fparams = {{-100,1}, {1001/1000,100}};
UnurTestDistrResultFile["uniform", isCONT, runfile, fparams, RunSampleSize];

(* Weibull *)
fparams = {{1/2,10},{1/100,100}};
UnurTestDistrResultFile["weibull", isCONT, runfile, fparams, RunSampleSize];

(* --- List of Discrete Distributions -------------------------------------- *)

(* Binomial *)
fparams = {{2,1000,isDISCR}, {1/1000,999/1000}};
UnurTestDistrResultFile["binomial", is_DISCR, runfile, fparams, RunSampleSize];

(* Geometric *)
fparams = {{1/1000,999/1000}};
UnurTestDistrResultFile["geometric", is_DISCR, runfile, fparams, RunSampleSize];

(* Hypergeometric; order of parameters is different *)
hgd[N_,M_,n_] = HypergeometricDistribution[n,M,N];
fparams = {{100,1000,isDISCR}, {1,99,isDISCR}, {1,99,isDISCR}};
UnurTestDistrResultFile["hypergeometric", is_DISCR, runfile, fparams, RunSampleSize, hgd];

(* Logarithmic *)
fparams = {{1/1000,999/1000}};
UnurTestDistrResultFile["logarithmic", is_DISCR, runfile, fparams, RunSampleSize, LogSeriesDistribution];

(* NegativeBinomial; order of parameters is different *)
nb[p_,n_] = NegativeBinomialDistribution[n,p];
fparams = {{1/1000,999/1000}, {1,100,isDISCR}};
UnurTestDistrResultFile["negativebinomial", is_DISCR, runfile, fparams, RunSampleSize, nb];

(* Poisson *)
fparams = {{1/100,100}};
UnurTestDistrResultFile["poisson", is_DISCR, runfile, fparams, RunSampleSize];


(* --- Done ---------------------------------------------------------------- *)

(* close C file *)
UnurTestRunClose[runfile];

(* === Exit ================================================================ *)

Quit[]
