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

   

(* Load statistic packages *)
<<Statistics`ContinuousDistributions`
<<Statistics`NormalDistribution`


(* format output string *)
UnurTestDistrResult[stream_, distr_, fparams__, x_] := Module[ {i},
	(*	stream  ... output stream
	   	distr   ... distribution
		fparams ... parameters for distribution
		x       ... argument where `distr' should be evaluated
	*)

	(* number of parameters *)
	WriteString[stream, Length[fparams]," "];

	(* list of parameters *)
	Table[
		WriteString[stream, CForm[ N[ fparams[[i]], 20]]," "],
		{i,1,Length[fparams]} ];

	(* argument *)
	WriteString[stream, CForm[ N[x,20] ]," "];

	(* CDF *)
	WriteString[stream, CForm[ N[ CDF[ Apply[distr,fparams], x ], 20]]," "];
	
	(* PDF *)
	WriteString[stream, CForm[ N[ PDF[ Apply[distr,fparams], x ], 20]]," "];
	
	(* derivative of PDF *)
	WriteString[stream, CForm[ N[ D[ PDF[ Apply[distr,fparams], t ], t] /. t->x, 20]]," "];
	


	(* end *)
	WriteString[stream, "\n"];
];


(* make a file with results *)
UnurTestDistrResultFile[dname_, distr_, fparbd__, size_] := Module[ 
	(*	dname  ... UNURAN name of distribution
	   	distr  ... distribution
		fparbd ... bounds for parameters for distribution
		size   ... size of sample
	*)
	{stream, filename, i, j, nfparams, fparams, x},

	(* compose a file name for outpit *)
	filename = "t_distr_" <> dname <> ".data";
	Print[ filename ];
	
	(* open output stream *)
	stream = OpenWrite[filename];

	(* number of parameters for distribution *)
	nfparams = Length[fparbd];

	Do[
		(* make parameters for distribution at random *)
		fparams = 
			Table[ fparbd[[j,1]] + Round[(fparbd[[j,2]]-fparbd[[j,1]])*1000*Random[]]/1000,
			{j,nfparams} ];	

		(* make an argument at random (use importance sampling) *)
		x = Round[1000 * Random[ Apply[distr,fparams] ]] / 1000;

		(* print line with this parameters *)
		UnurTestDistrResult[stream,distr,fparams,x],

	{i,size} ];

	(* close output stream *)
	Close[stream];
];

(* now make files for distributions *)

(* Beta *)
UnurTestDistrResultFile["beta",BetaDistribution,{{1.,10.},{1.,100.}},50];

(* Gamma *)
UnurTestDistrResultFile["gamma",GammaDistribution,{{0.5,10.},{0.01,100.}},50];

(* Normal *)
UnurTestDistrResultFile["normal",NormalDistribution,{{0.5,10.},{0.01,100.}},50];

Quit[]

