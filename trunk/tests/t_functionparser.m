

filename ="t_functionparser.data"; 
(* Eingeben der zu testenden Funktion im Format : 
      {Funktion_als_string",{kleinster_Funktionswert,maximaler_Funktionswert,\
Anzahl_Funktionswerte}       *)
Testsample = {{"3+5", {-2, 2, 5}},
      		{"3*x", {-2, 2, 5}},
      		{"exp[-4*X]", {-2, 2, 5}},
      		{"exp[-x^2]+Log[2,4]-Pi*Sin[x+x*2]", {-5, 2, 4}},
                      {"Sin[x]*3*Ln[x]",   {2, 4, 2}  },
      	      {"exp[x^2]*(cos[x]<1)", {-3, 8, 5}},
                      {"abs[x]-3*x", {-2, 2, 5}}};





 (* Definition der Konstanten *)
pi := Pi;
e := E;

(* Definition der parser - funktionen *)
exp[x_] := Exp[x];
ln[x_] := Log[x];
log[x_, y_] := Log[x, y];
sin[x_] := Sin[x];
cos[x_] := Cos[x];
tan[x_] := Tan[x];
sqrt[x_] := Sqrt[x];
abs[x_] := Abs[x];
sgn[x_] := Sgn[x];

exprlog[x_] := If[x > 1, 3, 4];

(* Umdefinition der Ableitungsfunktion analog der Bedeutung im Parser *)
(* 
  der_abs[x_] := If[x < 0, D[x, x
          (Unprotect[D];
            D[abs[x_], x_] := der_abs[x];
            Protect[D]);
        D[abs[x], x] /. x -> 2 *)

(* Ableitung > *)
(Unprotect[Derivative];
    Derivative[1, 0][Greater][x_, y_] := 0;
    Derivative[0, 1][Greater][x_, y_] := 0;
    Protect[Derivative]);

(* Ableitung < *)
(Unprotect[Derivative];
    Derivative[1, 0][Less][x_, y_] := 0;
    Derivative[0, 1][Less][x_, y_] := 0;
    Protect[Derivative]);

(* Ableitung 'abs' *)
(Unprotect[Derivative];
    Derivative[1][Abs][x_] := If[x < 0, -1, 1];
    Protect[Derivative]);

(*
  D[Greater[x, 4], x]
          D[x > 1, x]
          D[Sin[x]*(2*x > x), x] /. x -> 2 /. {True -> 1, False -> 0} // N
  *)


(* Berechnung und Ausgeben der Argumente, 
  Funktionswerte sowie der ersten Ableitung       *)

Do[  
  Testfunktion = Testsample[[i]];
  fstr = ToLowerCase[Testfunktion[[1]]];
  fstrmin = Testfunktion[[2]][[1]];
  fstrmax = Testfunktion[[2]][[2]];
  fstrstep = (fstrmax - fstrmin)/(Testfunktion[[2]][[3]] - 1);
  
  		
  (* Parse - 
            string in Datei schreiben; [, ] wird durch (, ) ersetzt sowie \
alle Funktionen klein geschrieben *)
  
  WriteString[filename, 
    OutputForm[
      StringJoin[ {"function=", 
          ToLowerCase[  StringReplace[ fstr, {"[" -> "(", "]" -> ")"}  ]  ] , 
          "\n"}  ] ]  ] ;
  
  
  (* Argument, Funktionswert und erste Ableitung in Datei schreiben *)
  Do[
     WriteString[filename,
        StringJoin[ToString[N[j]],
                                 "\t",
          		    
          ToString[ 
            CForm[N[ToExpression[fstr] /. { x -> j} /. {True -> 1, 
                    False -> 0}]]], "\t" ,
                                 
          ToString[
            CForm[N[  
                D[ToExpression [fstr] , x]    /. x -> j /. {True -> 1, 
                    False -> 0}]]]       , 
          "\n"                               ]];		
    , {j, fstrmin, fstrmax, fstrstep}  ];
  
  WriteString[filename, "\n"] ;
  
  (* naechste Funktion einlesen *)
  , {i, 1, Length[Testsample], 1}]

Quit[]

