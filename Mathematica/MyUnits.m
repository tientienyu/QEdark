(* ::Package:: *)

BeginPackage["MyUnits`"]


(* ::Subsection:: *)
(*Constants*)


alpha=1/137.035999139;
me=0.5109989 "MeV";
mproton=938.2725 "MeV";
mneutron=939.5654 "MeV";
mmuon=105.6583745 "MeV";
mtau=1776.82 "MeV";
Gfermi=1.1664 10^-5 ("GeV")^-2;
sin\[Theta]Wsq=0.2397;
mPl=1.221 10^19 "GeV";(*Planck mass*)
MPl=2.435 10^18 "GeV";(*reduced Planck mass*)
HubbleConstant=72"km" ("s")^-1 ("Mpc")^-1;
hc=1240 "nm" "eV"; (*deBroglie wavelength*)
c=299792458 "m"/"s"; (*speed of light*)
NAvogadro=6.02214086 10^23 ("Mole")^-1; (* Avogadro's constant*)
rbohr=0.52917721067*10^-10 "m"; (* Bohr radius *)
hbar = 1.0545718 10^-34 ("m")^2 "kg"/"s"; (* hbar *)


(* ::Subsection:: *)
(*Unit Conversion*)


toGeV={
(*energy*)
"eV"->10^-9 "GeV","keV"->10^-6 "GeV","MeV"->10^-3 "GeV","TeV"->10^3 "GeV","PeV"->10^6 "GeV",
"J"->6.2415 10^9 "GeV",
"erg"->624.150974"GeV",
(*mass*)
"g"->5.609589 10^23 "GeV", "kg"->5.609589 10^26 "GeV", "tonne"->5.609589 10^29 "GeV","amu"->0.931494088"GeV",
"Msolar"->1.1157 10^57 "GeV",
(*distance*)
"fm"->5.067731("GeV")^-1,"nm"->5.067731 10^6 ("GeV")^-1,"mm"->5.067731 10^12 ("GeV")^-1,"cm"->5.067731 10^13 ("GeV")^-1,"m"->5.067731 10^15 ("GeV")^-1,"km"->5.067731 10^18 ("GeV")^-1,
"pc"->1.5637 10^32 ("GeV")^-1,"kpc"->1.5637 10^35 ("GeV")^-1,"Mpc"->1.5637 10^38 ("GeV")^-1,"AU"->7.581214085918`*^26/("GeV"),
"Rsolar"->3.522073045`*^24("GeV")^-1,
    "angstrom"->5.067731 10^5("GeV")^-1,
(*time*)
"ps"->1.519268 10^12 ("GeV")^-1,"ns"->1.519268 10^15 ("GeV")^-1,"\[Mu]s"->1.519268 10^18 ("GeV")^-1,"ms"->1.519268 10^21 ("GeV")^-1,"s"->1.519268 10^24 ("GeV")^-1, "day"->(24 60 60 1.519268 10^24) ("GeV")^-1, "year"->(365.25 24 60 60 1.519268 10^24)("GeV")^-1, "yr"->(365.25 24 60 60 1.519268 10^24)("GeV")^-1, "Myr"->(365.25 24 60 60 1.519268 10^30)("GeV")^-1, "Gyr"->(365.25 24 60 60 1.519268 10^33)("GeV")^-1,
(*cross-section*)
"b"->2.5681897 10^3 ("GeV")^-2,"mb"->2.5681897 ("GeV")^-2,"pb"->2.5681897 10^-9 ("GeV")^-2,"fb"->2.5681897 10^-12 ("GeV")^-2,
(*temperature*)
"K"->8.6173 10^-14 "GeV",
(* bohr *)
"bohr"->268173 "GeV",
(*magnetic field conversion*)
(*CAUTION: depends on system of units*)
(*here using canonical particle physics units*)
(*i.e. energy density Subscript[U, B]=1/2 B^2 ;  electron gyroradius Subscript[r, g]=p/(e B)= p/(Sqrt[4 \[Pi] \[Alpha]]B)*)
"nG"->1.95 10^-29 ("GeV")^2,"\[Mu]G"->1.95 10^-26 ("GeV")^2,"mG"->1.95 10^-23 ("GeV")^2,"G"->1.95 10^-20 ("GeV")^2,"T"->1.95 10^-16 ("GeV")^2
};


(*
take expression of form:
  number "units"
and return:
  new_number "new_units"
*)
Units[units_][expression_]:=Module[{GeV},
Simplify[(expression/units /.toGeV)/."GeV"-> GeV,Assumptions->GeV>0]units/.GeV->"GeV"]


(*
take expression of form:
  number "units"
and return:
  number
*)
RemoveUnits[expression_]:=If[NumberQ[expression],expression,expression[[1]]]


EndPackage[]


(* ::Subsection:: *)
(*Examples*)


(*use rule /.toGeV to convert all units into GeV^#*)
1"kg"/.toGeV
10"m"/.toGeV


(*use function Units[newunits][x] to convert x to the desired units*)
(*NB this is noticeably slower than using /.toGeV when used at high frequency*)
(* quicker to use /.toGeV//Units["year"^-1] *)
(0.3"GeV")/("cm")^3 1/(10"MeV") 10^-40 ("cm")^2//Units[("year")^-1]


(*use function RemoveUnits[x] to chop the units off x*)
10"m"//Units["s"]
10"m"//Units["s"]//RemoveUnits;
 



