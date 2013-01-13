(* ::Package:: *)

BeginPackage["AvMLJ`"]
GetProbs::usage="S, hbar, DG, \[Lambda]o, kT, nmax, nmax,{\[Epsilon]\[Rule] 4, d\[Rule] 10, F\[Rule] 0.005} (optional)"
RateMLJ::usage="J_, \[Lambda]i_,h\[Omega]_, S_, dE_,kT";
TermMLJe::usage="[S_,h\[Omega]_, dE_, \[Lambda]i_,kT_,{FROM i_, TO n_}"
P::usage="[{A_,S}]"
Prefactor::usage="[J_, \[Lambda]i_,kT_]"
Begin["Private`"]
(*the Franck Condond factors*)
Clear[f, WaveFun];
WaveFun[x_, n_] := 1/Sqrt[2^n n!] Pi^(-1/4) Exp[-(x^2/2)]HermiteH[n,x];

TermsPolaron[\[Lambda]_, i_,j_] :=Module[{terms},
Exp[-(\[Lambda]^2/2)]Sum[(-\[Lambda])^k/k! \[Lambda]^(i-j+k)/(i-j+k)! Sqrt[j!i!]/(j-k)!,{k,0,j}]
]
f[a_,b_] := TermsPolaron[\[Lambda],a,b]

(*the normal MLJ expression from the vibrational ground state*)
tollMLJ=10^-3;
Clear[TermMLJ, SumMLJ,TermMLJe, SumMLJe];
TermMLJ[S_,h\[Omega]_, dE_, \[Lambda]i_,kT_,n_]:= S^n/n! Exp[-((dE+\[Lambda]i+n h\[Omega])^2/(4 \[Lambda]i kT))] Exp[-S]
SumMLJ[S_,h\[Omega]_, dE_, \[Lambda]i_,kT_]:= SumMLJ[S,h\[Omega], dE, \[Lambda]i,kT,0,TermMLJ[S,h\[Omega], dE, \[Lambda]i,kT,0],1]
SumMLJ[S_,h\[Omega]_, dE_, \[Lambda]i_,kT_, valsofar_,dv_,n_] := {valsofar+dv,n} /; (dv/(valsofar+dv) < tollMLJ && n > -Ceiling[dE/h\[Omega]])
SumMLJ[S_,h\[Omega]_, dE_, \[Lambda]i_,kT_, valsofar_,dv_,n_] := SumMLJ[S,h\[Omega], dE, \[Lambda]i,kT, valsofar+dv,TermMLJ[S,h\[Omega], dE, \[Lambda]i,kT,n],n+1] /; (dv/(valsofar+dv) >=  tollMLJ  || n > -Ceiling[dE/h\[Omega]])
\[CapitalGamma]MLJ[pref_, \[Lambda]i_,h\[Omega]_, S_, dE_,kT_] :=Module[{exp} ,
exp =  First[SumMLJ[S, h\[Omega], dE,\[Lambda]i,kT]];
pref*exp
]
RateMLJ[J_, \[Lambda]i_,h\[Omega]_, S_, dE_,kT_]:= \[CapitalGamma]MLJ[Prefactor[J, \[Lambda]i,kT], \[Lambda]i,h\[Omega],S,dE,kT]

(*the semi classical expression*)
Msm[J_, \[Lambda]_, \[CapitalDelta]G_, kT_]:=J^2/hbar Sqrt[Pi/(\[Lambda] kT)] Exp[-((\[Lambda] +\[CapitalDelta]G)^2/(4 \[Lambda] kT))]
(*the "normal" MLJ expression from the vibrational excited state i*)
TermMLJe[S_,h\[Omega]_, dE_, \[Lambda]i_,kT_,{i_,n_}]:= ((f[i,n])/.(\[Lambda]-> Sqrt[S] ))^2 Exp[-((dE+\[Lambda]i+(n-i) h\[Omega])^2/(4 \[Lambda]i kT))]

SumMLJe[S_,h\[Omega]_, dE_, \[Lambda]i_,kT_, i_]:= SumMLJe[S,h\[Omega], dE, \[Lambda]i,kT,i,0.,TermMLJe[S,h\[Omega], dE, \[Lambda]i,kT,{i,0}],1]

SumMLJe[S_,h\[Omega]_, dE_, \[Lambda]i_,kT_, i_ ,valsofar_,dv_,n_] := valsofar+dv /; (dv/(valsofar+dv) < tollMLJ  && n> i)

SumMLJe[S_,h\[Omega]_, dE_, \[Lambda]i_,kT_, i_ ,valsofar_,dv_,n_] := SumMLJe[S,h\[Omega], dE, \[Lambda]i,kT, i, valsofar+dv,TermMLJe[S,h\[Omega], dE, \[Lambda]i,kT,{i,n}],n+1] /; (dv/(valsofar+dv) >=  tollMLJ || n <=  i)

\[CapitalGamma]MLJe[pref_, \[Lambda]i_,h\[Omega]_, S_, dE_,kT_,i_] :=Module[{exp} ,
exp =  SumMLJe[S, h\[Omega], dE,\[Lambda]i,kT,i];
pref*exp
]
hbar = 6.57 10^-16 (*in eV s^-1*)
Prefactor[J_, \[Lambda]i_,kT_]:= J^2/hbar Sqrt[Pi/(\[Lambda]i kT)]

P[{A_,S_},{\[Mu]_,\[Nu]_}]:= 2^-(\[Mu]+\[Nu]) KroneckerDelta[A+S, \[Mu]+\[Nu]] ( Sum[Binomial[\[Mu],m] Binomial[\[Nu],n] (-1)^(\[Nu]-n) Sqrt[((m+n)!(\[Mu]+\[Nu]-m-n)!)/(\[Mu]!\[Nu]!)]KroneckerDelta[S, m+n],{m,0,\[Mu]}, {n,0,\[Nu]}])^2
P[{A_,S_}] := P[{A,S},{#,A+S-#}]&/@Range[0, A+S]
Clear [GetProbs];
GetProbs[S_, \[HBar]\[Omega]_, \[CapitalDelta]G_,\[Lambda]o_, kT_, nmax_] := GetProbs[S, \[HBar]\[Omega], \[CapitalDelta]G,\[Lambda]o, kT, nmax,{\[Epsilon]-> 4, d-> 10, F-> 0.005}] 

GetProbs[S_, \[HBar]\[Omega]_, \[CapitalDelta]G_,\[Lambda]o_, kT_, nmax_, NRGOptions_]:=Module[{All\[CapitalGamma],\[CapitalDelta]GCorr},
	\[CapitalDelta]GCorr = \[CapitalDelta]G - (27.21 *0.529)/(\[Epsilon] d)- F d /. NRGOptions;
	All\[CapitalGamma] = TermMLJ[S,\[HBar]\[Omega],\[CapitalDelta]GCorr, \[Lambda]o,kT,#]&/@ Range[0,nmax];
{SymmProb-> All\[CapitalGamma]/Total[All\[CapitalGamma]],LocalisedProbs->  Total[(All\[CapitalGamma]/Total[All\[CapitalGamma]])[[#]] PadRight[P[{0,#}],nmax+1] &/@ Range[1,nmax+1]]}]
End[]
EndPackage[]
