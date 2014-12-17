(* ::Package:: *)

(* Mathematica Package *)

(* :Author: Christian B. Mendl
	http://christian.mendl.net *)

(* :Copyright: Copyright (C) 2013, Christian B. Mendl
	All rights reserved.
	This program is free software; you can redistribute it and/or
	modify it under the terms of the Simplified BSD License
	http://www.opensource.org/licenses/bsd-license.php
*)

(* :Summary:
	Simple molecular dynamics simulation of a hard-point chain with alternating masses.
*)

(* :References:

	[1] Christian B. Mendl, Herbert Spohn
	Current fluctuations for anharmonic chains in thermal equilibrium
	arXiv:1412.4609

	[2] Christian B. Mendl, Herbert Spohn
	Equilibrium time-correlation functions for one-dimensional hard-point systems
	Phys. Rev. E 90, 012147 (2014), arXiv:1403.0213

	[3] Christian B. Mendl, Herbert Spohn
	Dynamic correlators of Fermi-Pasta-Ulam chains and nonlinear fluctuating hydrodynamics
	Phys. Rev. Lett. 111, 230601 (2013), arXiv:1305.1209

	[4] Herbert Spohn
	Nonlinear fluctuating hydrodynamics for anharmonic chains
	J. Stat. Phys. 154, 1191-1227 (2014), arXiv:1305.6412
*)


(* Simulation *)

FieldDifferences[q_List,L_]:=Differences[Append[q,q[[1]]+L]]

PredictCollision[r_,\[CapitalDelta]v_,{m0_,m1_},VdCore_]:={                 0, \[Infinity]           }/;\[CapitalDelta]v>=0
PredictCollision[r_,\[CapitalDelta]v_,{m0_,m1_},VdCore_]:={2 m0 m1 \[CapitalDelta]v/(m0+m1),(VdCore-r)/\[CapitalDelta]v}/;\[CapitalDelta]v<0

NextCollision[r_List,p_List,mass_,VdCore_]:=Module[{\[CapitalDelta]v=FieldDifferences[p/mass,0]},
	First[Sort[Table[Append[PredictCollision[r[[i]],\[CapitalDelta]v[[i]],mass[[{i,Mod[i+1,Length[r],1]}]],VdCore],i],{i,Length[r]}],#1[[2]]<#2[[2]]&]]]

UpdateField[{q_List,p_List,t_},L_,mass_,VdCore_]:=Module[{r=FieldDifferences[q,L],\[CapitalDelta]p,\[CapitalDelta]t,i,q1,p1},
	{\[CapitalDelta]p,\[CapitalDelta]t,i}=NextCollision[r,p,mass,VdCore];
	q1=q+\[CapitalDelta]t p/mass;
	i={i,Mod[i+1,Length[q],1]};
	p1=p;
	p1[[i]]+={\[CapitalDelta]p,-\[CapitalDelta]p};
	{q1,p1,t+\[CapitalDelta]t}]

AdvanceField[{q_List,p_List,t_},mass_,\[CapitalDelta]t_]:={q+\[CapitalDelta]t p/mass,p,t+\[CapitalDelta]t}

(* use 'AdvanceField' to ensure that processed collision is not registered again (r_i at potential boundary!) *)
FieldSequence[qpt0_,L_,mass_,VdCore_,nsteps_]:=NestList[UpdateField[AdvanceField[#,mass,0.0001],L,mass,VdCore]&,qpt0,nsteps]


(* Visualization *)

VisualizeField[q_List,relmass_,VdCore_,width_]:=Graphics[{
	Table[{If[relmass[[i]]<1,Green,Blue],Disk[{q[[i]],0},VdCore/2]},{i,Length[q]}],
	Red,Circle[{#,0},VdCore/2]&/@q,
	Gray,Line[{{0,0},{width,0}}]},ImageSize->Large,PlotRange->{{-1,width+2},All}]
