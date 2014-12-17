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
	Molecular dynamics simulation of a one-dimensional anharmonic chain of particles
	with a pairwise shoulder-shaped interaction potential.
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

PredictCollision[r_,\[CapitalDelta]p_,{VdCore_,VdPlat_,Vh_}]:={\[Infinity],      \[CapitalDelta]p,            \[Infinity],            NoColl}    /;(r>VdPlat \[And]\[CapitalDelta]p>=0)
PredictCollision[r_,\[CapitalDelta]p_,{VdCore_,VdPlat_,Vh_}]:={VdPlat,-Sqrt[\[CapitalDelta]p^2-4Vh],(VdPlat-r)/\[CapitalDelta]p,PlatIn}    /;(r>VdPlat \[And]\[CapitalDelta]p<0\[And]\[CapitalDelta]p^2-4Vh>=0)
PredictCollision[r_,\[CapitalDelta]p_,{VdCore_,VdPlat_,Vh_}]:={VdPlat, Sqrt[\[CapitalDelta]p^2+4Vh],(VdPlat-r)/\[CapitalDelta]p,PlatOut}   /;(r<=VdPlat\[And]\[CapitalDelta]p>=0)
PredictCollision[r_,\[CapitalDelta]p_,{VdCore_,VdPlat_,Vh_}]:={VdPlat,-\[CapitalDelta]p,            (VdPlat-r)/\[CapitalDelta]p,PlatBounce}/;(r>VdPlat \[And]\[CapitalDelta]p<0\[And]\[CapitalDelta]p^2-4Vh<0)
PredictCollision[r_,\[CapitalDelta]p_,{VdCore_,VdPlat_,Vh_}]:={VdCore,-\[CapitalDelta]p,            (VdCore-r)/\[CapitalDelta]p,Core}      /;(r<=VdPlat\[And]\[CapitalDelta]p<0)

NextCollision[r_List,p_List,Vpot_]:=Module[{\[CapitalDelta]p=FieldDifferences[p,0]},
	First[Sort[Table[Append[PredictCollision[r[[i]],\[CapitalDelta]p[[i]],Vpot],i],{i,Length[r]}],#1[[3]]<#2[[3]]&]]]

CollisionInfo[{q_List,p_List,t_},L_,Vpot_]:=Module[{r=FieldDifferences[q,L],\[CapitalDelta]p,\[CapitalDelta]t,type,i,inext,pmean,p1,Jp,Je},
	{\[CapitalDelta]p,\[CapitalDelta]t,type,i}=NextCollision[r,p,Vpot][[2;;-1]];
	inext=Mod[i+1,Length[q],1];
	pmean=Mean[p[[{i,inext}]]];
	p1=pmean+{-\[CapitalDelta]p/2,\[CapitalDelta]p/2};
	(*momentum and energy current*)
	Jp=(\[CapitalDelta]p-p[[inext]]+p[[i]])/2;
	Je=pmean Jp+If[type===PlatOut,Last[Vpot],If[type===PlatIn,-Last[Vpot],0]]/2;
	{t+\[CapitalDelta]t,type,i,{Jp,Je}}]

UpdateField[{q_List,p_List,t_},L_,Vpot_]:=Module[{r=FieldDifferences[q,L],\[CapitalDelta]p,\[CapitalDelta]t,i,q1,p1},
	{\[CapitalDelta]p,\[CapitalDelta]t,i}=NextCollision[r,p,Vpot][[{2,3,5}]];
	q1=q+\[CapitalDelta]t p;
	i={i,Mod[i+1,Length[q],1]};
	p1=p;
	p1[[i]]=Mean[p[[i]]]+{-\[CapitalDelta]p/2,\[CapitalDelta]p/2};
	{q1,p1,t+\[CapitalDelta]t}]

AdvanceField[{q_List,p_List,t_},\[CapitalDelta]t_]:={q+\[CapitalDelta]t p,p,t+\[CapitalDelta]t}

(* use 'AdvanceField' to ensure that processed collision is not registered again (Subscript[r, i] at potential boundary!) *)
FieldSequence[qpt0_,L_,Vpot_,nsteps_]:=NestList[UpdateField[AdvanceField[#,0.0001],L,Vpot]&,qpt0,nsteps]


(* Visualization *)

VisualizeField[q_List,{VdCore_,VdPlat_,Vh_},width_]:=Graphics[{
	Point[{#,0}]&/@q,
	Red,Circle[{#,0},VdCore/2]&/@q,
	Orange,Circle[{#,0},VdPlat/2]&/@q,
	Gray,Line[{{0,0},{width,0}}]},ImageSize->Large,PlotRange->{{-1,width+2},All}]
