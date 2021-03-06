#if 0
* hee-v1_All4_All5_Ind1_Ind2_Neu4.F
* form factors for {S[10, {hno}]} -> {F[2, {gno1}], -F[2, {gno2}]}
* generated by FormCalc 8.4 on 27-Feb-2015 17:20
* this file is part of FeynHiggs
#endif

#include "hee-vars.h"
#include "inline.h"

	subroutine hee_v1_All4_All5_Ind1_Ind2Neu4
	implicit none

#include "hee-vars.h"

#include "inline.h"

        Cloop(1) = 
     &   Cloop(1) + 
     &    1/(16.D0*Pi)*(AlfaMZ*ELMZ*
     &        (cI*(CB*CW*Sub547(All4,All5,Ind1,Ind2)*
     &             XHiggs(hno,3,uzext1)) + 
     &          Sub544(All4,All5,Ind1,Ind2))*
     &        (m3*pave26(PVC(cc2),All4,All5,Neu4)*
     &           Sub554(All4,All5,Neu4) + 
     &          m2*pave26(PVC(cc1),All4,All5,Neu4)*
     &           Sub555(All4,All5,Neu4) + 
     &          MNeu(Neu4)*pave26(PVC(cc0),All4,All5,Neu4)*
     &           Sub556(All4,All5,Neu4)))/
     &      (CB*CB2*CW*CW2*MW*MW2*SW*SW2)

#ifdef DETAILED_DEBUG
	DCOUP "v1_All4_All5_Ind1_Ind2_Neu4: Cloop(1) =", Cloop(1) ENDL
#endif

	end
