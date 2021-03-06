* Uncertainties.F
* Compute (actually, not completely invent) an error estimate
* for the Higgs masses and mixings
* this file is part of FeynHiggs
* last modified 2 Feb 17 th

#include "externals.h"
#include "types.h"
#include "debug.h"


	subroutine FHUncertainties(error,
     &    DeltaMHiggs_, DeltaSAeff_, DeltaUHiggs_, DeltaZHiggs_)
	implicit none

#include "FH.h"
#include "looptools.h"
#define __SUBROUTINE__ "FHUncertainties"

	integer error
	RealType DeltaMHiggs_(NHiggs)
	ComplexType DeltaSAeff_
	ComplexType DeltaUHiggs_(NNeutral,NNeutral)
	ComplexType DeltaZHiggs_(NNeutral,NNeutral)

	RealType HiggsErr(NHiggsErr)
	RealType HiggsSave(NHiggsData)
	RealType HiggsBase(NHiggsErr)
	RealType oldmudim
	integer i, j, err, tMsave

	integer NSfData
	parameter (NSfData = NSf - 9)
* Sf(*,bTR0) minus the CKM copy in bTR0

	RealType SfSave(NSfData)

	RealType SfData(NSfData)
	equivalence (Sf(1,bTR0), SfData)

	RealType DeltaMHiggs(NHiggs)
	ComplexType DeltaSAeff
	ComplexType DeltaXHiggs(0:NNeutral,0:NNeutral,2)
	equivalence (HiggsErr(1), DeltaMHiggs)
	equivalence (HiggsErr(1+NHiggs), DeltaSAeff)
	equivalence (HiggsErr(1+NHiggs+2), DeltaXHiggs)

	do i = 1, NHiggs
	  DeltaMHiggs_(i) = 100000
	enddo
	DeltaSAeff_ = uncomputable
	do i = 1, NNeutral
	  do j = 1, NNeutral
	    DeltaUHiggs_(j,i) = uncomputable
	    DeltaZHiggs_(j,i) = uncomputable
	  enddo
	enddo

	if( looplevel .lt. 2 ) return

	call CheckHiggsCorr(error)

	if( error .gt. 0 ) return

* store the central value

	do i = 1, NHiggsData
	  HiggsSave(i) = HiggsData(i)
	enddo
	do i = 1, NHiggsErr
	  HiggsBase(i) = HiggsData(i)
	enddo
	do i = 1, NHiggsErr
	  HiggsErr(i) = uncomputable
	enddo

	err = 0

* a) variation of mudim from MT/2 to 2 MT

	if( debuglevel .ge. 2 )
     &    DHIGGS "Uncertainties: mudim = 1/2 MT" ENDL

	oldmudim = mudim

	mudim = .25D0*MT2
	tl_valid = 0
	eft_valid = 0
	call CalcHiggsCorr(err)

	if( err .ne. 0 ) then
	  Warning("Could not evaluate corrections at mudim = MT/2")
	  mudim = oldmudim
	  goto 999
	endif

	do i = 1, NHiggsErr
	  if( HiggsData(i) .ne. uncomputable .and.
     &        HiggsBase(i) .ne. uncomputable )
     &      HiggsErr(i) = abs(HiggsData(i) - HiggsBase(i))
	enddo

	if( debuglevel .ge. 2 )
     &    DHIGGS "Uncertainties: mudim = 2 MT" ENDL

	mudim = 4*MT2
	tl_valid = 0
	eft_valid = 0
	call CalcHiggsCorr(err)

	mudim = oldmudim

	if( err .ne. 0 ) then
	  Warning("Could not evaluate corrections at mudim = 2 MT")
	  goto 999
	endif

	do i = 1, NHiggsErr
	  if( HiggsData(i) .ne. uncomputable .and.
     &        HiggsErr(i) .ne. uncomputable )
     &      HiggsErr(i) = max(HiggsErr(i),
     &        abs(HiggsData(i) - HiggsBase(i)))
	enddo

	do i = 1, NSfData
	  SfSave(i) = SfData(i)
	enddo

* b) exclude higher-order MB resummation effects

	if( debuglevel .ge. 2 )
     &    DHIGGS "Uncertainties: no higher-order MB resummation" ENDL

	if( botResum .eq. 0 ) then
	  bM = bTR
	  bM0 = bTR0
	  tl_valid = 0
	  eft_valid = 0
	  call CalcHiggsCorr(err)
	  bM = 4
	  bM0 = 4
	  do i = 1, NHiggsErr
	    HiggsBase(i) = HiggsData(i)
	  enddo
	endif

	botResum = ibset(botResum, 8)
	tl_valid = 0
	eft_valid = 0
	call CalcHiggsCorr(err)
	botResum = ibclr(botResum, 8)

	if( err .ne. 0 ) then
	  Warning("Could not evaluate corrections for unresummed MB")
	  goto 998
	endif

	do i = 1, NHiggsErr
	  if( HiggsData(i) .ne. uncomputable .and.
     &        HiggsErr(i) .ne. uncomputable )
     &      HiggsErr(i) = HiggsErr(i) +
     &        abs(HiggsData(i) - HiggsBase(i))
	enddo

* c) use MTpole instead of MTrun in the two-loop corrections

	if( debuglevel .ge. 2 )
     &    DHIGGS "Uncertainties: MTpole in 2L corrections" ENDL

	tMsave = tM2

	if( tM2 .eq. tT .and. loglevel .eq. 1 ) then
* we re-do the computation with tT2 (i.e. MT(MT) incl EW corr),
* so that the error estimate for loglevel = 1 is not smaller than
* for loglevel = 2
	  tM2 = tT2
	  tM1 = tT2
	  tl_valid = 0
	  eft_valid = 0
	  call CalcHiggsCorr(err)
	  do i = 1, NHiggsErr
	    HiggsBase(i) = HiggsData(i)
	  enddo
	endif

	tM2 = tSelect(tT2,3,3,3, runningMT)
c	if( loglevel .eq. 0 ) tM1 = tM2
	tM1 = tM2

	tl_valid = 0
	eft_valid = 0
	call CalcHiggsCorr(err)

	tM2 = tMsave
	tM1 = tM2

	if( err .ne. 0 ) then
	  Warning("Could not evaluate corrections for MTrun <-> MTpole")
	  goto 998
	endif

	do i = 1, NHiggsErr
	  if( HiggsData(i) .ne. uncomputable .and.
     &        HiggsErr(i) .ne. uncomputable )
     &      HiggsErr(i) = HiggsErr(i) +
     &        abs(HiggsData(i) - HiggsBase(i))
	enddo

998	do i = 1, NSfData
	  SfData(i) = SfSave(i)
	enddo

999	do i = 1, NHiggsData
	  HiggsData(i) = HiggsSave(i)
	enddo
	tl_valid = valid
	eft_valid = valid
	higgs_valid = valid

	if( err .eq. 0 ) then
	  do i = 1, NHiggs
	    DeltaMHiggs_(i) = DeltaMHiggs(i)
	  enddo
	  DeltaSAeff_ = DeltaSAeff
	  do i = 1, NNeutral
	    do j = 1, NNeutral
	      DeltaUHiggs_(j,i) = DeltaXHiggs(j,i,1)
	      DeltaZHiggs_(j,i) = DeltaXHiggs(j,i,2)
	    enddo
	  enddo
	endif
	end
