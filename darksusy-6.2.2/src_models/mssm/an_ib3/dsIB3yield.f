*****************************************************************************
***   function dsIB3yield gives the yield resulting from from gluon internal 
***   bremsstrahlung (IB) for SUSY models, i.e. the *difference* to the 2-body
***   yield when taking into account qqg final states. 
***   Input: egev - energy in GeV
***          yieldk - which yield to calculate
***             (see dshaloyield for an explanation, 
***              currently only photon and antiproton yields are implemented)
***   Output: yield (annihilation)**1, differential also GeV**-1
***           istat is set as follows in case of errors
***   bit  decimal  reason
***     0        1  dsIBf_intdy (kinematic integration for 3-body rate) failed
***
*** Author: torsten.bringmann@fys.uio.no
*** Date: 2015-06-01
*****************************************************************************

      real*8 function dsIB3yield(egev,yieldk,istat)
      implicit none
      include 'dsibcom.h'

c------------------------ functions ------------------------------------

      real*8 dsIB3yieldone, dssigmav0tot, dssigmav0

c------------------------ variables ------------------------------------

      real*8 egev,yield,sv
      real*8 IBabr(12)
      integer IBch,istat,yieldk,itmp

c----------------------------------------------- set-up common variables
 
      istat=0
      yield=0.0d0

c... add IB3 contribution, using same channel numbering as for photon IB

      sv=dssigmav0tot()
      if (yieldk.eq.52.or.yieldk.eq.152.or.
     &    yieldk.eq.54.or.yieldk.eq.154) then 
        IBabr(7)=dssigmav0(2,-2)/sv    ! u u-bar
        IBabr(8)=dssigmav0(1,-1)/sv    ! d d-bar
        IBabr(9)=dssigmav0(4,-4)/sv    ! c c-bar
        IBabr(10)=dssigmav0(3,-3)/sv   ! s s-bar
        IBabr(11)=dssigmav0(6,-6)/sv   ! t t-bar
        IBabr(12)=dssigmav0(5,-5)/sv   ! b b-bar

        if (ibhow.eq.2) call dsibselect ! dynamically choose important channels
        do 200 IBch=7,12  
          if ((IBflag(IBch).ne.0).and.(IBabr(IBch).gt.0.0d0)) then  
            yield=yield+IBabr(IBch) * dsIB3yieldone(egev,IBch,yieldk,itmp)
            istat=or(istat,itmp)
          endif
  200   continue
      endif

      dsIB3yield=yield

      end
