*****************************************************************************
***   function dsIByield_fsr gives the photon yield from final state radiation 
***   (FSR) from the deacy of a hypothetical partical with mass 2*m0, such as
***   included in the Pythia runs. Just like in Pythia, only FSR from fermionic
***   final states is included
***   Input: egev - energy in GeV
***          yieldk - 52 for integrated or 152 for differential yield
***   Output: yield (annihilation)**1, differential also GeV**-1
***           istat is set as follows in case of errors
***   bit  decimal  reason
***     0        1  dsIBf_intdxdy failed
***     1        2  dsIBf_intdy failed
*** Author: Torsten Bringmann, troms@physto.se
*** Date: 2008-02-10
*****************************************************************************

      real*8 function dsIByield_fsr(egev,yieldk,istat)
      implicit none
 

c------------------------ functions ------------------------------------

      real*8 dsIByieldone_fsr,dsmwimp

c------------------------ variables ------------------------------------

      real*8 egev,yield,sv,dssigmav0tot,dssigmav0
      real*8 FSRabr(4:12)               ! only fermions
      integer FSRch,istat,yieldk,itmp


      istat=0
      yield=0.0d0

      sv=dssigmav0tot()
      if (yieldk.eq.52.or.yieldk.eq.152) then 
        FSRabr(4)  = dssigmav0(11,-11)/sv    ! e- e+
        FSRabr(5)  = dssigmav0(13,-13)/sv    ! mu- mu+
        FSRabr(6)  = dssigmav0(15,-15)/sv    ! tau- tau+
        FSRabr(7)  = dssigmav0(2,-2)/sv    ! u u-bar
        FSRabr(8)  = dssigmav0(1,-1)/sv    ! d d-bar
        FSRabr(9)  = dssigmav0(4,-4)/sv    ! c c-bar
        FSRabr(10) = dssigmav0(3,-3)/sv   ! s s-bar
        FSRabr(11) = dssigmav0(6,-6)/sv   ! t t-bar
        FSRabr(12) = dssigmav0(5,-5)/sv   ! b b-bar

c... now compute the contribution from each channel
        do 200 FSRch=4,12     
          if (FSRabr(FSRch).gt.0.0d0) then  
            yield=yield+FSRabr(FSRch)*
     &             dsIByieldone_fsr(dsmwimp(),egev,FSRch,yieldk,itmp)
           endif
  200   continue
      endif

      dsIByield_fsr=yield

      end
