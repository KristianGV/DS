************************************************************************
*** function dsIB2yield returns the yield from electroweak internal
*** bremsstrahlung (IB) for SUSY models. 
***
***   Input: egev    - energy in GeV
***          yieldk  - which yield to calculate
***                    (see dshaloyield for an explanation)
***          onshell - 0 subtracts the on-shell contributions ("2-body final 
***                      states")in the narrow width approximation [default]
***                      (note that the result can be negative!)
***                    1 returns the result including on-shell contributions
***   Output: yield [total tree-level(!) annihilation]**-1, 
***                  differential also [GeV]**-1
***           istat - 0 if everything went well
***
*** Author: Torsten Bringmann (torsten.bringmann@fys.uio.no) 
*** Date:   2014-03-03
************************************************************************

      real*8 function dsIB2yield(egev,yieldk,onshell,istat)
      implicit none
      include 'dsib2com.h'

c------------------------ functions ------------------------------------

      real*8 dsIB2yieldone

c------------------------ variables ------------------------------------

      real*8 egev,yield
      integer IB2ch,istat,yieldk,itmp, onshell
c-----------------------------------------------------------------------


      istat=0
      yield=0.d0
      dsIB2yield=0.d0

c... now sum over all IB2 channels:

      do 100 IB2ch=101,112    ! ffZ           
         if (IB2flag(IB2ch).eq.1) then  
            yield=yield+dsIB2yieldone(egev,IB2ch,yieldk,onshell,itmp)
            istat=or(istat,itmp)
         endif
 100  continue

      do 200 IB2ch=201,212      ! ffW; Note that exploiting CP invariance is not 
                                !      straightforward for charged particle yields!
         if (IB2flag(IB2ch).eq.1) then  
            yield=yield+dsIB2yieldone(egev,IB2ch,yieldk,onshell,itmp)
            istat=or(istat,itmp)
         endif
 200  continue

      do 300 IB2ch=301,312    ! ffh           
         if (IB2flag(IB2ch).eq.1) then  
            yield=yield+dsIB2yieldone(egev,IB2ch,yieldk,onshell,itmp)
            istat=or(istat,itmp)
         endif
 300  continue

      do 400 IB2ch=401,412    ! ffH           
         if (IB2flag(IB2ch).eq.1) then  
            yield=yield+dsIB2yieldone(egev,IB2ch,yieldk,onshell,itmp)
            istat=or(istat,itmp)
         endif
 400  continue

      do 500 IB2ch=501,512    ! ffA           
         if (IB2flag(IB2ch).eq.1) then  
            yield=yield+dsIB2yieldone(egev,IB2ch,yieldk,onshell,itmp)
            istat=or(istat,itmp)
         endif
 500  continue

      do 600 IB2ch=601,612     ! ffH+ 
         if (IB2flag(IB2ch).eq.1) then  
            yield=yield+dsIB2yieldone(egev,IB2ch,yieldk,onshell,itmp)
            istat=or(istat,itmp)
         endif
 600  continue


      dsIB2yield=yield
      
      return
      end
