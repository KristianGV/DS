*****************************************************************************
*** auxiliary function that selects integrand for integration routines
*** author: Torsten Bringmann 2008-03-12
*** update 2014-11-12: changed call to dsanyield_sim to PDG codes
*** update 2016-04-12: added light quark channels
*****************************************************************************

      real*8 function dsIBintsel2(xint) 

      implicit none
      include 'dsibcom.h'

      real*8 xint,xintcut
      real*8 dsIBf_intdy,dsIBf_intdy2,dsanyield_sim
      integer istat,pdg,yieldkk,yieldpdg,diff

c... intch = 1..12 for IB photons
c... intch = 101..112 for FSR photons
c... intch = -1...-12 for IB positrons

      if ((intch.ge.1.and.intch.le.12).or.
     &    (intch.ge.101.and.intch.le.112)) then
         dsIBintsel2 = 
     &      dsIBf_intdy(intch,xint,ibcom_mx,ibcom_mp1,ibcom_mp2)
      elseif (intch.eq.-4) then
         dsIBintsel2 =
     &      dsIBf_intdy2(4,xint,ibcom_mx,ibcom_mp1,ibcom_mp2)
      elseif ((intch.eq.-1).or.
     &        (intch.le.-5.and.intch.ge.-12)) then
c... set channels for dsanyield_sim 
         if (intch.eq.-1) then
              pdg=24 ! W+
           elseif (intch.eq.-4) then
              pdg=11 !e-
           elseif (intch.eq.-5) then
              pdg=13 !mu-
           elseif (intch.eq.-6) then
              pdg=15 !tau+
           elseif (intch.eq.-7) then
              pdg=2 ! up 
           elseif (intch.eq.-8) then
              pdg=1 ! down 
           elseif (intch.eq.-9) then
              pdg=4 ! charm 
           elseif (intch.eq.-10) then
              pdg=3 ! strange 
           elseif (intch.eq.-11) then
              pdg=6  ! top
           elseif (intch.eq.-12) then
              pdg=5  ! bottom
           else
              dsIBintsel2=0.d0
              return
         endif
         xintcut=xint                  ! for higher energies, take the IB 
         if (xint.ge.0.8) xintcut=0.8  ! contribution at xp=IRcut as a
                                       ! *conservative* estimate for the yield;
                                       ! a full treatment would have to include 
                                       ! virtual photons

c... map to PDG code input required by dsanyield_sim !

      yieldkk=mod(intyield,100)
      diff=intyield/100
      if (diff.ne.0.and.diff.ne.1) then
        write(*,*) 'ERROR in dsIBintsel2: unspoorted argument intyield =', intyield
        dsIBintsel2=0.0d0
        return
      endif
      
      if (yieldkk.eq.51) then
         yieldpdg = -11 ! positron yields
      elseif (yieldkk.eq.52) then 
         yieldpdg = 22  ! cont. gammas     
      elseif (yieldkk.eq.53) then 
         yieldpdg = 14  ! muon neutrinos     
      elseif (yieldkk.eq.54) then 
         yieldpdg = -2212 ! antiproton yields     
      elseif (yieldkk.eq.71) then 
         yieldpdg = 14  ! neutrino yields (same as 53)    
      elseif (yieldkk.eq.72) then 
         yieldpdg = 130072 ! muon yields at creation     
      elseif (yieldkk.eq.73) then 
         yieldpdg = 130073 ! integrated muon yields in ice     
      else
        write(*,*) 'ERROR in dsIBintsel2: unspoorted intyield =', intyield
        dsIBintsel2=0.0d0
        return
      endif

         dsIBintsel2 =
     &        dsIBf_intdy2(-intch,xintcut,ibcom_mx,ibcom_mp1,ibcom_mp2)*
     &          dsanyield_sim(xint*ibcom_mx,ibcom_x*ibcom_mx,
     &                    pdg,0,yieldpdg,diff,istat)
      else
         dsIBintsel2=0.d0
      endif
 
      end



