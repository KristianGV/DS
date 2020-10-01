**********************************************************************
*** wimpyields. This program calculates various stable particle yields
*** (positrons, neutrinos, gamma-rays,...) resulting from final states
*** of DM annihilation or decay. It is thus an example of how to use 
*** DarkSUSY without relying on a specific SUSY model (instead linking 
*** to the generic WIMP model).
*** Author: Joakim Edsjo, edsjo@fysik.su.se
*** Date: October 25, 2017
**********************************************************************

      program wimpyields
      implicit none

      real*8 dsanyield_sim, dsseyield_sim
      integer annpdg,yieldpdg,diff,ndec,ntot
      real*8 yield,mwimp,decmin,e
      character*1 hel
      integer i,istat
      character*80 filename


      call dsinit

c...Setup, see header of src/an_yield/dsanyield_sim.f for details
      annpdg=5                  ! PDG code of annihilation products
      yieldpdg=22               ! PDG code of yield of interest
      mwimp=100.d0              ! WIMP mass (GeV)
      hel='0'                   ! helicity state, '0'=unpolarized
      diff=1                    ! differential yields
      decmin=-5.d0              ! lowest energy, 10^-decmin * mwimp
      ndec=20                   ! number of bins per decade
                                ! total number of bins = ndec*(-decmin)
      filename='yield.dat'      

c...Preliminary calculations
      ntot=ndec*(-decmin)

c...Open file and perform calculation      
      open(unit=47,file=filename,status='unknown',form='formatted')
      if (diff.eq.1) then
         write(47,'(A)') '#  i        E [GeV]   yield [GeV^-1 ann^-1]'
      else
         write(47,'(A)') '#  i        E [GeV]   yield [ann^-1]'
      endif
      
      do i=1,ntot
         e=mwimp*10**(decmin-(dble(i)-0.5d0)/abs(ntot)*decmin)
         yield=dsanyield_sim(mwimp,e,annpdg,hel,yieldpdg,diff,istat)
         write(47,100) i,e,yield
 100     format(I4,1x,E14.6,1x,E14.6)
      enddo

      end
