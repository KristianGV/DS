*******************************************************************************
*** Prepares integration of Boltzmann equation and has to be called         ***
*** before dskdboltz (called by dskdtkd)                                    ***
***                                                                         ***
***  type : interface                                                       ***
***                                                                         ***
*** author: torsten bringmann (troms@physto.se), 2010-01-23                 ***
*** updates: 2013-06-11 (made model-dependence explicit)                    ***
*******************************************************************************

      subroutine dskdparticles
      implicit none

      include 'dsmssm.h'
      include 'dskdcom.h'

      real*8  ktmp
      integer i,j,n
      real*8 dsmwimp


c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('MSSM','dskdparticles')

      m0=dsmwimp()

c... are there any BSM scattering partners?
      nBSMscatt=0
      nBSMscattlight=0
      
c... set up potential resonances and order them
      nKDres=0
      n=0
 10     n=n+1
        ktmp=0.d0
        if (n.ge.10) goto 40 ! no more resonances
        if (n.eq.1) ktmp=(mass(ksnu_flav(1,1))**2-m0**2)/(2.*m0)
        if (n.eq.2) ktmp=(mass(ksnu_flav(2,1))**2-m0**2)/(2.*m0)
        if (n.eq.3) ktmp=(mass(ksnu_flav(3,1))**2-m0**2)/(2.*m0)
        if (n.eq.4) then 
           ktmp=(mass(ksl_flav(1,1))**2-m0**2-mass(ke)**2)/(2.*m0)
           if (ktmp.lt.mass(ke)/10.) goto 10
        endif
        if (n.eq.5) then
           ktmp=(mass(ksl_flav(1,2))**2-m0**2-mass(ke)**2)/(2.*m0)
           if (ktmp.lt.mass(ke)/10.) goto 10
        endif
        if (n.eq.6) then
           ktmp=(mass(ksl_flav(2,1))**2-m0**2-mass(kmu)**2)/(2.*m0)
           if (ktmp.lt.mass(kmu)/10.) goto 10
        endif
        if (n.eq.7) then
           ktmp=(mass(ksl_flav(2,2))**2-m0**2-mass(kmu)**2)/(2.*m0)
           if (ktmp.lt.mass(kmu)/10.) goto 10
        endif
        if (n.eq.8) then
           ktmp=(mass(ksl_flav(3,1))**2-m0**2-mass(ktau)**2)/(2.*m0)
           if (ktmp.lt.mass(ktau)/10.) goto 10
        endif
        if (n.eq.9) then
           ktmp=(mass(ksl_flav(3,2))**2-m0**2-mass(ktau)**2)/(2.*m0)
           if (ktmp.lt.mass(ktau)/10.) goto 10
        endif

        if (ktmp.gt.m0/10.) goto 10    ! only keep resonances for relevant 
                                       ! (small) SM particle energies
        nKDres=nKDres+1
        j=1
 20     if (j.lt.nKDres) then
          if (resk(j).le.ktmp) then
            j=j+1
            goto 20
          endif
        endif 
        do 30 i=0,nKDres-j-1
          resk(nKDres-i)=resk(nKDres-i-1)
 30     continue
        resk(j)=ktmp
        goto 10 ! check next resonance

 40   return

      end





        

