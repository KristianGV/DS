***********************************************************************
*** dskdm2simp returns the scattering amplitude for DM at zero 
*** momentum transfer squared in the geenric WIMP model, SUMMED over both initial 
*** and final spin and other internal states, and then divided by the DM 
*** internal degrees of freedom. This is only valid in the limit of 
*** relativistic scattering partners with small energies omega, where we 
*** can expand as
***                                                                         
***    |M|**2 = cn*(omega/m0)**n + O( (omega/m0)**(n+1) )
*** 
***  type : interface                                                       
***                                                                         
***  input: SMtype   - SM scattering partners:
***                    7,8,9 - u,d,s quarks
***                    4,5,6 - e,m,t leptons
***                    1,2,3 - e,m,t neutrinos
***
*** author: torsten.bringmann@fys.uio.no, 2015-06-22
*** updated 2019-12-13 : added photon scattering
***********************************************************************

      subroutine dskdm2simp(SMtype,cn,n)
      implicit none
      include 'dsgeneric_wimp.h'
      include 'dsmpconst.h'

      integer n, SMtype
      real*8  cn

      real*8 s, dsmwimp,dsmass, mf2

c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('generic_wimp','dskdm2simp')

c... assume constant amplitude (contact interaction) and that dominant scattering
c... channel is the same as dominant annihilation channel
      n=0
      cn=0.0d0

c... FIXME: handle this conversion table centrally...      
      if ((SMtype.eq.1.and.svch.eq.12).or.    ! nu_e
     &    (SMtype.eq.2.and.svch.eq.14).or.    ! nu_mu
     &    (SMtype.eq.3.and.svch.eq.16).or.    ! nu_tau
     &    (SMtype.eq.4.and.svch.eq.11).or.    ! e
     &    (SMtype.eq.5.and.svch.eq.13).or.    ! mu
     &    (SMtype.eq.6.and.svch.eq.15).or.    ! tau
     &    (SMtype.eq.7.and.svch.eq.1).or.     ! d
     &    (SMtype.eq.8.and.svch.eq.2).or.     ! u
     &    (SMtype.eq.9.and.svch.eq.3).or.     ! s
     &    (SMtype.eq.13.and.svch.eq.22)) then ! photon

        if (SMtype.le.3) mf2=0d0           ! neutrinos
        if (SMtype.eq.4) mf2=dsmass(11)**2 ! PDG code for electrons
        if (SMtype.eq.5) mf2=dsmass(13)**2 ! PDG code for muons
        if (SMtype.eq.6) mf2=dsmass(15)**2 ! PDG code for taus
        if (SMtype.eq.7) mf2=dsmass(2)**2  ! PDG code for up quarks
        if (SMtype.eq.8) mf2=dsmass(1)**2  ! PDG code for down quarks
        if (SMtype.eq.9) mf2=dsmass(3)**2  ! PDG code for strange quarks
        if (SMtype.eq.13) mf2=0.0d0

        s=4*dsmwimp()**2
        cn=4*pi*s/sqrt(1.0d0-4*mf2/s)*sva/gev2cm3s
        cn=cn*kdof(kdm)  ! DM d.o.f.
        
c... take into account particle and anti-particle scattering
        if (SMtype.ge.4) cn=cn*2.d0

c (color factor already taken into account in conversion from annihilation cross section)   
     
      endif


      return
      end





        

