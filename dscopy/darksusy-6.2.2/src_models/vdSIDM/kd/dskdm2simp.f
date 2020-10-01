***********************************************************************
*** dskdm2simp returns the scattering amplitude for DM at zero 
*** momentum transfer squared in the EMPTY model, SUMMED over both initial 
*** and final spin and other internal states, and then divided by the DM 
*** internal degrees of freedom. This is only valid in the limit of 
*** relativistic scattering partners with small energies omega, where we 
*** can expand as
***                                                                         
***    |M|**2 = cn*(omega/m0)**n + O( (omega/m0)**(n+1) )
*** 
***  type : interface                                                       
***                                                                         
***  input: Stype   - SM scattering partners:
***                    7,8,9 - u,d,s quarks
***                    4,5,6 - e,m,t leptons
***                    1,2,3 - e,m,t neutrinos
***
***         Stype    > 100 - non-SM scattering partners
***                          (NB: massless MUST be listed first!)
***             "    = 101 - DR particle
***             "    = 102 - mediator particle    
***
*** author: torsten.bringmann@fys.uio.no, 2015-06-22
***********************************************************************

      subroutine dskdm2simp(Stype,cn,n)
      implicit none
      include 'dsvdSIDM.h'
      include 'dsmpconst.h'

      integer n, Stype
      real*8  cn

c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('vdSIDM','dskdm2simp')

      n=0
      cn=0.0d0

      if (Stype.ne.101.and.Stype.ne.102) return    ! only scattering with DR or mediators!

c... DM fermion, vector mediator
      if (DMtype.eq.1.and.Mediatortype.eq.2) then        

        if (Stype.eq.101) then
c... This is Table II (B27) in 1603.04884, valid for mdm >> mmed >> omega
          n  = 2
          cn = 256*gDM**2*gDR**2*mass(kdm)**4/3.0d0/mass(kmed)**4 !sum over DR antiparticles included!
        else
          n  = 0
          cn = 64./3.*gDM**4 ! Table I (B5)
        endif


c... DM fermion, scalar mediator
      elseif (DMtype.eq.1.and.Mediatortype.eq.3) then

        if (Stype.eq.101) then
c... This is Table II (B24) in 1603.04884, valid for mdm >> mmed >> omega
          n  = 2
          cn = 512*gDM**2*gDR**2*mass(kdm)**4/3.0d0/mass(kmed)**4 !sum over DR antiparticles included!
        else
          n  = 0
          cn = 16./3.*gDM**4 ! Table I (B3)
        endif

            
      else
        write(*,*) 'ERROR in dskdm2simp: called with unsupported '
        write(*,*) 'DM/DR/mediator spin combination: ',  
     &              DMtype, Mediatortype, DRtype   
        write(*,*) 'program stopping...'
        stop
      endif


      return
      end





        

