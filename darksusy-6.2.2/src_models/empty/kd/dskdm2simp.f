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
***  desc : Scattering amplitude squared for zero momentum transfer
***
***  input: SMtype   - SM scattering partners:
***                    7,8,9 - u,d,s quarks
***                    4,5,6 - e,m,t leptons
***                    1,2,3 - e,m,t neutrinos
***
*** author: torsten.bringmann@fys.uio.no, 2014-05-09
***********************************************************************

      subroutine dskdm2simp(SMtype,cn,n)
      implicit none

      integer n, SMtype
      real*8  cn

c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('empty','dskdm2simp')

      n=2
      cn=0.d0

      return

      end





        

