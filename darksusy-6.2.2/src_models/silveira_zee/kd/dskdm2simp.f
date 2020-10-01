***********************************************************************
*** dskdm2simp returns the scattering amplitude squared, averaged over t.
*** Here, |M|^2 is SUMMED over both initial and final spin and 
*** other internal states, and then divided by the DM internal degrees of 
*** freedom. This is only valid in the limit of relativistic scattering 
*** partners with small energies omega, where we can expand as
*** and averaged over the 
***                                                                         
***    <|M|**2>_t = cn*(omega/m0)**n + O( (omega/m0)**(n+1) )
*** 
***  type : INTERFACE                                                        
***                                                                         
***  input: SMtype   - SM scattering partners:
***                    4,5,6 - e,m,t leptons
***                    1,2,3 - e,m,t neutrinos
***
*** author: torsten.bringmann@fys.uio.no, 2016-06-29
***********************************************************************

      subroutine dskdm2simp(SMtype,cn,n)
      implicit none
      include 'dssilveira_zee.h'

      integer n, SMtype
      real*8  cn

      real*8 mf,mh, dsmwimp, dsmass


c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('silveira_zee','dskdm2simp')

      n=0
      cn=0
      if (SMtype.ge.7) return
      
      if (SMtype.le.3) mf=0d0        ! neutrinos
      if (SMtype.eq.4) mf=dsmass(11) ! PDG code for electrons
      if (SMtype.eq.5) mf=dsmass(13) ! PDG code for muons
      if (SMtype.eq.6) mf=dsmass(15) ! PDG code for taus
      
      mh = mass(khsm)

      n = 2
      cn = 4./3.d0 * lambda**2 * mf**2 * dsmwimp()**2 / mh**4

c... take into account particle and anti-particle scattering
      if (SMtype.ge.4) cn=cn*2.d0
      
c      if ((SMtype.ge.7).and.(SMtype.le.9)) cn = cn*3.

      return
      end





        

