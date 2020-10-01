***********************************************************************
*** dskdm2 returns the full scattering amplitude squared, SUMMED over 
*** both initial and final spin and other internal states, and then 
*** divided by the DM internal degrees of freedom. 
*** The returned value is not evaluated at zero momentum tranasfer, but
*** averaged over t.
***                                                                         
***  type : interface                                                       
***                                                                         
***  desc : Full scattering amplitude squared, averaged over momentum transfer
***                                                                         
***  input: omega -- CMS MOMENTUM of scattering partner
***  output: omega -- CMS ENERGY of scattering partner
***
***  input:   Stype - SM scattering partners
***             "    = 10,11,12 - c,b,t quarks
***                    7,8,9 - u,d,s quarks
***                    4,5,6 - e,m,t leptons
***                    1,2,3 - e,m,t neutrinos
***
***           Stype  > 100 - non-SM scattering partners
***                          (NB: massless MUST be listed first!)
***             "    = 101 - DR particle
***             "    = 102 - mediator particle    
***
*** author: torsten.bringmann@fys.uio.no, 2018-05-14
***********************************************************************

      real*8 function dskdm2(omega,Stype)
      implicit none
      include 'dsvdSIDM.h'

      real*8  omega
      integer Stype

      real*8 cn
      integer n

c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('vdSIDM','dskdm2')

      dskdm2=0.0d0
      if (Stype.eq.101) then
        omega=sqrt(omega**2+mass(kdr))**2 ! convert from DR momentum (input) 
                                          ! to DR energy (output)  
      elseif (Stype.eq.102) then
        omega=sqrt(omega**2+mass(kmed))**2 ! mediator scattering
      else
        return ! no other scattering in this model!
      endif    
      

c... so far, only leading order scattering contribution implemented:
      call dskdm2simp(Stype,cn,n)


      dskdm2 = cn*(omega/mass(kdm))**n

      
      return
      end

