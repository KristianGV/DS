***********************************************************************
*** dskdm2 returns the full scattering amplitude for DM at zero 
*** momentum transfer squared in the EMPTY model, SUMMED over both initial 
*** and final spin and other internal states, and then divided by the DM 
*** internal degrees of freedom. 
***                                                                         
***  type : interface                                                       
***                                                                         
***  desc : Full scattering amplitude squared at zero momentum transfer
***                                                                         
***  input: omega -- CMS MOMENTUM of scattering partner
***  output: omega -- CMS ENERGY of scattering partner
***
***  input:   SMtype - SM scattering partners
***             "    = 7,8,9 - u,d,s quarks
***                    4,5,6 - e,m,t leptons
***                    1,2,3 - e,m,t neutrinos
***
*** author: torsten.bringmann@fys.uio.no, 2014-05-09
***********************************************************************

      real*8 function dskdm2(omega,SMtype)
      implicit none

      real*8  omega
      integer SMtype

c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('empty','dskdm2')
      
      dskdm2=0.d0
    
      return

      end

