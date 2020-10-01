***********************************************************************
*** dskdm2 returns the full scattering amplitude for DM at zero 
*** momentum transfer squared in the generic WIMP model, SUMMED over both initial 
*** and final spin and other internal states, and then divided by the DM 
*** internal degrees of freedom.
***                                                                         
***  type : interface                                                       
***                                                                         
***  desc : Full scattering amplitude squared, at zero momentum transfer
***                                                                         
***  input: omega -- CMS MOMENTUM of scattering partner
***  output: omega -- CMS ENERGY of scattering partner
***
***  input:   SMtype - SM scattering partners
***             "    = 7,8,9 - u,d,s quarks
***                    4,5,6 - e,m,t leptons
***                    1,2,3 - e,m,t neutrinos
***                  = 13    - photons
***
*** author: torsten.bringmann@fys.uio.no, 2015-06-22
*** updated 2018-05-14 : moved momentum->energy conversion to src_models/
*** updated 2019-12-13 : added photon scattering
***********************************************************************

      real*8 function dskdm2(omega,SMtype)
      implicit none

      real*8  omega,mf
      integer SMtype

      real*8 dsmwimp,dsmass,cn
      integer n

c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('generic_wimp','dskdm2')

      dskdm2 = 0.0d0
      if ((SMtype.lt.1).or.(SMtype.gt.13)) return

      if (SMtype.le.3) mf=0d0         ! neutrinos
      if (SMtype.eq.4) mf=dsmass(11)  ! PDG code for electrons
      if (SMtype.eq.5) mf=dsmass(13)  ! PDG code for muons
      if (SMtype.eq.6) mf=dsmass(15)  ! PDG code for taus
      if (SMtype.eq.7) mf=dsmass(1)   ! PDG code for up quarks
      if (SMtype.eq.8) mf=dsmass(2)   ! PDG code for down quarks
      if (SMtype.eq.9) mf=dsmass(3)   ! PDG code for strange quarks
      if (SMtype.eq.10) mf=dsmass(4)  ! PDG code for charm quarks
      if (SMtype.eq.11) mf=dsmass(5)  ! PDG code for bottom quarks
      if (SMtype.eq.12) mf=dsmass(6)  ! PDG code for top quarks
      if (SMtype.eq.13) mf=0.0d0      ! photon mass

      
      omega = sqrt(mf**2+omega**2) ! convert from momentum (input) 
                                   ! to energy (output)


c... For the simplified model, only the expanded amplitude is available
c... (here exceptionally also for light quarks)
      call dskdm2simp(SMtype,cn,n)
      
      dskdm2=cn*(omega/dsmwimp())**n
                
      return
      end

