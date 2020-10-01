***********************************************************************
*** dskdm2 returns the full scattering amplitude squared, SUMMED over 
*** both initial and final spin and other internal states, and then 
*** divided by the DM internal degrees of freedom. 
*** The returned value is not evaluated at zero momentum tranasfer, but
*** averaged over t.
***                                                                         
***  type : INTERFACE                                                        
***                                                                         
***  desc : Full scattering amplitude squared, averaged over momentum transfer
***                                                                         
***  input: omega -- CMS MOMENTUM of scattering partner
***  output: omega -- CMS ENERGY of scattering partner
***
***  input:   SMtype - SM scattering partners
***             "    = 7,8,9,10,11,12 - u,d,s,c,b,t quarks
***                    4,5,6 - e,m,t leptons
***                    1,2,3 - e,m,t neutrinos
***
*** author: torsten.bringmann@fys.uio.no, 2016-12-21
*** updated 2018-05-18 : moved momentum->energy conversion to src_models
***********************************************************************

      real*8 function dskdm2(omega,SMtype)
      implicit none
      include 'dssilveira_zee.h'

      real*8  omega
      integer SMtype

      integer kf
      real*8 mf2,mh2, res, kk2, kk2cm, dsmwimp, dsmass


c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('silveira_zee','dskdm2')

      if (SMtype.le.3) mf2=0d0           ! neutrinos
      if (SMtype.eq.4) mf2=dsmass(11)**2 ! PDG code for electrons
      if (SMtype.eq.5) mf2=dsmass(13)**2 ! PDG code for muons
      if (SMtype.eq.6) mf2=dsmass(15)**2 ! PDG code for taus
      if (SMtype.eq.7) mf2=dsmass(1)**2  ! PDG code for up quarks
      if (SMtype.eq.8) mf2=dsmass(2)**2  ! PDG code for down quarks
      if (SMtype.eq.9) mf2=dsmass(3)**2  ! PDG code for strange quarks
      if (SMtype.eq.10) mf2=dsmass(4)**2  ! PDG code for charm quarks
      if (SMtype.eq.11) mf2=dsmass(5)**2  ! PDG code for bottom quarks
      if (SMtype.eq.12) mf2=dsmass(6)**2  ! PDG code for top quarks

      omega = sqrt(mf2+omega**2) ! convert from momentum (input) 
                                 ! to energy (output)

      
      mh2 = mass(khsm)**2

      kk2 = omega**2 - mf2
      kk2cm = kk2 / (1. + 2.*omega/dsmwimp() + mf2/dsmwimp()**2)
      
      if (kk2.lt.0d0) then
        write(*,*) 'dskdm2: ERROR in assignment of particle energies!' ,omega, sqrt(mf2)
      endif  


      if (kk2cm.gt.mh2/1.d4) then  ! full expression from 1706.07433
        res = 4*kk2cm*(2*kk2cm-2*mf2+mh2)/(4*kk2cm+mh2)
     &      -(mh2-2*mf2)*dlog(1.+4*kk2cm/mh2) 
        res = res*lambda**2*mf2/kk2**2/8.d0
      else ! expand log for small kk2 to avoid negative Msq due to numerical precision
        res = mf2/mh2 + (2./3.)*kk2cm/mh2*(1.-8*mf2/mh2) 
     &                - 4*(kk2cm/mh2)**2*(1.-6*mf2/mh2)
        res = res*2*lambda**2*mf2/mh2*kk2cm/kk2
      endif   
     
      if (res.lt.0.d0) then
        write(*,*) 'ERROR in dskdm2: negative |M|^2: ' ,SMtype, res
        res=0.d0
      endif  

c... take into account particle and anti-particle scattering
      if (SMtype.ge.4) res=res*2.d0

      if ((SMtype.ge.7).and.(SMtype.le.12)) res = res*3.d0 ! color factor for quarks

      dskdm2 = res

c... TB test      
c      call dskdm2simp(SMtype,cn,n)
c      dskdm2 = cn*(omega/dsmwimp())**n
                
      return
      end

