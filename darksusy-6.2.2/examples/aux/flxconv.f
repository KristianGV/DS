      program flxconv
      implicit none
c_______________________________________________________________________
c
c  Program to find how limits from neutrino telescopes derived with a
c  given energy threshold converts to limits of a given other energy
c  threshold or a diffrent kind of limit, like a limit on the cross section.
c  This program relies on the results from the WimpSim code (WimpAnn
c  and WimpEvent). WimpSim takes care of annihilation of WIMPs in the centre
c  of the Sun/Earth and simulates neutrino interactions and oscillations out
c  of the Sun/Earth and close to the detector.
c  Large simulations with this program have been summarized
c  as data tables that are included in DarkSUSY together with tools to 
c  access and interpolate in the tables. Hence, you should link
c  this code to DarkSUSY to use the tables. DarkSUSY also includeds routines
c  to calculate the capture rate in the Sun and Earth and hence conversion
c  can be made not only between different neutrino telescope yields/fluxes
c  but also scattering cross sections. In principle, one can change the 
c  velocity distributions etc DarkSUSY uses, but the defaults should be
c  good enough for most uses.
c  If you want to produce conversion factors for many different cases,
c  the script scr/flxconvbatch.pl is available on an as-is-basis.      
c  Author: Joakim Edsjo, edsjo@fysik.su.se
c  Date: 1999-12-03
c  Modified: 2000-06-14
c  Modified: 2008-04-11 to work with new WIMP annihilation results
c    as implemented in DarkSUSY rev >= 215.
c  Modified: 2009-02-xx 
c  Modified: 2009-04-23 to include Kaluza-Klein channel
c  Modified: 2015-09-25 to work properly with latest DarkSUSY 5.1.2.
c    Uses new capture rate routines (Sun) with more solar elements and
c    a numerical integration over form factors.
c  Modified 2016-06-29 to use DarkSUSY 6.   
c
c=======================================================================

      include 'dsver.h'
      include 'dssecom.h'
      include 'dshmcom.h'

c------------------------------------------------------------ Functions

      real*8 dsseget

c------------------------------------------------------------ Variables

      real*8 mneu,denscorr,cf2a,ca2f,cf2f
      real*8 e1,flux1,th1
      real*8 e2,flux2,th2
      real*8 arateea,aratesu
      real*8 gps,gns,gpa,gna
      real*8 sigsip0,sigsin0,sigsdp0,sigsdn0
      real*8 rho
      integer kind1,type1,kind2,type2
      integer ch,istat,i
      character*2 wh
      character*4 type(6)
      character*80 filename
      data type/'soft','soft','soft','hard','hard','hard'/
c      data denscorr/0.92d0/  ! correction for sims at 1 g/cm^3
      character*17 chname(15)
      data chname/'d d-bar','u u-bar','s s-bar','c c-bar',
     &  'b b-bar','t t-bar','glue glue','W+ W-','Z0 Z0',
     &  'mu+ mu- (zero)','tau+ tau-',
     &  'nu_e nu_e-bar','nu_mu nu_mu-bar','nu_tau nu_tau-bar',
     &  'KK (UED)'/
      character*33 typename(0:31)
      data typename/
     &  'annihilation rate in Earth/Sun',
c...here comes simulated channels
     &  'nu_e at a plane in det.',
     &  'nu_e-bar at a plane in det.',
     &  'nu_mu at a plane in det.',
     &  'nu_mu-bar at a plane in det.',
     &  'nu_tau at a plane in det.',
     &  'nu_tau-bar at a plane in det.',
     &  'e- at neutrino-nucleon vertex',
     &  'e+ at neutrino-nucleon vertex',
     &  'mu- at neutrino-nucleon vertex',
     &  'mu+ at neutrino-nucleon vertex',
     &  'tau- at neutrino-nucleon vertex',
     &  'tau+ at neutrino-nucleon vertex',
     &  'mu- at a plane in the detector',
     &  'mu+ at a plane in the detector',
     &  'hadr. sh. from nu_e CC int.',
     &  'hadr. sh. from nu_e-bar CC int.',
     &  'hadr. sh. from nu_mu CC int.',
     &  'hadr. sh. from nu_mu-bar CC int.',
     &  'hadr. sh. from nu_tau CC int.',
     &  'hadr. sh. from nu_tau-bar CC int.',
     &  'hadr. sh. from nu_e NC int.',
     &  'hadr. sh. from nu_e-bar NC int.',
     &  'hadr. sh. from nu_mu NC int.',
     &  'hadr. sh. from nu_mu-bar NC int.',
     &  'hadr. sh. from nu_tau NC int.',
     &  'hadr. sh. from nu_tau-bar NC int.',
c....here comes summed channels
     &  'mu- + mu+ at neutrino-nucl. vert.',
     &  'mu- + mu+ at a plane in the det.',
     &  'tau- + tau+ at neutr-nucl. vert.',
c....here comes cross sections
     &  'spin-independent cross section',
     &  'spin-dependent cross section'/
      integer typedim(0:31) ! 2=(area)^-1, 3=(volume)^-1, 4=pb
      data typedim/0,6*2,6*3,2*2,12*3,3,2,3,4,4/
      character*11 units(0:4) ! index is typedim above
      data units/
     &  'ann. s^-1',
     &  ' ',
     &  'km^-2 yr^-1',
     &  'km^-3 yr^-1',
     &  'pb'/


c---------------------------------------------------------------- Setup

c...Initialize DarkSUSY
      call dsinit 

c      call dssenu_set('nummed')       ! if you want some other than default

c   ******************************
c   ********** Settings **********
c   ******************************
c...secalcmet determines how the capture rate calculation is perfomed
c...secalcmet=4, full integration over Sun/Earth radius and velocity
c...distribution but analytic form factor integration (exponential)
c...secalcmet=5, full numeric integration over radius, velocity distribution
c...and form factor with many elements included (Earth still assumes
c...exponential form factor).
c...The best most up-to-date calculation is secalcmet=5 (it takes some 
c...time per point though).
      secalcmet=5

c...When we put limits on SI cross sections, we can set limits on
c...the SI scatteriing on protons or neutrons separately or assume that
c...they are the same. This is decided here. If both are set to 1 pb, they
c...are assumed to be equal (default)
      sigsip0=1.d-36  ! cm^2
      sigsin0=1.d-36  ! cm^2

c...When we put limits on SD cross sections, we can set limits on
c...the SD scatteriing on protons or neutrons separately or assume that
c...they are the same. This is decided here. If both are set to 1 pb, they
c...are assumed to be equal. The SD neutron scattering cross section only
c...matters for the Sun and with secalcmet=5 where we can actually scatter
c...on elements with odd neutrons. The default is to only include SD scattering
c...on protons.
      sigsdp0=1.d-36  ! cm^2
      sigsdn0=0.d0    ! cm^2

c...Determine velocity distribution to use for Earth
c...Eventually we want to use the more up to date velocity distributions
c...from Sivertsson & Edsjo, 2012, arXiv:1201.1895, but this is not yet
c...in the public version of DarkSUSY. In the online script, the default has
c...(until Sep 24, 2015) been to use the 'sdbest' distribution from Lundberg
c...and Edsjo, 2004, but this is too pessimistic. Hence we here choose to
c...use a gaussian as the default instead. It is slightly too optimistic, but
c...not by too much (see arXiv:1201.1895, fig. 2 for details).
c      veldfearth='sdbest' ! default in 'old' online script and DS
      veldfearth='gauss' ! better default

c...Default local halo density      
      rho=0.3d0

      write(*,*)
      write(*,*) '                 FLXCONV 2.4'
      write(*,'(A,A)') '       Linked to DarkSUSY version ',dsversion
      write(*,*) ' '
      write(*,*) '       A program to convert from a flux'
      write(*,*) '     to an annihilation rate or another flux'
      write(*,*) '    from WIMP annihilation in the Earth/Sun.'
      write(*,*) ' '
      write(*,*) '                 June 29, 2016'
      write(*,*) '                   J. Edsjo'   
      write(*,*) '              edsjo@fysik.su.se'
      write(*,*) ' '


 10   write(*,*)
      write(*,*) 'Enter WIMP mass [GeV] (0 to quit):'
      read(*,*) mneu
      write(*,*) mneu
      if (mneu.eq.0.0d0) then
        close(15)
        stop
      endif

      write(*,*) 'Annihilation in the Earth(''ea'') or the Sun(''su''):'
      read(*,'(a)') wh
      write(*,*) wh

      write(*,*) 'Annihilation channel?'
      do i=1,15
         write(*,111) i,chname(i)
 111     format(I3,' = ',A)
      enddo
      read(*,*) ch
      write(*,*) ch

      write(*,*) 'Density of target detector material (g/cm^3):'
      write(*,*) '   [0.92 g/cm^3 for ice]'
      read(*,*) denscorr
      write(*,*) denscorr

      write(*,*) 'Now enter STARTING fluxes: '

      write(*,*) 'What type of flux do you have?'
      write(*,112) 0,typename(0)
      do i=1,16
         if (i.lt.16) then
            write(*,112) i,typename(i),i+16,typename(i+16)
         else
            write(*,112) i,typename(i)
         endif
      enddo

      read(*,*) type1
      write(*,*) type1

      if (type1.gt.0.and.type1.lt.30) then
         write(*,*) 'Enter kind of flux:'
         write(*,*) '  1 = integrated above E_min and below theta_max'
         write(*,*) '  2 = differential in energy and angle'
         read(*,*) kind1
         write(*,*) kind1


         if (kind1.eq.1) then
            write(*,*) 'Enter energy threshold:'
         else
            write(*,*) 'Enter energy:'
         endif
         read(*,*) e1
         write(*,*) e1

         write(*,*) 'Enter opening angle (half-aperture) in degrees:'
         read(*,*) th1
         write(*,*) th1
      else
         kind1=0
         e1=0.d0
         th1=0.d0
      endif
       


      if (type1.gt.0) then ! flux of some kind
         if (type1.eq.27) then
           flux1 = dsseget(mneu,e1,th1,ch,wh,kind1,9,istat)
     &       + dsseget(mneu,e1,th1,ch,wh,kind1,10,istat)
         elseif (type1.eq.28) then
           flux1 = dsseget(mneu,e1,th1,ch,wh,kind1,13,istat)
     &       + dsseget(mneu,e1,th1,ch,wh,kind1,14,istat)
         elseif (type1.eq.29) then
           flux1 = dsseget(mneu,e1,th1,ch,wh,kind1,11,istat)
     &       + dsseget(mneu,e1,th1,ch,wh,kind1,12,istat)
         elseif (type1.eq.30) then ! spin-indep scattering
           if (secalcmet.ne.5) then ! old capture routines
              call dssenu_annrate(mneu,rho,sigsip0,0.d0,1.d-10,
     &            secalcmet,arateea,aratesu)
           else ! new routines
c...Below we set sigma_SI^p = sigma_SI^n
              call dsddsig2g(mneu,sigsip0,sigsin0,0.d0,0.d0,
     &           gps,gns,gpa,gna)
       write(*,*) 'Stand by for calculation (can take some time)...'
              call dssenu_annrateff(mneu,rho,gps,gns,gpa,gna,1.d-10,
     &           arateea,aratesu)
           endif
           if (wh.eq.'ea'.or.wh.eq.'EA') then
              flux1 = (arateea*1.d24/3.15d7)
           else
              flux1 = (aratesu*1.d24/3.15d7)
           endif
         elseif (type1.eq.31) then ! spin-dep scattering
           if (secalcmet.ne.5) then ! old routines
              call dssenu_annrate(mneu,rho,0.d0,sigsdp0,1.d-10,
     &            secalcmet,arateea,aratesu)
           else ! new routines
c...Below we use sigma_SD^p, but set sigma_SD^n to zero 
             call dsddsig2g(mneu,0.d0,0.d0,sigsdp0,sigsdn0,
     &          gps,gns,gpa,gna)
             call dssenu_annrateff(mneu,rho,gps,gns,gpa,gna,1.d-10,
     &           arateea,aratesu)
           endif
           if (wh.eq.'ea'.or.wh.eq.'EA') then
              flux1 = (arateea*1.d24/3.15d7)
           else
              flux1 = (aratesu*1.d24/3.15d7)
           endif
         else
           flux1 = dsseget(mneu,e1,th1,ch,wh,kind1,type1,istat)
         endif

c...First calculate conversion ann-rate -> flux
         if (typedim(type1).eq.2) then
            cf2a=3.16d13              ! m^-2 s^-1 -> km^-2 yr^-1
     &        *flux1                   ! yield 
     &        *1.0d-30                ! factor in yield
            cf2a=1.d0/cf2a
         elseif (typedim(type1).eq.3) then
            cf2a=3.16d13              ! m^-2 s^-1 -> km^-2 yr^-1
     &        *flux1                   ! yield 
     &        *1.0d-30                ! factor in yield
            cf2a=cf2a*1000.0d0     ! Another m^-1 -> km^-1 for the conv. rate
            cf2a=cf2a*denscorr     ! tables for dens = 1 g/cm^3
            cf2a=1.d0/cf2a
         elseif (typedim(type1).eq.4) then
            cf2a=flux1
         endif


      else ! annihilation rate directly

         cf2a=1.0d0

      endif


c======================================================================
c======================================================================
c======================================================================

c...Now convert to some other flux
      write(*,*) 
      write(*,*) 'Now enter ENDING fluxes: '

      write(*,*) 'What type of flux do you have?'
      write(*,112) 0,typename(0)
      do i=1,16
         if (i.lt.16) then
            write(*,112) i,typename(i),i+16,typename(i+16)
         else
            write(*,112) i,typename(i)
         endif
 112     format(I3,' = ',A33,' ',I3,' = ',A33)
      enddo

      read(*,*) type2
      write(*,*) type2

      if (type2.gt.0.and.type2.lt.30) then
         write(*,*) 'Enter kind of flux:'
         write(*,*) '  1 = integrated above E_min and below tehta_max'
         write(*,*) '  2 = differential in energy and angle'
         read(*,*) kind2
         write(*,*) kind2


         if (kind2.eq.1) then
            write(*,*) 'Enter energy threshold:'
         else
            write(*,*) 'Enter energy:'
         endif
         read(*,*) e2
         write(*,*) e2

         write(*,*) 'Enter opening angle (half-aperture) in degrees:'
         read(*,*) th2
         write(*,*) th2

      else
         kind2=0
         e2=0.d0
         th2=0.d0
      endif

      if (type2.gt.0) then ! flux of some kind
         if (type2.eq.27) then
           flux2 = dsseget(mneu,e2,th2,ch,wh,kind2,9,istat)
     &       + dsseget(mneu,e2,th2,ch,wh,kind2,10,istat)
         elseif (type2.eq.28) then
           flux2 = dsseget(mneu,e2,th2,ch,wh,kind2,13,istat)
     &       + dsseget(mneu,e2,th2,ch,wh,kind2,14,istat)
         elseif (type2.eq.29) then
           flux2 = dsseget(mneu,e2,th2,ch,wh,kind2,11,istat)
     &       + dsseget(mneu,e2,th2,ch,wh,kind2,12,istat)
         elseif (type2.eq.30) then ! spin-indep scattering
           if (secalcmet.ne.5) then ! old capture routines
              call dssenu_annrate(mneu,rho,sigsip0,0.d0,1.d-10,
     &          secalcmet,arateea,aratesu)
           else ! new routines
c...Below we set sigma_SI^p = sigma_SI^n
              call dsddsig2g(mneu,sigsip0,sigsin0,0.d0,0.d0,
     &           gps,gns,gpa,gna)
       write(*,*) 'Stand by for calculation (can take some time)...'
              call dssenu_annrateff(mneu,rho,gps,gns,gpa,gna,1.d-10,
     &           arateea,aratesu)
           endif
           if (wh.eq.'ea'.or.wh.eq.'EA') then
              flux2 = 1.d0/(arateea*1.d24/3.15d7)
           else
              flux2 = 1.d0/(aratesu*1.d24/3.15d7)
           endif
         elseif (type2.eq.31) then ! spin-dep scattering
           if (secalcmet.ne.5) then ! old routines
              call dssenu_annrate(mneu,rho,0.d0,sigsdp0,1.d-10,
     &          secalcmet,arateea,aratesu)
           else ! new routines
c...Below we use sigma_SD^p, but set sigma_SD^n to zero 
             call dsddsig2g(mneu,0.d0,0.d0,sigsdp0,sigsdn0,
     &         gps,gns,gpa,gna)
             call dssenu_annrateff(mneu,rho,gps,gns,gpa,gna,1.d-10,
     &           arateea,aratesu)
           endif
           if (wh.eq.'ea'.or.wh.eq.'EA') then
              flux2 = 1.d0/(arateea*1.d24/3.15d7)
           else
              flux2 = 1.d0/(aratesu*1.d24/3.15d7)
           endif
         else
           flux2 = dsseget(mneu,e2,th2,ch,wh,kind2,type2,istat)
         endif

c...First calculate conversion ann-rate -> flux
         if (typedim(type2).eq.2) then
            ca2f=3.16d13              ! m^-2 s^-1 -> km^-2 yr^-1
     &        *flux2                   ! yield 
     &        *1.0d-30                ! factor in yield
         elseif (typedim(type2).eq.3) then
            ca2f=3.16d13              ! m^-2 s^-1 -> km^-2 yr^-1
     &        *flux2                   ! yield 
     &        *1.0d-30                ! factor in yield
            ca2f=ca2f*1000.0d0     ! Another m^-1 -> km^-1 for the conv. rate
            ca2f=ca2f*denscorr     ! tables for dens = 1 g/cm^3
         elseif (typedim(type2).eq.4) then
            ca2f=flux2
         endif

      else ! annihilation rate directly

         ca2f=1.0d0

      endif

c...Now put the two conversions together
      cf2f=cf2a*ca2f

 114  format(A17,16x,A17)
 115  format(A12,F8.2,A5,10x,A12,F8.2,A5)
 116  format(A12,I2,21x,A12,I2)
 117  format(A10,A11,14x,A10,A11)
      write(*,*) ' '
      write(*,*) '***** RESULTS *****'
      write(*,*) ' '
      write(*,113) '  mneu  = ',mneu,' GeV    wh = ',wh,'    ch = ',ch,
     &  '    dens = ',denscorr,' g/cm^3'
 113  format(A,F10.2,A,A2,A,I2,A,F8.4,A)
      write(*,*) ' '
      write(*,114) 'STARTING flux','ENDING flux'
      write(*,114) '-------------','-----------'
      write(*,115) '  e1    = ',e1,' GeV ','  e2    = ',e2,' GeV '
      write(*,115) '  th1   = ',th1,' deg.','  th2   = ',th2,' deg.'
      write(*,116) '  kind1 = ',kind1,'  kind2 = ',kind2
      write(*,116) '  type1 = ',type1,'  type2 = ',type2
      write(*,117) '  Unit: ',units(typedim(type1)),
     &   'Unit: ',units(typedim(type2))
      write(*,*) ' '
      write(*,*) '  Conversion from STARTING flux to',
     &   ' annihilation rate (ann. s^-1):'
      write(*,*) '     ====> Conversion factor cf2a = ',cf2a

      write(*,*) ' '
      write(*,*) '  Conversion from annihilation rate (ann. s^-1)',
     &   ' to ENDING flux:'
      write(*,*) '     ====> Conversion factor ca2f = ',ca2f

      write(*,*) ' '
      write(*,*) '  Conversion from STARTING flux to ENDING flux:'
      write(*,*) '     ====> Conversion factor cf2f = ',cf2f

      goto 10 

      end


      real*8 function dsseget(mneu,e1,th1,ch,wh,kind1,type1,istat)
**********************************************************************
*** dsseget is a wrapper routine that call dsseyield_sim to get the
*** flux/yield and properly sums for the KK (UED) channel
*** Author: Joakim Edsjo, edsjo@fysik.su.se
*** Date: 2009-04-23
*** Modified to work with new se_yield in DS 6 by J. Edsjo, 2016-06-29      
**********************************************************************
      implicit none
      real*8 mneu,e1,th1,yield,dsseyield_sim
      character*2 wh
      integer ch,kind1,type1,istat,i

      real*8 mtop
      parameter(mtop=172.4d0) ! should match what WimpSim used

      real*8 kkbr(14)
      data kkbr/.007d0,.113d0,.007d0,.113d0,.007d0,.113d0,
     &  0.0d0,0.0d0,0.0d0,0d0,0.203d0,0.012d0,0.012d0,0.012d0/
c...Note mu+ mu- set to zero above as it gives zero neutrinos
      save kkbr

      integer ch2pdg(14)
      data ch2pdg/1,2,3,4,5,6,21,24,23,13,15,12,14,16/
      save ch2pdg

      real*8 kkbrre(14),kkbrsum0,kkbrsum

      if (ch.le.14) then 
         dsseget=dsseyield_sim(mneu,e1,th1,ch2pdg(ch),'0',wh,
     &     kind1,type1,istat)
      elseif (ch.eq.15) then ! KK (UED)

c...sum up kinematically allowed channels and rescale kkbr properly
         do i=1,14
            if (mneu.ge.mtop) then ! no rescaling
               kkbrre(i)=kkbr(i)
            else
               if (i.eq.6) then
                  kkbrre(i)=0.d0
               else
                  kkbrre(i)=kkbr(i)/(1.d0-kkbr(6))
               endif
            endif
         enddo

         yield=0.d0
         do i=1,14
            if (kkbrre(i).gt.0.d0) then
               yield=yield+kkbrre(i)*
     &              dsseyield_sim(mneu,e1,th1,ch2pdg(i),'0',wh,
     &                kind1,type1,istat)
            endif
         enddo
         dsseget=yield
      else
         write(*,*) 'Incorrect channel in dsseget, ch=',ch
         dsseget=0.d0
      endif

      return
      end


C======================================================================

      subroutine dsddsig2g(mwimp,sigsip,sigsin,sigsdp,sigsdn,
     &  gps,gns,gpa,gna)
c_______________________________________________________________________
c  Calculate WIMP-nucleon four-fermion couplings from cross sections
c
c...NOTE: By calculating the cross section, we have lost the information
c...on the sign of the g's. We will here assume that they are positive and
c...of the same sign, but this does not need to be the case. Hence, only
c...use this routine if you do not have access to the g's directly in any
c...other way.
c
c    'dssusy.h' - file with susy common blocks
c  input:
c    sigsip, sigsin : proton and neutron spin-independent cross sections
c    sigsdp, sigsdn : proton and neutron spin-dependent cross sections
c    units: cm^2
c  output:
c    gps,gns,gpa,gna is as in dsddgpgn
c  author: paolo gondolo (gondolo@lpthe.jussieu.fr) 1994,1995,2002
c     13-sep-94 pg no drees-nojiri twist-2 terms
c     22-apr-95 pg important bug corrected [ft -> ft mp/mq]
c     06-apr-02 pg drees-nojiri treatment added
c     2010-05-28 added this routine to go 'backwards'
c=======================================================================
      implicit none
      include 'dsmpconst.h'

      real*8 sigsip,sigsdp,sigsin,sigsdn,mwimp
      real*8 fkinp,fkinn,gps,gns,gpa,gna

      fkinp = 4./pi*(m_p*mwimp/(m_p+mwimp))**2
      fkinn = 4./pi*(m_n*mwimp/(m_n+mwimp))**2

c...NOTE: By calculating the cross section, we have lost the information
c...on the sign of the g's. We will here assume that they are positive and
c...of the same sign, but this does not need to be the case. Hence, only
c...use this routine if you do not have access to the g's directly in any
c...other way.

      gps=2.0d0*sqrt(sigsip/fkinp/gev2cm2)
      gns=2.0d0*sqrt(sigsin/fkinn/gev2cm2)

      gpa=sqrt(sigsdp/gev2cm2/fkinp*4.d0/3.d0)
      gna=sqrt(sigsdn/gev2cm2/fkinn*4.d0/3.d0)


      return
      end
