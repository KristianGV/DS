      program DMhalo_los
c
c     This program is an example DarkSUSY main program to demonstrate 
c     how to perform line-of-sight integrations for non-trivial
c     regions of interest (ROI) in the Milky Way halo.
c     For simplicity, this is only done for temporary halo models,
c     i.e. we don't store the halo parameters in the halo database.
c     (as, e.g., in DMhalo_predef)
c     A variety of pre-implemented ROI shapes is specified in the function 
c     mask(l,b) at the end of the file.
c
c     Author: Torsten Bringmann, 2018-11-25
c     
c-----This line is 72 columns long--------------------------------------
c
      implicit none      
      include 'dsidtag.h'   ! to have access to which particle module is used
      include 'dsmpconst.h' ! to have access to pi
      include 'dsio.h'      ! to control prtlevel

      real*8 tstart,tfinish                    ! aux
      integer opt
      character charopt
      character*12 halolabel                   ! halo label
      real*8 cospsi0, psi0, theta0             ! direction of l.o.s. integration
      real*8 jpsi                              ! J-factor in direction psi
      real*8 rhos,rs,distance,rinn,rmax,alpha,rho0  ! halo parameters
      integer nminhx,nmaxhx, niterhx,ierrhx    ! healpix
      real*8 epshx,epsabshx      

c... functions
      real*8 jfunc, mask, dsjfactor
      external jfunc, mask

      call CPU_TIME(tstart)
      call dsinit
      prtlevel=2 ! to make sure we get warning messages 
                 ! from the HEALPIX integration 
                 ! (2 also shows how the integration converges)

c... The following are standard halo parameters that we just set by hand 
c... (for galactic center J-factors no need to change anything here)
      distance=8.5d0 ! kpc, distance of halo to earth
      rinn=1.d-3     ! kpc, inner cutoff (only relevant for very cuspy profiles)
      rmax=50.d0     ! kpc,  radial size of DM halo
      rho0=0.4d0     ! GeV/cm^3, local DM density

c... Settings for the HEALPIX integration. For singular profiles and ROIs 
c... much smaller than ~1º very high accuracies may be needed. 
      nminhx=6        ! minimal HEALPIX resolution (level) used
      nmaxhx=13       ! maximal HEALPIX resolution (level) used 
      epshx=1.d-3     ! relative accuracy goal
      epsabshx=1.d-5  ! absolute accuracy goal

      
100   write(*,*)
      write(*,*) '-------------------------------------------------------'
      write(*,*) 'Choose a Milky Way halo parameterization'
      write(*,*) '-------------------------------------------------------'
      write(*,*) 'a) NFW     -> rho(r) = rhos/((r/rs)*(1+r/rs)^2)' 
      write(*,*) 'b) Einasto -> rho(r) = rhos*exp(-2/alpha*((r/rs)^alpha-1))' 
      write(*,*) 'c) Burkert -> rho(r) = rhos/(1+(r/rs))/(1+(r/rs)^2)' 
      write(*,*) '[rhos in all cases fixed by local rho0=0.4 GeV/cm^3]'
      read(*,*) charopt
      if (charopt.ne.'a'.and.charopt.ne.'b'.and.charopt.ne.'c') then
        write(*,*) 'ERROR: unsupported option -- quitting program.'
        goto 500
      endif
      write(*,*) 'Scale radius rs [kpc]:'
      read (*,*) rs
      if (charopt.eq.'a') then
         rhos = rho0*(distance/rs)*(1+distance/rs)**2
c         write(*,*) 'Scale density rhos [GeV cm^-3]:'
c         read (*,*) rhos
     
c... NB: Even if we don't store the halo in the database, we still must specify 
c... a label. It must contains the identifying strings 'nfw', 'bur' or 'ein'.
         halolabel='dummynfw'
         call sethalo_nfw(halolabel,rhos,rs,distance,rinn,rmax)
      endif  
      if (charopt.eq.'b') then
         alpha=0.17
c         write(*,*) 'alpha:'
c         read (*,*) alpha        
         rhos = rho0*exp(2./alpha*((distance/rs)**alpha-1.))
c         write(*,*) 'Scale density rhos [GeV cm**-5]:'
c         read (*,*) rhos
         halolabel='dummyein'
         call sethalo_ein(halolabel,rhos,rs,alpha,distance,rinn,rmax)
      endif  
      if (charopt.eq.'c') then
         rhos = rho0*(1.+distance/rs)*(1+(distance/rs)**2)
c         write(*,*) 'Scale density rhos [GeV cm**-5]:'
c         read (*,*) rhos
         halolabel='dummybur'
         call sethalo_bur(halolabel,rhos,rs,distance,rinn,rmax)
      endif  

      write(*,*)
      write(*,*) '-------------------------------------------------------'
      write(*,*) 'Line-of-sight-integration'
      write(*,*) '-------------------------------------------------------'
      write(*,*) '  1 = calculate J-factor for a circle around the GC (fast)'
      write(*,*) '  2 = calculate J-factor for an arbitrary region (slow)'
      write(*,*) '      [as specified in function mask(l,b)]'
      read (*,*) opt
      write(*,*)
      
      if (opt.eq.1) then
         psi0=0.0d0 ! Direction to GC [in degree]
         write(*,*) '(Half) opening angle of observation cone [deg]:'
         read(*,*) theta0
         theta0=theta0 * pi/180.
         jpsi=dsjfactor(halolabel,psi0,theta0) 

      elseif (opt.eq.2) then
        call jfunclbset(halolabel) ! initializes jfunc(l,b)=j(l,b)*mask(l,b)
                                   ! see below for def. of ROI via mask(l,b)
        call dshealpixint(jfunc,nminhx,nmaxhx,epshx,epsabshx,jpsi,
     &                    niterhx,ierrhx)
c        write (*,*) 'healpix level actually used = ',niterhx
      else
        goto 500
      endif

      write(*,*)
      write(*,*) 'J_factor = ', jpsi, ' kpc GeV^2 cm^-6'
      write(*,*) '         = ', jpsi*kpc*1.0d21, ' GeV^2 cm^-5'
      write(*,*)

500   call CPU_TIME(tfinish)
      write (*,*)
      write (*,*) '-------------------------------------------------------'
      write (*,*) 'The DarkSUSY example program has finished successfully.'
      write (*,*) 'Total time needed (in seconds): ',tfinish-tstart
      write (*,*) 'Particle module that was used: ', moduletag
      write (*,*) '-------------------------------------------------------'
      write (*,*)
      stop
 999  end ! end of main program



********************************************************************
***   Subroutines used by the main program.                      ***
********************************************************************

      real*8 function mask(l,b)
c... Function mask returns 1.0d0 inside, and 0.0d0 outside an ROI
c... defined in terms of galactic longitude l [rad] and latitude b [rad]
        implicit none
        include 'dsmpconst.h'
        real*8 l,b
        character*12 ROItype
        real*8 b1,b2,l1,l2, psimin, psimax, alpha
      
        mask = 0.0d0

c CHOOSE ROI TYPE HERE
        ROItype = 'ring'        
        
        if (ROItype.eq.'ring') then
c------------------------------------
c SPECIFY 'ring' ROI PARAMETERS HERE
c------------------------------------
           psimin = 0.0d0               ! ROI defined by inner and outer angle
           psimax = 2.5d0 * pi/180.d0
c------------------------------------
           if ((cos(b)*cos(l)).ge.cos(psimax).and. ! ('Pythagoras on a sphere')
     &        (cos(b)*cos(l)).le.cos(psimin)) 
     &        mask = 1.0d0
        
        elseif (ROItype.eq.'ringnoplane') then 
c------------------------------------
c SPECIFY 'ringnoplane' ROI PARAMETERS HERE
c------------------------------------
           b1 = 0.3d0 * pi/180.      ! ROI defined by b > |b1| inside psimax
           psimax = 1.0d0 * pi/180.d0
c------------------------------------
           if (cos(b).ge.cos(b1).and.(cos(b)*cos(l)).ge.cos(psimax)) 
     &        mask = 1.0d0

        elseif (ROItype.eq.'hourglass') then  
c------------------------------------
c SPECIFY 'hourglas' ROI PARAMETERS HERE
c------------------------------------
           psimax = 2.5d0 * pi/180.     ! ROI defined by upper and lower
c------------------------------------   ! 45º quadrant, out to angle psimax
           if (abs(b).ge.abs(l).and.
     &         (cos(b)*cos(l)).ge.cos(psimax)) mask = 1.0d0
        
        endif        
        return
      end ! mask

      
      real*8 function jfunc(l,b)
c... returns j*mask; l,b given in radians
      implicit none
      include 'dsdmdcom.h'  ! psannihi,dmdradouttr
      include 'dslosidrvrcom.h' ! ilosigasph
      include 'dsmpconst.h'
      real*8 l,b
      real*8 cospsi0,theta0in,outradius,dslosisph, mask
      integer how,power
      character*12 labhalojbl
      logical loadlab
      common/labhalojblcom/labhalojbl,loadlab
      jfunc=0.0d0
      if (mask(l,b).eq.0.0d0) return
c      jfunc=1.0d0
c      return
      cospsi0=dcos(l)*dcos(b)
      theta0in=dacos(cospsi0)
      outradius=dmdradouttr
      how=ilosigasph  ! prompt 
      power=2  ! annihilations
      jfunc=dslosisph(theta0in,outradius,how,power,labhalojbl,
     &       loadlab)
      return
      end ! jfunc

      subroutine jfunclbset(labhalo)    
c... needed to initialze jfunc(l,b)    
      implicit none
      character*12 labhalo, labhalojbl
      logical loadlab
      common/labhalojblcom/labhalojbl,loadlab
      labhalojbl=labhalo
      call dsdmdselect_halomodel(labhalo)
      loadlab=.false.  ! do not load again the halo in dsomlosisph
      return
      end
      




********************************************************************
***   profile setting routines from DMhalo_predef.               ***
***   See there for detailed comments!                           ***   
********************************************************************
      subroutine sethalo_nfw(tag,rhos,rs,distance,rinn,rmax)
      implicit none
      character*12 tag
      real*8 rhos,rs,distance,rinn,rmax
      logical match
      integer nin
      real*8 rein(5)
      character*10 chin(5)
      call dslabcheck(12,tag,3,'nfw',match)
      if(.not.match) then
        write(*,*) 'DS: call to sethalo_nfw with tag = ',tag
        write(*,*) 'DS: not containg the string - nfw'
        write(*,*) 'DS: program stopped'
        stop
      endif  
      chin(1)='objdist'
      rein(1)=distance
      chin(2)='radintr'
      rein(2)=rinn
      chin(3)='radouttr'
      rein(3)=rmax
      chin(4)='rhosnfw'
      rein(4)=rhos
      chin(5)='rsnfw'
      rein(5)=rs
      nin=5
      call dsdmdset_halomodel(tag,nin,chin,rein)
      return
      end
      
      subroutine sethalo_bur(tag,rhos,rs,distance,rinn,rmax)
      implicit none
      character*12 tag
      real*8 rhos,rs,distance,rinn,rmax
      logical match
      integer nin
      real*8 rein(5)
      character*10 chin(5)
      call dslabcheck(12,tag,3,'bur',match)
      if(.not.match) then
        write(*,*) 'DS: call to sethalo_bur with tag = ',tag
        write(*,*) 'DS: not containg the string - bur'
        write(*,*) 'DS: program stopped'
        stop
      endif  
      chin(1)='objdist'
      rein(1)=distance
      chin(2)='radintr'
      rein(2)=rinn
      chin(3)='radouttr'
      rein(3)=rmax
      chin(4)='rhosbur'
      rein(4)=rhos
      chin(5)='rsbur'
      rein(5)=rs
      nin=5
      call dsdmdset_halomodel(tag,nin,chin,rein)
      return
      end

      subroutine sethalo_ein(tag,rhos,rs,alpha,distance,rinn,rmax)
      implicit none
      character*12 tag
      real*8 rhos,rs,alpha,distance,rinn,rmax
      logical match
      integer nin
      real*8 rein(6)
      character*10 chin(6)
      call dslabcheck(12,tag,3,'ein',match)
      if(.not.match) then
        write(*,*) 'DS: call to sethalo_ein with tag = ',tag
        write(*,*) 'DS: not containg the string - ein'
        write(*,*) 'DS: program stopped'
        stop
      endif  
      chin(1)='objdist'
      rein(1)=distance
      chin(2)='radintr'
      rein(2)=rinn
      chin(3)='radouttr'
      rein(3)=rmax
      chin(4)='rhosein'
      rein(4)=rhos
      chin(5)='rsein'
      rein(5)=rs
      chin(6)='alphaein'
      rein(6)=alpha
      nin=6
      call dsdmdset_halomodel(tag,nin,chin,rein)
      return
      end
