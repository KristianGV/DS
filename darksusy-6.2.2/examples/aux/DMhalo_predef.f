************************************************************************
*** The program DMhalo_predef demonstrates initialization and basic usage of 
*** the DM density profiles already pre-defined in parametric forms in 
*** the DS release. Quick reading guide:
***
*** line 55 (117, 137) ff : initialization of NFW (Burkert, Einasto) profile
*** line 153 ff : compute Jfactors and gamma-ray fluxes
*** line 217 ff : demonstrate how to scan over halo parameters
*** line 249 ff : reaidng out profile values for given halo
*** line 279 ff : example on how to re-define input parameters to
***               characterize a halo (here: local vs. scale density of NFW)
*** line 322 ff : auxiliary routines needed by main programme 
***
************************************************************************
      program DMhalo_predef
      implicit none
ccc
      real*8 mwimp,sv,SI  ! generic wimp parameters
      integer pdgann
      logical selfconj      
ccc
      character*12 tag1,tag2,tag3,tag4,tag
      real*8 rhos,rs,distance,rinn,rmax,alpha
ccc
      real*8 psi0,alphas,jfactor,dsjfactor,egev,xi,flux,dsgafluxsph
      integer diff,istat
ccc
      real*8 dsrnd1
      integer idum,ii      
      real*8 radius,res1,res2,res3,dsdmdsph
      real*8 rho0,r0
      logical check
ccc      
ccc Initialize DarkSUSY
      call dsinit
ccc      
ccc sample generic WIMP model setup:
ccc      
      mwimp=100.d0              ! WIMP mass in GeV
      selfconj=.true.           ! self-conjugated WIMP
      sv=3.d-27                 ! WIMP pair annihilation cm^3/s
      pdgann=5                  ! b b-bar
      SI=1.d-5                  ! WIMP-proton scattering (pb)
                                ! irrelevant for this program
ccc      
      call dsgivemodel_generic_wimp(mwimp,selfconj,sv,pdgann,SI)
ccc
      write(*,*) 'Sample results in case of a WIMP of mass = ',mwimp
      write(*,*) 'annihilating into the b-bbar channel'

ccc
ccc example on how to define a dark matter profile using one of the
ccc three parametric forms available in the DS release:
ccc      
ccc 1) a spherical NFW profile is defined by:
ccc
ccc   rho(r) = rhos/((r/rs)*(1+r/rs)^2)
ccc
ccc so you must specify the two parameters in this expression:     
ccc
      rhos=0.73d0 ! GeV cm^-3
      rs=1.8d0    ! kpc, NOTE: keep these units to have all rates
                  ! computed consistently!
ccc
ccc plus three extra parameters you need to set for any spherical halo:
ccc
      distance=76.d0 ! kpc, how far is the center of the DM profile
                     ! from the observer
      rinn=1.d-3  ! kpc, some inner cutoff radius, this is not  
                  ! very important except for very cuspy profiles
      rmax=5.d0 ! kpc, the radial size of the DM halo, the exact     
                ! value is not very important except for slowly
                ! decling density profiles
ccc
ccc you then need to specify a tag for this model: a tag can be up to 12
ccc characters, and substrings of characters within it are used to provide
ccc information to the code. 
ccc   'nfw' -> you tell the code that this is a NFW profile
ccc   'mw'-> you tell the code this profile refers to the Milky Way.
ccc          NOTE: this is an important info to specify, because  
ccc          otherwise the code will assume it cannot be used to compute
ccc          DM rates that are Milky Way specific!
ccc   'hdb' -> you tell the code to store this model in the halo
ccc          database. NOTE: only such hdb models are saved for the 
ccc          whole run and can be used to read/write tabulated
ccc          quantities entering the DM rates. Other models (except for
ccc          those with 'def' see below) are treated as temporary
ccc          models; only one temporary model can be active at any time,
ccc          when switching to another, the previous model - as well as
ccc          stored quantities relative to that model - are overwritten!
ccc   'def' -> these are the DS default models, you can add more of them
ccc          and they will be treated as 'hdb' models.
ccc the full list of strings can be found in the include files       
ccc 'dsdmdcom.h' & 'dsdmddrvrcom.h'
ccc
      tag1='nfwhdb1'
ccc
ccc you then need to create this profile: this is done via a call
ccc to the subroutine dsdmdset_halomodel; this involves a proper match
ccc between parameter values and parameter identification tags, which is
ccc provided in the auxiliary routine sethalo_nfw included below:
ccc
      call sethalo_nfw(tag1,rhos,rs,distance,rinn,rmax)
ccc
ccc a tag containing 'hdb' or 'def' must be uniquely defined, since DM 
ccc rates for that profile can be computed simply using that tag as an
ccc input argument of a function/subroutine, so the attempt to store 
ccc another model with the same tag returns an error (it returns an 
ccc error even when attempting to store the same model twice). So after 
ccc changing the database model, we need to change the label. E,g,
ccc      
      rhos=rhos*2.d0
      tag2='nfwhdb2'
      call sethalo_nfw(tag2,rhos,rs,distance,rinn,rmax)
ccc      
ccc Another possibility:
ccc 2) a spherical Burkert profile is defined by:
ccc
ccc   rho(r) = rhos/(1+(r/rs))/(1+(r/rs)^2)
ccc
ccc so you must specify again rhos and rs and the three extra
ccc parameters:
ccc
      rhos=7.18d0  ! GeV cm^-3
      rs=0.57d0    ! kpc
      distance=76.d0 ! kpc
      rinn=1.d-3  ! kpc
      rmax=5.d0 ! kpc
ccc
ccc and initialize it with the analogous routine sethalo_bur (see below) after 
ccc defining a tag which contains the string 'bur'
ccc      
      tag3='burhdb1'
      call sethalo_bur(tag3,rhos,rs,distance,rinn,rmax)
ccc      
ccc The third possibility is:
ccc 3) a spherical Einasto profile defined by:
ccc
ccc   rho(r) = rhos*exp(-2/alpha*((r/rs)^alpha-1))
ccc
ccc just proceed as above, just make sure 'ein' appears in the tag:
ccc
      rhos=0.156d0  ! GeV cm^-3
      rs=1.8d0    ! kpc
      alpha=0.18d0
      distance=76.d0 ! kpc
      rinn=1.d-3  ! kpc
      rmax=5.d0 ! kpc
      tag4='einhdb1'
      call sethalo_ein(tag4,rhos,rs,alpha,distance,rinn,rmax)

ccc      
ccc Use the models that have been just created for DM rates:
ccc   e.g. compute Jfactors and gamma-ray fluxes:      
ccc      
      psi0=0.d0     ! degrees, direction of observation
      alphas=0.5d0  ! degrees, aperture of the instrument acceptance
      write(*,*)
      write(*,*) 'Compute Jfactors in the direction psi = ',psi0
     &  ,' [deg]'
      write(*,*) 'within a cone of aperture theta = ',alphas
     &  ,' [deg]'
      psi0=psi0*4.d0*datan(1.d0)/180.d0 ! in rad
      alphas=alphas*4.d0*datan(1.d0)/180.d0 ! in rad
ccc      
      jfactor=dsjfactor(tag1,psi0,alphas) ! GeV^2 cm^-6 kpc sr
      jfactor=jfactor*3.08d21*2.247d-7 ! Msun^2 kpc^-5 sr
      write(*,*) 'For model: ',tag1,
     &  '  Log_10(J(0.5 deg)/[Msun^2 kpc^-5 sr]) = ',
     &  dlog10(jfactor)
ccc      
      jfactor=dsjfactor(tag2,psi0,alphas) ! GeV^2 cm^-6 kpc sr
      jfactor=jfactor*3.08d21*2.247d-7 ! Msun^2 kpc^-5 sr
      write(*,*) 'For model: ',tag2,
     &  '  Log_10(J(0.5 deg)/[Msun^2 kpc^-5 sr]) = ',
     &  dlog10(jfactor)
ccc      
      jfactor=dsjfactor(tag3,psi0,alphas) ! GeV^2 cm^-6 kpc sr
      jfactor=jfactor*3.08d21*2.247d-7 ! Msun^2 kpc^-5 sr
      write(*,*) 'For model: ',tag3,
     &  '  Log_10(J(0.5 deg)/[Msun^2 kpc^-5 sr]) = ',
     &  dlog10(jfactor)
ccc      
      jfactor=dsjfactor(tag4,psi0,alphas) ! GeV^2 cm^-6 kpc sr
      jfactor=jfactor*3.08d21*2.247d-7 ! Msun^2 kpc^-5 sr
      write(*,*) 'For model: ',tag4,
     &  '  Log_10(J(0.5 deg)/[Msun^2 kpc^-5 sr]) = ',
     &  dlog10(jfactor)
ccc      
ccc the computation of the gamma-ray flux needs J-factors calculations;
ccc if psi0 and alphas are not changed, values stored in dsjfactor are 
ccc used rather than recomputing them:
ccc      
      egev=1.d0 ! energy at which the flux is computed [GeV]
      diff=1    ! differential flux
      xi=1.d0   ! rescaling factor
      write(*,*)
      write(*,*) 'Compute the corresponding differential gamma-ray'
      write(*,*) 'flux at the energy = ',egev,' [GeV]'
      write(*,*) 'in the direction psi = ',psi0,' [deg]'
      write(*,*) 'within a cone of aperture theta = ',alphas
     &  ,' [deg]'
      flux=dsgafluxsph(egev,diff,xi,tag1,psi0,alphas,istat)
      write(*,*) 'For model: ',tag1,' flux = ',flux,
     & ' [cm^-2 s^-1 GeV^-1]'
      flux=dsgafluxsph(egev,diff,xi,tag2,psi0,alphas,istat)
      write(*,*) 'For model: ',tag2,' flux = ',flux,
     & ' [cm^-2 s^-1 GeV^-1]'
      flux=dsgafluxsph(egev,diff,xi,tag3,psi0,alphas,istat)
      write(*,*) 'For model: ',tag3,' flux = ',flux,
     & ' [cm^-2 s^-1 GeV^-1]'
      flux=dsgafluxsph(egev,diff,xi,tag4,psi0,alphas,istat)
      write(*,*) 'For model: ',tag4,' flux = ',flux,
     & ' [cm^-2 s^-1 GeV^-1]'

ccc
ccc When peforming, e.g., a scan one needs to generate a large number of
ccc models and there is no need to store all of them in memory. The
ccc rules to create TEMPORARY models are the same as for the models 
ccc above except that 'hdb' or 'def' MUST NOT appear in the model tag.
ccc The tag is needed to tell the code the rules to generate the model,
ccc but it does not need anymore to be uniquely defined, nor it is 
ccc needed to retrive a model: there is only one temporary model active
ccc at any give time and calling a DM rate with any halo tag not to be
ccc searched in the halo model database returns a result valid for the 
ccc temporary profile in its latest initialization.      
ccc      
ccc sample scan:      
ccc
      tag='bur'
      distance=76.d0 ! kpc
      rinn=1.d-3  ! kpc
      rmax=5.d0 ! kpc
      ii=0
      idum=-4173044
 10   ii=ii+1
      rhos=6.d0+2.d0*dsrnd1(idum)  ! GeV cm^-3
      rs=0.4d0+0.8d0*dsrnd1(idum)    ! kpc
      call sethalo_bur(tag,rhos,rs,distance,rinn,rmax)
ccc
      flux=dsgafluxsph(egev,diff,xi,'last',psi0,alphas,istat)
      write(*,*) 'For temporary Burkert model with rhos = ',rhos,
     & ' [GeV cm^-3]' 
      write(*,*) 'and rs = ',rs,' [kpc]','  flux = ',flux,
     & ' [cm^-2 s^-1 GeV^-1]'
      if(ii.lt.10) goto 10 

ccc
ccc While functions like dsgafluxsph and dsjfactor have a tag among
ccc their arguments to allow to connect explicitly to models in the halo
ccc database, you can also call in the main file other quantities
ccc depending on the halo model, in which the dependence is implicit
ccc and the computation is carried out for the currently active profile
ccc (as for dsgafluxsph above in case of a temporary profile); you then
ccc have to make sure to load the proper profile: this is done via the
ccc subroutine dsdmdselect_halomodel. E.g. to compute the DM profile at
ccc a given radius for 3 of the models defined above, you need to
ccc proceed as follows:
ccc
      write(*,*)
      write(*,1000) 'radius      '
     &             ,'NFW profile '
     &             ,'Burkert prof'
     &             ,'Einasto prof'
      do ii=0,20
        radius=dexp(dlog(rinn)+(dlog(rmax*0.99d0)-dlog(rinn))/20.d0*ii) 
        call dsdmdselect_halomodel(tag1)
        res1=dsdmdsph(radius)
        call dsdmdselect_halomodel(tag3)
        res2=dsdmdsph(radius)
        call dsdmdselect_halomodel(tag4)
        res3=dsdmdsph(radius)
        write(*,1001) radius,res1,res2,res3
      enddo
 1000 format(8(3x,a12))
 1001 format(10(1x,E14.6))      

ccc
ccc Finally notice that DM density profiles referring to the Milky Way
ccc can be all defined replacing the scale density with the local halo
ccc density. An example of how to do that is provided in the subroutine
ccc sethalo_nfwmw; below we define a profile which is just two times
ccc larger than the default NFW profile in DS 'mwnfwdef' (loaded here 
ccc through the call to dsinit) 
ccc      
      rho0=0.3d0 ! local halo density [GeV cm^-3]
      rho0=2.d0*rho0
      rs=20.d0   ! scale radius [kpc]
      r0=8.d0    ! local Galactocentric distance [kpc]
      rinn=1.d-5 ! inner cutoff [kpc]
      rmax=2.d2  ! outer truncation radius [kpc]
      tag='mwnfw1' ! temporary NFW to be applied to the Milky Way
      call sethalo_nfwmw(tag,rho0,rs,r0,rinn,rmax)
ccc
      check=.false.
      do ii=0,20
        radius=dexp(dlog(rinn)+(dlog(rmax*0.99d0)-dlog(rinn))/20.d0*ii) 
        call sethalo_nfwmw(tag,rho0,rs,r0,rinn,rmax)
        res1=dsdmdsph(radius)
        call dsdmdselect_halomodel('mwnfwdef')
        res2=dsdmdsph(radius)
        res2=res2*2.d0
        if(dabs(res1-res2).gt.1.d-8*dabs(res1+res2)) check=.true.
      enddo
      write(*,*)
      if(check) then
        write(*,*) ' Check on the generation of Milky Way profile ',
     &       'not passed'
      else
        write(*,*) ' Check on the generation of Milky Way profile ',
     &       'passed'
      endif  
      stop
      end ! main program 


c_______________________________________________________________________
c_______________________ AUXILIARY ROUTINES ____________________________
c_______________________________________________________________________


      subroutine sethalo_nfw(tag,rhos,rs,distance,rinn,rmax)
c_______________________________________________________________________
c
c sample routine with proper parameter-label assignment for the
c intialization of a spherical NFW profile:
c
c   rho(r) = rhos/((r/rs)*(1+r/rs)^2)
c
c inputs:
c   - tag: label for the model, it must contain the string 'nfw'
c   - rhos, rs: the parameters in the formula above
c   - distance: distance from the center of the object
c   - rinn: inner truncation radius, i.e. rho(r)=rho(rinn) for r<rinn
c   - rmax: maximum radius, i.e. rho(r)=0 for r>rmax   
c distances to be given in kpc, rhos in GeV cm^-3
c_______________________________________________________________________
      implicit none
      character*12 tag
      real*8 rhos,rs,distance,rinn,rmax
      logical match
      integer nin
      real*8 rein(5)
      character*10 chin(5)
ccc
ccc check that tag contains 'nfw', otherwise the initialization
ccc will fail      
ccc
      call dslabcheck(12,tag,3,'nfw',match)
      if(.not.match) then
        write(*,*) 'DS: call to sethalo_nfw with tag = ',tag
        write(*,*) 'DS: not containg the string - nfw'
        write(*,*) 'DS: program stopped'
        stop
      endif  
ccc
ccc the order in which each of the 5 parameters is assigned is arbitrary
ccc you just need to be careful to assign each predefined parameter
ccc label chin(i) to the parameter value rein(i), with i going from 1 
ccc to 5:      
ccc
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
ccc initialization:
      call dsdmdset_halomodel(tag,nin,chin,rein)
      return
      end



      
      subroutine sethalo_bur(tag,rhos,rs,distance,rinn,rmax)
c_______________________________________________________________________
c
c sample routine with proper parameter-label assignment for the
c intialization of a spherical Burkert profile:
c
c   rho(r) = rhos/(1+(r/rs))/(1+(r/rs)^2)
c
c inputs:
c   - tag: label for the model, it must contain the string 'bur'
c   - rhos, rs: the parameters in the formula above
c   - distance: distance from the center of the object
c   - rinn: inner truncation radius, i.e. rho(r)=rho(rinn) for r<rinn    
c   - rmax: maximum radius, i.e. rho(r)=0 for r>rmax   
c distances to be given in kpc, rhos in GeV cm^-3
c_______________________________________________________________________
      implicit none
      character*12 tag
      real*8 rhos,rs,distance,rinn,rmax
      logical match
      integer nin
      real*8 rein(5)
      character*10 chin(5)
ccc
ccc check that tag contains 'bur', otherwise the initialization
ccc will fail      
ccc
      call dslabcheck(12,tag,3,'bur',match)
      if(.not.match) then
        write(*,*) 'DS: call to sethalo_bur with tag = ',tag
        write(*,*) 'DS: not containg the string - bur'
        write(*,*) 'DS: program stopped'
        stop
      endif  
ccc
ccc the order in which each of the 5 parameters is assigned is arbitrary
ccc you just need to be careful to assign each predefined parameter
ccc label chin(i) to the parameter value rein(i), with i going from 1 
ccc to 5:      
ccc
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
ccc initialization:
      call dsdmdset_halomodel(tag,nin,chin,rein)
      return
      end



      
      subroutine sethalo_ein(tag,rhos,rs,alpha,distance,rinn,rmax)
c_______________________________________________________________________
c
c sample routine with proper parameter-label assignment for the
c intialization of a spherical Burkert profile:
c
c   rho(r) = rhos*exp(-2/alpha*((r/rs)^alpha-1))
c
c inputs:
c   - tag: label for the model, it must contain the string 'bur'
c   - rhos, rs, alpha: the parameters in the formula above
c   - distance: distance from the center of the object
c   - rinn: inner truncation radius, i.e. rho(r)=rho(rinn) for r<rinn    
c   - rmax: maximum radius, i.e. rho(r)=0 for r>rmax   
c distances to be given in kpc, rhos in GeV cm^-3
c_______________________________________________________________________
      implicit none
      character*12 tag
      real*8 rhos,rs,alpha,distance,rinn,rmax
      logical match
      integer nin
      real*8 rein(6)
      character*10 chin(6)
ccc
ccc check that tag contains 'ein', otherwise the initialization will 
ccc fail:     
ccc
      call dslabcheck(12,tag,3,'ein',match)
      if(.not.match) then
        write(*,*) 'DS: call to sethalo_ein with tag = ',tag
        write(*,*) 'DS: not containg the string - ein'
        write(*,*) 'DS: program stopped'
        stop
      endif  
ccc
ccc the order in which each of the 6 parameters is assigned is arbitrary
ccc you just need to be careful to assign each predefined parameter
ccc label chin(i) to the parameter value rein(i), with i going from 1 
ccc to 6:      
ccc
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
ccc initialization:
      call dsdmdset_halomodel(tag,nin,chin,rein)
      return
      end


      subroutine sethalo_nfwmw(tag,rho0,rs,r0,rinn,rmax)
c_______________________________________________________________________
c
c sample routine with proper parameter-label assignment for the
c intialization of a spherical NFW profile as applicable to the Milky
c Way:
c
c   rho(r) = rho0*((r0/rs)*(1+r0/rs)^2)/((r/rs)*(1+r/rs)^2)
c
c inputs:
c   - tag: label for the model, it must contain the string 'nfw'
c   - rho0, rs: the parameters in the formula above
c   - r0: local galactocentric distance
c   - rinn: inner truncation radius, i.e. rho(r)=rho(rinn) for r<rinn
c   - rmax: maximum radius, i.e. rho(r)=0 for r>rmax   
c distances to be given in kpc, rhos in GeV cm^-3
c_______________________________________________________________________
      implicit none
      character*12 tag
      real*8 rho0,rs,r0,rinn,rmax
      logical match
      integer nin
      real*8 rein(5)
      character*10 chin(5)
ccc
ccc check that tag contains 'nfw', otherwise the initialization
ccc will fail      
ccc
      call dslabcheck(12,tag,3,'nfw',match)
      if(.not.match) then
        write(*,*) 'DS: call to sethalo_nfw with tag = ',tag
        write(*,*) 'DS: not containg the string - nfw'
        write(*,*) 'DS: program stopped'
        stop
      endif  
ccc
ccc check that tag contains 'mw', otherwise the initialization
ccc will fail      
ccc
      call dslabcheck(12,tag,2,'mw',match)
      if(.not.match) then
        write(*,*) 'DS: call to sethalo_abgmw with tag = ',tag
        write(*,*) 'DS: not containg the string - mw'
        write(*,*) 'DS: program stopped'
        stop
      endif  
ccc
ccc the order in which each of the 5 parameters is assigned is arbitrary
ccc you just need to be careful to assign each predefined parameter
ccc label chin(i) to the parameter value rein(i), with i going from 1 
ccc to 5:      
ccc
      chin(1)='objdist'
      rein(1)=r0
      chin(2)='radintr'
      rein(2)=rinn
      chin(3)='radouttr'
      rein(3)=rmax
      chin(4)='rho0'
      rein(4)=rho0
      chin(5)='rsnfw'
      rein(5)=rs
      nin=5
ccc initialization:
      call dsdmdset_halomodel(tag,nin,chin,rein)
      return
      end
      
