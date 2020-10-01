**********************************************************************
*** Program DMhalo_table shows how to load a tabulated DM density profile
*** (radius vs density). In this example, the profile is spherical, and 
*** interpolated with cubic splines or linear interpolation, in linear or
*** logarithmic scales. The program also contains examples about how to call
*** rates for Milky Way cosmic ray fluxes for (such) temporary profiles   
*** Quick reading guide:
***
*** line 58 ff  : define (NFW-like) halo model based on table from file
***               + sample check against standard implemented (NFW) halo
*** line 122 ff : define (Burkert-like) halo model based on table 
***               + sample check against standard implemented (Burkert) halo
*** line 197 ff : detailed examples for how to compute CR fluxes for first 
***               numerical halo model
***               + comparison to standard implemented (NFW) halo
*** line 366 ff : Short version for how to compute CR fluxes for second 
***               numerical halo model
***               + comparison to standard implemented (Burkert) halo
***
**********************************************************************
      program DMhalo_table
      implicit none
ccc
      real*8 mwimp,sv,SI  ! generic wimp parameters
      integer pdgann
      logical selfconj      
ccc
      integer ii,jj,kk      
      real*8 radius,res1,res2,dsdmdsph
      real*8 res3,res4,tp,phiin,dspbdphidtaxi,dsdbdphidtaxi,
     & dsepdphidpaxi
      integer ntnfw,ntbur
      real*8 xnfw(1000),ynfw(1000),xbur(1000),ybur(1000)
ccc
      character*12 tag
      integer how,nt
      real*8 distance,rinn,rmax,x(1000),y(1000),rlow,rup
      integer testunit
      data testunit/98/
      character*100 scr
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
      SI=1.d-5                  ! WIMP-proton scattering (pb) [irrelevant here]

ccc      
      call dsgivemodel_generic_wimp(mwimp,selfconj,sv,pdgann,SI)
ccc
      write(*,*) 'Sample results in case of a WIMP of mass = ',mwimp
      write(*,*) 'annihilating into the b-bbar channel'

ccc
ccc load a table of values radius vs density; this particular example table 
ccc was created with the default MW NFW profile in DS, and with tag 'mwnfwdef'     
ccc      
      open(unit=testunit,file='mwnfwtab.dat',status='unknown',
     &     form='formatted')
      read(testunit,*) scr ! get rid of one line header
      ii=0
 10   ii=ii+1
      read(testunit,*,end=20) x(ii),y(ii)
      goto 10
 20   close(testunit)      
ccc
ccc specify the other parameters which are needed to define the profile
ccc      
      nt=ii-1 ! number of entries
      how=3   ! log(r) vs log(rho), cubic spline interpolation      
      rinn=1.d-5     ! other settings matching the definition '
      rmax=2.d2      ! of 'mwnfwdef
      distance=8.d0
ccc
ccc define a tag: to use the interpolation from the table the tag
ccc MUST contain the string 'num'
ccc
      tag='mwnum1'
ccc
ccc initialize the tabulated profile via the sample routine provided
ccc below, see the header for more details on the input parameters:
ccc      
      call sethalo_num(tag,distance,rinn,rmax,how,x,y,nt)
ccc
ccc save tables for later use:
ccc
      ntnfw=nt
      do ii=1,nt
        xnfw(ii)=x(ii)
        ynfw(ii)=y(ii)
      enddo  
ccc
ccc sample check of the tabulated profile versus the default MW NFW 
ccc profile:
ccc      
      write(*,*)
      write(*,*) ' sample check on default NFW Milky Way profile'
      write(*,1000) 'radius      '
     &             ,'param. prof.'
     &             ,'tabul. prof.'
     &             ,'% difference'
      rlow=0.1d-2
      rup=50.d0
      do ii=0,10
        radius=dexp(dlog(rlow)+(dlog(rup)-dlog(rlow))/10.d0*ii) 
        call dsdmdselect_halomodel('mwnfwdef')
        res1=dsdmdsph(radius)
        call dsdmdselect_halomodel('mwnum1')        
        res2=dsdmdsph(radius)
        write(*,1001) radius,res1,res2,dabs(res1-res2)/dabs(res1+res2)
      enddo
 1000 format(8(3x,a12))
 1001 format(10(1x,E14.6))      

ccc
ccc Let us now consider a second tabulated profile, based on the 
ccc default MW Burkert profile:
ccc      
      call dsdmdselect_halomodel('mwburdef')
ccc 
ccc This example demonstrates that the table {(x,y)} containing the
ccc profile does NOT have to be an ORDERED list (but must contain
ccc physical values, i.e. x(ii)>0 for all ii). We thus first generate
ccc a table with the 'wrong' order of entries:
ccc      
      rinn=1.d-5
      rmax=2.d2
      distance=8.d0
      jj=0
      do ii=200,0,-2
        radius=dexp(dlog(rinn)+(dlog(rmax)-dlog(rinn))/200*ii) 
        res1=dsdmdsph(radius)
        jj=jj+1
        x(jj)=radius
        y(jj)=res1
      enddo
      do ii=1,199,2
        radius=dexp(dlog(rinn)+(dlog(rmax)-dlog(rinn))/200*ii) 
        res1=dsdmdsph(radius)
        jj=jj+1
        x(jj)=radius
        y(jj)=res1
      enddo
      tag='mwnum2'
      how=3
      nt=jj
ccc
ccc and load it as a table to be interpolated:
ccc      
      call sethalo_num(tag,distance,rinn,rmax,how,x,y,nt)
ccc
ccc save tables for later use:
ccc
      ntbur=nt
      do ii=1,nt
        xbur(ii)=x(ii)
        ybur(ii)=y(ii)
      enddo  
ccc
ccc repeat the sample check of the tabulated profile
ccc      
      write(*,*)
      write(*,*) ' sample check on default Burkert Milky Way profile'
      write(*,1000) 'radius      '
     &             ,'param. prof.'
     &             ,'tabul. prof.'
     &             ,'% difference'
      rlow=0.1d-2
      rup=50.d0
      do ii=0,10
        radius=dexp(dlog(rlow)+(dlog(rup)-dlog(rlow))/10.d0*ii) 
        call dsdmdselect_halomodel('mwburdef')
        res1=dsdmdsph(radius)
        call dsdmdselect_halomodel('mwnum2')        
        res2=dsdmdsph(radius)
        write(*,1001) radius,res1,res2,dabs(res1-res2)/dabs(res1+res2)
      enddo
ccc
ccc NOTE:
ccc   1) you can decide case by case how many points you need,
ccc   however in the current setup 999 points is the maximum, change
ccc   in 'dsdmddrvrcom.h' if you need a larger number.      
ccc   2) tabulated profiles are treated as temporary profiles, there is
ccc   only one available at any given time and they overwrite previous
ccc   temporary profiles, even those defined in parametric form.      
ccc      


      write (*,*)
      write (*,*)
      write (*,*) 'Computing cosmic ray fluxes for the tabulated MW models',
     &     ' introduced so far'
      write (*,*) 'Start with the NFW profile'
      write (*,*)
ccc
ccc first of all you need to reload the temporary one:      
ccc
      rinn=1.d-5
      rmax=2.d2
      distance=8.d0
      tag='mwnum1'
      how=3
      call sethalo_num(tag,distance,rinn,rmax,how,xnfw,ynfw,ntnfw)
ccc
ccc the computation of the antiproton flux for a profile which is
ccc singular towards the Galactic center is in general CPU
ccc consuming; for a large scan over DM models and computations
ccc at several different energies a tabulation of the effective
ccc "confinement time" is convenient:      
ccc
      write (*,*)      
      write (*,*) 'Calculating antiproton fluxes...'
      tp=1.d0      ! kinetic energy GeV
      phiin=0.32d0 ! solar modulation parameter (forced-field) GV
      how=4        ! use previously computed table for default NFW
      write (*,*) 'Loading from the default (pre-computed) NFW table' 
      res1=dspbdphidtaxi(tp,phiin,how,'mwnfwdef')
      how=1        ! directly compute for default NFW
      write (*,*) 'Directly computing for the default NFW at',
     & ' tp+phiin = ',tp+phiin
      res2=dspbdphidtaxi(tp,phiin,how,'mwnfwdef')
ccc
ccc for a temporary profile you cannot store tables because temporary
ccc profiles are not uniquely defined; you can generate a table but it
ccc will not be saved to a file, nor can it be retrieved if you switch to 
ccc another temporary profile (even forcing the call with how=4, this
ccc parameter is internally changed to how=2, namely create the table
ccc but do not save it a file). In general it is not convenient to
ccc generate tabulations for temporary halo models; better to use
ccc direct calculations
ccc      
      how=1                     ! directly compute for default NFW
      write (*,*) 'Directly computing for the tabulated NFW at',
     & ' tp+phiin = ',tp+phiin
      res3=dspbdphidtaxi(tp,phiin,how,'mwnum1')
      write(*,1000) 'Kin. energy '
     &             ,'result      '
     &             ,'% diff. 1   '
     &             ,'% diff. 2   '
      write(*,1001) tp,res1,dabs(res1-res2)/dabs(res1+res2),
     & dabs(res2-res3)/dabs(res2+res3)   

      write (*,*)
      write (*,*) 'Calculating antiproton fluxes for a different ',
     & 'energy'
      tp=2.d0      ! kinetic energy GeV
      phiin=0.32d0 ! solar modulation parameter (forced-field) GV
      how=4        ! use preaviously computed table for default NFW
      write (*,*) 'Loading from the default NFW table' 
      res1=dspbdphidtaxi(tp,phiin,how,'mwnfwdef')
      how=1        ! directly compute for default NFW
      write (*,*) 'Directly computing for the default NFW at',
     & ' tp+phiin = ',tp+phiin
      res2=dspbdphidtaxi(tp,phiin,how,'mwnfwdef')
      how=1                     ! directly compute for tabulated NFW
      write (*,*) 'Directly computing for the tabulated NFW at',
     & ' tp+phiin = ',tp+phiin
      res3=dspbdphidtaxi(tp,phiin,how,'mwnum1')
      write(*,1000) 'Kin. energy '
     &             ,'result      '
     &             ,'% diff. 1   '
     &             ,'% diff. 2   '
      write(*,1001) tp,res1,dabs(res1-res2)/dabs(res1+res2),
     & dabs(res2-res3)/dabs(res2+res3)   
ccc
ccc Note however that in a scan over DM models connecting to the same
ccc values of "tp+phiin", the effective "confinement time" is not
ccc recomputed but values previously stored are used:  
ccc      
      write (*,*) 'Compute flux for tabulated NFW in a model scan'
      write(*,1000) 'WIMP mass   '
     &             ,'Kin. energy '
     &             ,'pbar flux   '
      do ii=1,10
        mwimp=100.d0+20.d0*ii ! WIMP mass in GeV
        call dsgivemodel_generic_wimp(mwimp,selfconj,sv,pdgann,SI)
        how=1
        res3=dspbdphidtaxi(tp,phiin,how,'mwnum1')
        write(*,1001) mwimp,tp,res3
      enddo
ccc
ccc back to WIMP mass of 100 GeV
ccc
      mwimp=100.d0 ! WIMP mass in GeV
      call dsgivemodel_generic_wimp(mwimp,selfconj,sv,pdgann,SI)
ccc 
ccc the case for antideuteron fluxes is specular:
ccc      
      write(*,*)    
      write(*,*) 'Calculating antideuteron fluxes...'
      tp=1.d0      ! kinetic energy GeV
      phiin=0.32d0 ! solar modulation parameter (forced-field) GV
      how=4        ! use preaviously computed table for default NFW
      write (*,*) 'Loading from the default NFW table' 
      res1=dsdbdphidtaxi(tp,phiin,how,'mwnfwdef')
      how=1        ! directly compute for default NFW
      write (*,*) 'Directly computing for the default NFW at',
     & ' tp+phiin = ',tp+phiin
      res2=dsdbdphidtaxi(tp,phiin,how,'mwnfwdef')
      how=1                     ! directly compute for tabulated NFW
      write (*,*) 'Directly computing for the tabulated NFW at',
     & ' tp+phiin = ',tp+phiin
      res3=dsdbdphidtaxi(tp,phiin,how,'mwnum1')
      write(*,1000) 'Kin. energy '
     &             ,'result      '
     &             ,'% diff. 1   '
     &             ,'% diff. 2   '
      write(*,1001) tp,res1,dabs(res1-res2)/dabs(res1+res2),
     & dabs(res2-res3)/dabs(res2+res3)   

      
      write (*,*)
      write (*,*) 'Calculating antideuteron fluxes for a different ',
     & 'energy'
      tp=2.d0      ! kinetic energy GeV
      phiin=0.32d0 ! solar modulation parameter (forced-field) GV
      how=4        ! use preaviously computed table for default NFW
      write (*,*) 'Loading from the default NFW table' 
      res1=dsdbdphidtaxi(tp,phiin,how,'mwnfwdef')
      how=1        ! directly compute for default NFW
      write (*,*) 'Directly computing for the default NFW at',
     & ' tp+phiin = ',tp+phiin
      res2=dsdbdphidtaxi(tp,phiin,how,'mwnfwdef')
      how=1                     ! directly compute for tabulated NFW
      write (*,*) 'Directly computing for the tabulated NFW at',
     & ' tp+phiin = ',tp+phiin
      res3=dsdbdphidtaxi(tp,phiin,how,'mwnum1')
      write(*,1000) 'Kin. energy '
     &             ,'result      '
     &             ,'% diff. 1   '
     &             ,'% diff. 2   '
      write(*,1001) tp,res1,dabs(res1-res2)/dabs(res1+res2),
     & dabs(res2-res3)/dabs(res2+res3)   
ccc
ccc the computation of the positron flux requires the tabulation of
ccc a Green function: only one table is available at any given time
ccc hence for temporary halo models (for which again tables cannot
ccc be saved to disc) the table is lost whenever you change profile
ccc      
      write (*,*)
      write (*,*) 'Calculating positron fluxes...'
      tp=1.d0
      phiin=0.32d0
      how=4     ! use preaviously computed table for default NFW
      write (*,*) 'Loading from the default NFW table' 
      res1=dsepdphidpaxi(tp,phiin,how,'mwnfwdef')
      how=2     ! use preaviously computed table for tabulated NFW
      write (*,*) 'Compute a table for the temporary tabulated NFW'
      res2=dsepdphidpaxi(tp,phiin,how,'mwnum1')
      write(*,1000) 'energy      '
     &             ,'result      '
     &             ,'% diff. 1   '
      write(*,1001) tp,res1,dabs(res1-res2)/dabs(res1+res2)




      write (*,*)
      write (*,*)
      write (*,*) 'Now switch to the Burkert profile'
      write (*,*)
      write (*,*) 'In this case all computations are much faster' 
ccc
ccc first of all you need to reload the temporary Burkert:      
ccc
      rinn=1.d-5
      rmax=2.d2
      distance=8.d0
      tag='mwnum2'
      how=3
      call sethalo_num(tag,distance,rinn,rmax,how,xbur,ybur,ntbur)

      write (*,*)
      write (*,*) 'Calculating antiproton fluxes...'
      tp=1.d0      ! kinetic energy GeV
      phiin=0.32d0 ! solar modulation parameter (forced-field) GV
      how=4        ! use preaviously computed table for default NFW
      write (*,*) 'Loading from the default Burkert table' 
      res1=dspbdphidtaxi(tp,phiin,how,'mwburdef')
      how=1        ! directly compute for default Burkert
      write (*,*) 'Directly computing for the default Burkert at',
     & ' tp+phiin = ',tp+phiin
      res2=dspbdphidtaxi(tp,phiin,how,'mwburdef')
      how=1                     ! directly compute for tabulated Burkert
      write (*,*) 'Directly computing for the tabulated Burkert at',
     & ' tp+phiin = ',tp+phiin
      res3=dspbdphidtaxi(tp,phiin,how,'mwnum2')
      write(*,1000) 'Kin. energy '
     &             ,'result      '
     &             ,'% diff. 1   '
     &             ,'% diff. 2   '
      write(*,1001) tp,res1,dabs(res1-res2)/dabs(res1+res2),
     & dabs(res2-res3)/dabs(res2+res3)   


      write(*,*)    
      write(*,*) 'Calculating antideuteron fluxes...'
      tp=1.d0      ! kinetic energy GeV
      phiin=0.32d0 ! solar modulation parameter (forced-field) GV
      how=4        ! use preaviously computed table for default NFW
      write (*,*) 'Loading from the default Burkert table' 
      res1=dsdbdphidtaxi(tp,phiin,how,'mwburdef')
      how=1        ! directly compute for default Burkert
      write (*,*) 'Directly computing for the default Burkert at',
     & ' tp+phiin = ',tp+phiin
      res2=dsdbdphidtaxi(tp,phiin,how,'mwburdef')
      how=1                   ! directly compute for tabulated Burkert
      write (*,*) 'Directly computing for the tabulated Burkert at',
     & ' tp+phiin = ',tp+phiin
      res3=dsdbdphidtaxi(tp,phiin,how,'mwnum2')
      write(*,1000) 'Kin. energy '
     &             ,'result      '
     &             ,'% diff. 1   '
     &             ,'% diff. 2   '
      write(*,1001) tp,res1,dabs(res1-res2)/dabs(res1+res2),
     & dabs(res2-res3)/dabs(res2+res3)   


      write (*,*)
      write (*,*) 'Calculating positron fluxes...'
      tp=1.d0
      phiin=0.32d0
      how=4     ! use preaviously computed table for default Burkert
      write (*,*) 'Loading from the default Burkert table' 
      res1=dsepdphidpaxi(tp,phiin,how,'mwburdef')
      how=2     ! use preaviously computed table for tabulated Burkert
      write (*,*) 'Compute a table for the temporary tab. Burkert'
      res2=dsepdphidpaxi(tp,phiin,how,'mwnum2')
      write(*,1000) 'energy      '
     &             ,'result      '
     &             ,'% diff. 1   '
      write(*,1001) tp,res1,dabs(res1-res2)/dabs(res1+res2)

      
      stop
      end ! main program


c_______________________________________________________________________
c_______________________ AUXILIARY ROUTINES ____________________________
c_______________________________________________________________________
  


      subroutine sethalo_num(tag,distance,rinn,rmax,how,x,y,nt)
c_______________________________________________________________________
c
c sample routine with proper parameter-label assignment for the
c intialization of a spherical model given as a table interpolation
c
c inputs:
c   - tag: label for the model, it must contain the string 'num'
c   - distance: distance from the center of the object
c   - rinn: inner truncation radius, i.e. rho(r)=rho(rinn) for r<rinn 
c   - rmax: maximum radius, i.e. rho(r)=0 for r>rmax
c   - how: interpolation options, i.e.:
c     how=1: tabulation radius vs density and spline interpolation     
c     how=2: tabulation radius vs density and linear interpolation     
c     how=3: tabulation log(r) vs log(rho) and spline interpolation     
c     how=4: tabulation log(r) vs log(rho) and linear interpolation     
c   - x: vector with values of the radius   
c   - y: vector with corresponding values of density   
c   - nt: number of entries for x & y
c distances to be given in kpc, densities in GeV cm^-3
c_______________________________________________________________________
      implicit none
      character*12 tag
      integer how,nt
      real*8 distance,rinn,rmax,x(nt),y(nt)
ccc      
      logical match
      integer ii,kk
      real*8 rein(2*nt+5)
      character*10 chin(2*nt+5),xlabel,ylabel
ccc
ccc check that tag contains 'num', otherwise the initialization will 
ccc fail:     
ccc
      call dslabcheck(12,tag,3,'num',match)
      if(.not.match) then
        write(*,*) 'DS: call to sethalo_num with tag = ',tag
        write(*,*) 'DS: not containg the string - num'
        write(*,*) 'DS: program stopped'
        stop
      endif  
ccc
ccc the order in which each of the parameters is assigned is arbitrary
ccc you just need to be careful to assign each predefined parameter
ccc label chin(i) to the parameter value rein(i), with i going from 1 
ccc to 5+2*nt:      
ccc
      chin(1)='objdist'
      rein(1)=distance
      chin(2)='radintr'
      rein(2)=rinn
      chin(3)='radouttr'
      rein(3)=rmax
      chin(4)='tabnum'
      rein(4)=dble(how)
ccc
ccc the labels for each entry in x & y are automatically generated with
ccc the subroutine dsdmdnumlab, the same called in dsdmddriver_num      
ccc      
      kk=4
      do ii=1,nt
        call dsdmdnumlab('x',ii,xlabel)
        call dsdmdnumlab('y',ii,ylabel)
        kk=kk+1
        chin(kk)=xlabel
        rein(kk)=x(ii)
        kk=kk+1
        chin(kk)=ylabel
        rein(kk)=y(ii)
      enddo
      call dsdmdnumlab('0',0,xlabel) ! storing the # of entries in x & y
      kk=kk+1
      chin(kk)=xlabel
      rein(kk)=dble(nt)
ccc initialization:
      call dsdmdset_halomodel(tag,kk,chin,rein)
      return
      end



      
