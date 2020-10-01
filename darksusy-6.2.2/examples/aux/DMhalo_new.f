**********************************************************************
*** program DMhalo_new shows how to a define a new parametric DM density
*** profile and make it consistently available to the halo driver routines
*** and hence all DS functionality connected to halo density profiles.
***
*** WARNING: unlike the example files DMhalo_predef.f and DMhalo_table.f, 
*** this is an example for 'expert users'!
***
*** Concretely, we demonstrate how to still use the default dsdmddriver 
*** and just add one extra option, for the example of the Zhao 
*** ('alpha-beta-gamma') profile. For that purpose one needs to
*** define a new driver, namely dsdmddriver_abg specific for the Zhao
*** profile, and overwrite the subroutine dsdmddriver_choice      
**********************************************************************
      program DMhalo_new
      implicit none
ccc
      real*8 mwimp,sv,SI  ! generic wimp parameters
      integer pdgann
      logical selfconj      
ccc
      integer ii
      real*8 radius,res1,res2,dsdmdsph
      real*8 rho0,rs,al,be,ga,r0,rinn,rmax
      character*12 tag1,tag2
ccc
      character*12 tagtmp,tagvec(7)
      character*2 n2label(2)
      real*8 res(7,11)
      integer jj
      real*8 flux,dsgafluxsph,egev,diff,xi,psi0,alphas,degtorad
      integer istat
ccc      
ccc Initialize DarkSUSY
ccc
      call dsinit
ccc      
ccc sample generic WIMP model setup:
ccc      
      mwimp=100.d0              ! WIMP mass in GeV
      selfconj=.true.           ! self-conjugated WIMP
      sv=3.d-27                 ! WIMP pair annihilation cm^3/s
      pdgann=5                  ! b b-bar
      SI=1.d-5                  ! WIMP-proton scattering (pb)
                                ! irrelevant
ccc      
      call dsgivemodel_generic_wimp(mwimp,selfconj,sv,pdgann,SI)
ccc
      write(*,*) 'Sample results in case of a WIMP of mass = ',mwimp
      write(*,*) 'annihilating into the b-bbar channel'
ccc
ccc the scheme to implement to add a new parametric profile is
ccc essentially the same as for the other parametric profiles already 
ccc available in DS library:
ccc      
ccc 1) define a set of rules to identify the parameters of that
ccc    specific parametric profiles. For the spherical Zhao profile:
ccc
ccc     rho(r) = rhosabg/((r/rsabg)^gamma
ccc                       *(1+(r/rsabg)^alpha)^(beta-gamma)/alpha)
ccc
ccc    one possibility is:
ccc
ccc ******* start of the block to be uncommented and integrated in code
ccc   
ccc setting the tag and the corresponding integer identifier to select
ccc to this model; both of them have to be different from any of the
ccc others already in use in the code which are:
ccc   'nfw', 'bur', 'ein' & 'num' for the tag      
ccc   infwsuff=1, ibursuff=2, ieinsuff=3 & inumsuff=4 for the integer
ccc see the include file dsdmddrvrcom.h for details. While in principle
ccc just one of the two is needed, it is more practical to have one 
ccc tag as external identifier to initialize the model to be integrated
ccc in the halo model tag, and one integer for fast internal linker in
ccc the code runs:      
c      character(*) abgsuff ! tag to uniquely identify that profile
c      integer nchabgsuff ! number of characters for the abg tag 
c      integer iabgsuff ! integer number uniquely identifying the abgtag
c      parameter(abgsuff='abg',nchabgsuff=3,iabgsuff=5) 
ccc setting for each parameter a tag and an integer which is used as
ccc an index tracking the entry at which the value of that parameter is
ccc stored in the halo database when the parameter value is assigned.
ccc The order in which parameters are assigned is arbitrary, as well as
ccc of course the parameter tags which just need to be different between
ccc different parameters and from other tags which are associated to
ccc that model: e.g. 'rho0' is already taken for the local DM halo
ccc density (no authomatic check for label doubling in the present DS
ccc release), check dsdmddrvrcom.h for details. Also since these extra
ccc quantities are already associated to an integer, the assignment
ccc of integers below needs to be unique and different from those
ccc already taken, in the present DS release from iwhichobject to
ccc iradouttr, namely from 1 to 7; it also MUST be smaller or equal
ccc to nhaloindmax, equal to 12 in the current DS release, check
ccc dsdmddrvrcom.h for details. While nhaloindmax is kept minimal since
ccc memory allocations scale up accordindly, would you need more
ccc entries you need to change the hardcoded value of nhaloindmax in
ccc dsdmddrvrcom.h and recompile. 
c      integer irhosabg  ! refence density
c      character(*) chrhosabg
c      parameter(irhosabg=8,chrhosabg='rhosabg')
c      integer irsabg  ! scale length
c      character(*) chrsabg
c      parameter(irsabg=9,chrsabg='rsabg')
c      integer ialphaabg  ! alpha in the abg profile
c      character(*) chalphaabg
c      parameter(ialphaabg=10,chalphaabg='alphaabg')
c      integer ibetaabg  ! beta in the abg profile
c      character(*) chbetaabg
c      parameter(ibetaabg=11,chbetaabg='betaabg')
c      integer igammaabg  ! gamma in the abg profile
c      character(*) chgammaabg
c      parameter(igammaabg=12,chgammaabg='gammaabg')
ccc
ccc ******* end of the block to be uncommented and integrated in code
ccc
ccc 2) define a driver taking care of: a) storing of parameters and
ccc    eventual generation of labels for storing tabulated quantities,
ccc    such as, e.g., the antiproton "confinement time"; b) perform
ccc    the computation of the profile; c) print info on the currently
ccc    used profile - optional, a) and b) are compulsory. An example of
ccc    such driver is the subroutine dsdmddriver_abg below; the block
ccc    of definitions at 1) MUST be included in this subroutine
ccc
ccc 3) update the routine selecting among available DM halo profiles
ccc    options, namely the subroutine dsdmddriver_choice. Again an
ccc    example of such update is included below; the block of 
ccc    definitions at 1) MUST be included in this subroutine. The
ccc    provided routine MUST have exactly this name and exactly the same
ccc    input/output arguments, since including it in the main file or
ccc    linking to it when compiling this main file replaces the same 
ccc    name routine in the DS library!
ccc
ccc 4) properly initialize the profile, either in the main file or via 
ccc    some function; the subroutine sethalo_abgmw is an example on how
ccc    to do it, analogous to the subroutine sethalo_nfwmw in DMhalo_predef.f
ccc
ccc Since the NFW profile can be considered as a special case of the
ccc Zhao profile, to check that this setting works, we generate a model
ccc exactly matching the DS default NFW mwnfwdef:      
ccc      
      rho0=0.3d0 ! GeV cm^-3
      rs=20.d0   ! kpc
      r0=8.d0    ! kpc
      rinn=1.d-5 ! kpc    ! other settings matching the definition '
      rmax=2.d2  ! kpc    ! of 'mwnfwdef
      al=1.d0
      be=3.d0
      ga=1.d0
      tag1='mwabghdb1'
ccc
ccc the rules for halo model assignments are the same as those explained
ccc in DMhalo_predef.f, in particular the strings:
ccc     'mw' -> Milky Way,
ccc     'abg' -> Zhao profile, see the block of definitions above
ccc     'hdb' -> include this model in the DM halo database
ccc      
      call sethalo_abgmw(tag1,rho0,rs,al,be,ga,r0,rinn,rmax)

      tag2='mwnfwdef'
      write(*,*)
      write(*,*) 'Testing the Zhao profile against the default NFW,'
      write(*,*) 'having intitialized the two to be identical:'
      write(*,1000) 'radius      '
     &          ,'NFW profile '
     &          ,'Zhao profile'
     &          ,'% difference'
      do ii=0,20
        radius=dexp(dlog(rinn)+(dlog(rmax*0.99d0)-dlog(rinn))/20.d0*ii) 
        call dsdmdselect_halomodel(tag2)
        res1=dsdmdsph(radius)
        call dsdmdselect_halomodel(tag1)
        res2=dsdmdsph(radius)
        write(*,1001) radius,res1,res2,dabs(res1-res2)/dabs(res1+res2)
      enddo
 1000 format(8(3x,a12))
 1001 format(10(1x,E14.6))      
ccc
ccc after this initialization the profile can be used for any of the
ccc DM rates available in DS: e.g. gamma-ray angular profile for
ccc different gamma DM density profiles:
ccc      
ccc  rho(r) = rhos/(r/rs)^gamma/(1+(r/rs))^(3-gamma)
ccc 
      alphas=0.5d0  ! degrees, aperture of the instrument acceptance
      degtorad=4.d0*datan(1.d0)/180.d0 ! degrees to rad conversion
      alphas=alphas*degtorad ! rad
      egev=1.d0 ! energy at which the flux is computed [GeV]
      diff=1    ! differential flux
      xi=1.d0   ! rescaling factor
ccc      
      do jj=0,6
        ga=0.7d0+0.1d0*jj
        tagtmp='mwabg' 
ccc
ccc NOTE: this is a temporary halo model, there is no need to define
ccc it uniquely, however we do it for clearness in printing the output
ccc        
        write(n2label,1002) jj+7
 1002   format(i2)
        if(jj+7.lt.10) call dscharadd(tagtmp,'0')
        call dscharadd(tagtmp,n2label)
        call sethalo_abgmw(tagtmp,rho0,rs,al,be,ga,r0,rinn,rmax)
        do ii=0,10
          psi0=1.d0*ii ! degrees, direction of observation
          psi0=psi0*degtorad
          flux=dsgafluxsph(egev,diff,xi,tagtmp,psi0,alphas,istat)
ccc
ccc save the output:
ccc          
          res(jj+1,ii+1)=flux
        enddo
        tagvec(jj+1)=tagtmp
      enddo

      write(*,*)
      write(*,*) 'Gamma-ray angular profile for different gamma ',
     &  'DM density profiles'
      write(*,*)
      write(*,1000) 'angle (deg) ',(tagvec(jj),jj=1,7)
      do ii=1,11
        write(*,1001) 1.d0*(ii-1),(res(jj,ii),jj=1,7)
      enddo  
      
      stop
      end ! main program


c_______________________________________________________________________
c_______________________ AUXILIARY ROUTINES ____________________________
c_______________________________________________________________________
        



      subroutine sethalo_abgmw(tag,rho0,rs,al,be,ga,r0,rinn,rmax)
c_______________________________________________________________________
c
c sample routine with proper parameter-label assignment for the
c intialization of a spherical Zhao profile:
c
c   rho(r) = rho0*((r0/rs)^ga*(1+(r0/rs)^al)^(be-ga)/al)
c                /((r/rs)^ga*(1+(r/rs)^al)^(be-ga)/al)
c
c inputs:
c   - tag: label for the model, it must contain the string 'abg' & 'mw'
c   - rho0,rs,al,be,ga: the parameters in the formula above
c   - r0: local galactocentric distance
c   - rinn: inner truncation radius, i.e. rho(r)=rho(rinn) for r<rinn
c   - rmax: maximum radius, i.e. rho(r)=0 for r>rmax   
c distances to be given in kpc, rhos in GeV cm^-3
c_______________________________________________________________________
      implicit none
      character*12 tag
      real*8 rho0,rs,al,be,ga,r0,rinn,rmax
      logical match
      integer nin
      real*8 rein(8)
      character*10 chin(8)
ccc
ccc check that tag contains 'abg', otherwise the initialization
ccc will fail      
ccc
      call dslabcheck(12,tag,3,'abg',match)
      if(.not.match) then
        write(*,*) 'DS: call to sethalo_abgmw with tag = ',tag
        write(*,*) 'DS: not containg the string - abg'
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
ccc the order in which each of the 8 parameters is assigned is arbitrary
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
      chin(5)='rsabg'  ! NOTE: these are from the extra terms in
      rein(5)=rs       ! addition to dsdmddrvrcom.h
      chin(6)='alphaabg'
      rein(6)=al
      chin(7)='betaabg'
      rein(7)=be
      chin(8)='gammaabg'
      rein(8)=ga
      nin=8
ccc initialization:
      call dsdmdset_halomodel(tag,nin,chin,rein)
      return
      end

      
      

      subroutine dsdmddriver_choice(iwin,nrein,rein,nchin,chin,labin,
     &  reout)
c_______________________________________________________________________
c subroutine looping over available halo model drivers;
c NOTE: this version is obatined from the default version adding the
c Zhao profile
c_______________________________________________________________________
      implicit none
      include 'dsdmdcom.h'
      include 'dsdmddrvrcom.h'
      include 'dsdvcom.h'
ccc      
ccc *** included block see the explanation in the main above
ccc   
      character(*) abgsuff ! tag to uniquely identify that profile
      integer nchabgsuff ! number of characters for the abg tag 
      integer iabgsuff ! integer number uniquely identifying the abg tag
      parameter(abgsuff='abg',nchabgsuff=3,iabgsuff=5) 
      integer irhosabg  ! refence density
      character(*) chrhosabg
      parameter(irhosabg=8,chrhosabg='rhosabg')
      integer irsabg  ! scale length
      character(*) chrsabg
      parameter(irsabg=9,chrsabg='rsabg')
      integer ialphaabg  ! alpha in the abg profile
      character(*) chalphaabg
      parameter(ialphaabg=10,chalphaabg='alphaabg')
      integer ibetaabg  ! beta in the abg profile
      character(*) chbetaabg
      parameter(ibetaabg=11,chbetaabg='betaabg')
      integer igammaabg  ! gamma in the abg profile
      character(*) chgammaabg
      parameter(igammaabg=12,chgammaabg='gammaabg')
ccc      
ccc *** end of included block
ccc
      integer iwin
      character(*) labin
      integer nrein,nchin
      real*8 rein(nrein),reout
      character(*) chin(nchin)
      logical match
      integer ii
c
c if you are calling this function to create a profile you need first
c to identify which profile 
c      
      if(iwin.eq.idmdcreate) then
c
c is it a spherical nfw? if 'labin' contains the string 'nfwsuff' it
c is assumed to be so:
c        
        call dslabcheck(nchhalotag,labin,nchnfwsuff,nfwsuff,match)
        if(match) then
c store that it is a spherical nfw:         
          ihalotag(iwhichprof,ihalocur)=infwsuff
c store the coordinate type that the spherical nfw assumes:         
          ihalotag(itypeprof,ihalocur)=itydvsph
          goto 100
        endif
c
c is it a spherical burkert? if 'labin' contains the string 'bursuff' it
c is assumed to be so:
c        
        call dslabcheck(nchhalotag,labin,nchbursuff,bursuff,match)
        if(match) then
c store that it is a spherical burkert:         
          ihalotag(iwhichprof,ihalocur)=ibursuff
c store the coordinate type that the spherical burkert assumes:         
          ihalotag(itypeprof,ihalocur)=itydvsph
          goto 100
        endif
c
c is it a spherical einasto? if 'labin' contains the string 'einsuff' it
c is assumed to be so:
c        
        call dslabcheck(nchhalotag,labin,ncheinsuff,einsuff,match)
        if(match) then
c store that it is a spherical einasto:         
          ihalotag(iwhichprof,ihalocur)=ieinsuff
c store the coordinate type that the spherical einasto assumes:         
          ihalotag(itypeprof,ihalocur)=itydvsph
          goto 100
        endif
c
c is it a spherical model to be interpolated from a table? if 'labin'
c contains the string 'numsuff' it is assumed to be so:
c        
        call dslabcheck(nchhalotag,labin,nchnumsuff,numsuff,match)
        if(match) then
c store that it is a radially tabulated profile:         
          ihalotag(iwhichprof,ihalocur)=inumsuff
c store the coordinate type that this radial profile assumes:         
          ihalotag(itypeprof,ihalocur)=itydvsph
          goto 100
        endif
ccc NOTE: new case added here:
c
c is it a spherical Zhao profile? if 'labin' contains the string
c 'abgsuff' it is assumed to be so:
c        
        call dslabcheck(nchhalotag,labin,nchabgsuff,abgsuff,match)
        if(match) then
c store that it is a spherical nfw:         
          ihalotag(iwhichprof,ihalocur)=iabgsuff
c store the coordinate type that the spherical Zhao profile assumes:
          ihalotag(itypeprof,ihalocur)=itydvsph
          goto 100
        endif
ccc
ccc you called this function with a label for which we cannot find a 
ccc match a message is written and the program is stopped      
ccc      
        write(*,*) 'DS: connecting to dsdmddriver_choice with option'
        write(*,*) 'DS: idmdcreate and labin = ',labin
        write(*,*) 'DS: not containing any of the defined suffices'
        write(*,*) 'DS: namely:',nfwsuff,' & ',bursuff,' & ',einsuff
        write(*,*) 'DS: program stopped'
        stop
 100    continue
      endif
ccc
      if(ihalotag(iwhichprof,ihalocur).eq.infwsuff) then
        call dsdmddriver_nfw(iwin,nrein,rein,nchin,chin,labin,reout)
      elseif(ihalotag(iwhichprof,ihalocur).eq.ibursuff) then
        call dsdmddriver_bur(iwin,nrein,rein,nchin,chin,labin,reout)
      elseif(ihalotag(iwhichprof,ihalocur).eq.ieinsuff) then
        call dsdmddriver_ein(iwin,nrein,rein,nchin,chin,labin,reout)
ccc NOTE: new case added here:
      elseif(ihalotag(iwhichprof,ihalocur).eq.iabgsuff) then
        call dsdmddriver_abg(iwin,nrein,rein,nchin,chin,labin,reout)
      endif
ccc
      return
      end
      


      

      subroutine dsdmddriver_abg(iwin,nrein,rein,nchin,chin,labin,reout)
c_______________________________________________________________________
c this is the sample driver to initialize, print or connect corresponding
c functions in case of a dark matter density profile described by the
c spherical Zhao profile:
c
c     rho(r) = rhosabg/((r/rsabg)^gamma
c                       *(1+(r/rsabg)^alpha)^(beta-gamma)/alpha)
c
c inputs:
c   - iwin: integer setting the action of the driver; it can be equal
c     to:
c       idmdcreate -> profile creation is completed here (first settings
c         are in the dsdmddriver subroutine)      
c       idmdparprint -> the current halo profile is printed
c         depending on the input labin, see below      
c       idmddensph -> the driver connects to the spherical density profile
c   - labin: external label of the current halo model: you got here in
c       case of labin containing the string 'abg'; if it contains also
c       the string 'mw' the model is intepreted as referring to the
c       Milky Way, you can define it assigning the local halo density
c       rather than the reference density 'rhosabg' and an internal
c       label for tabulation files is automatically created      
c   - rein(nrein) real*8 vector for different inputs according to
c      different iwin values       
c   - chin(nchin) character*10 vector for input parameter labels
c output:         
c   - reout: real*8 output for the different function connecting as
c     specified by the input value iwin
c NOTE: this function is written in analogy to dsdmddriver_nfw, the
c parts which are changed compared to that are enlighted by NOTE
c_______________________________________________________________________
      implicit none
      include 'dsdmdcom.h'
      include 'dsdmddrvrcom.h'
      include 'dsdvcom.h'
      include 'dsmpconst.h'
ccc      
ccc *** included block see the explanation in the main above
ccc   
      character(*) abgsuff ! tag to uniquely identify that profile
      integer nchabgsuff ! number of characters for the abg tag 
      integer iabgsuff ! integer number uniquely identifying the abg tag
      parameter(abgsuff='abg',nchabgsuff=3,iabgsuff=5) 
      integer irhosabg  ! refence density
      character(*) chrhosabg
      parameter(irhosabg=8,chrhosabg='rhosabg')
      integer irsabg  ! scale length
      character(*) chrsabg
      parameter(irsabg=9,chrsabg='rsabg')
      integer ialphaabg  ! alpha in the abg profile
      character(*) chalphaabg
      parameter(ialphaabg=10,chalphaabg='alphaabg')
      integer ibetaabg  ! beta in the abg profile
      character(*) chbetaabg
      parameter(ibetaabg=11,chbetaabg='betaabg')
      integer igammaabg  ! gamma in the abg profile
      character(*) chgammaabg
      parameter(igammaabg=12,chgammaabg='gammaabg')
ccc      
ccc *** end of included block
ccc
      integer iwin
      character(*) labin
      integer nrein,nchin
      real*8 rein(nrein),reout
      character(*) chin(nchin)
      character*10 labelout
      real*8 dsdmdsetlabel
      external dsdmdsetlabel
ccc      
      integer ii,jj,indexout,ndmdstart
      real*8 rs,rhos,alpha,beta,gamma,dsdmdprof_abg ! NOTE
      real*8 r,x,radouttr,radintr,normx,rho0,units
      logical mwlab,match,setrho
ccc
      if(iwin.eq.idmdcreate.or.iwin.eq.idmdparprint) then
c
c set first the list of parameters which needs to be associated to
c an internal tag (via dsdmdsetlabel.f & dsdmdgetlabel.f):        
c the reference density needs to be the first entry since in case of
c Milky Way profiles you can set instead the value of the density at
c the local position in the Galaxy, see below.
c        
ccc
ccc NOTE: the following is specific for the Zhao profile
ccc        
        idmdlab(1)=irhosabg
        chdmdlab(1)=chrhosabg
        idmdlab(2)=irsabg
        chdmdlab(2)=chrsabg
        idmdlab(3)=ialphaabg
        chdmdlab(3)=chalphaabg
        idmdlab(4)=ibetaabg
        chdmdlab(4)=chbetaabg
        idmdlab(5)=igammaabg
        chdmdlab(5)=chgammaabg
        idmdlab(6)=iradintr
        chdmdlab(6)=chradintr
        idmdlab(7)=iradouttr
        chdmdlab(7)=chradouttr
        ndmdlab=7  ! this number must be equal to the number of 
                   ! entries so far
        sufflab=abgsuff
c
c add then extra parameters to be also stored for this halo profile:
c
        idmdlab(8)=iobjdist
        chdmdlab(8)=chobjdist
        ndmdlabtot=8 ! again this number must be equal to the number of 
                     ! entries so far
      endif
c
c then go thorugh the list of iwin options:
c      
      if(iwin.eq.idmdcreate) then
c
c load whether it is Milky Way or not
c        
        mwlab=.false.
        if(ihalotag(iwhichobject,ihalocur).eq.imwsuff) mwlab=.true. 
c        
c complete the creation of the dark matter density profile: if it is
c not the Milky Way load all label parameters, if it is the Milky Way
c skip rhosabg at first        
c
        ndmdstart=1
        if(mwlab) ndmdstart=2 
        do ii=ndmdstart,ndmdlabtot
          match=.false. 
          do jj=1,nchin
            if(chin(jj).eq.chdmdlab(ii)) then
              call dsdmdsetind(indexout,rein(jj))
              ihalotag(idmdlab(ii),ihalocur)=indexout
              match=.true.
            endif
          enddo
          if(.not.match) then
            write(*,*) 'DS: in dsdmddriver_abg missing initialization'
            write(*,*) 'DS: of the halo parameter : ',chdmdlab(ii)
            write(*,*) 'DS: when setting the halo model : ',labin
            write(*,*) 'DS: program stopped '
            stop
          endif
        enddo
        if(mwlab) then
c
c adding rhosabg, possibly in terms of the local halo density
c
          setrho=.false.
          do jj=1,nchin
            if(chin(jj).eq.chrho0) then      
              r=dmdparvec(ihalotag(iobjdist,ihalocur))
              rs=dmdparvec(ihalotag(irsabg,ihalocur))
              alpha=dmdparvec(ihalotag(ialphaabg,ihalocur)) ! NOTE
              beta=dmdparvec(ihalotag(ibetaabg,ihalocur))
              gamma=dmdparvec(ihalotag(igammaabg,ihalocur))
              x=r/rs
              normx=dsdmdprof_abg(x,alpha,beta,gamma) ! NOTE: this
                                         ! function is given below
              normx=rein(jj)/normx
              call dsdmdsetind(indexout,normx)
              ihalotag(irhosabg,ihalocur)=indexout
              call dsdmdsetind(indexout,rein(jj))
              ihalotag(irho0,ihalocur)=indexout
              setrho=.true.
            endif
          enddo
          do jj=1,nchin
            if(chin(jj).eq.chrhosabg) then
              if(setrho) then
                write(*,*) 'DS: inconsistency !!! In dsdmddriver_abg'
                write(*,*) 'DS: you are setting both rho0 and rhos'
                stop
              endif
              call dsdmdsetind(indexout,rein(jj))
              ihalotag(irhosabg,ihalocur)=indexout
              r=dmdparvec(ihalotag(iobjdist,ihalocur))
              rs=dmdparvec(ihalotag(irsabg,ihalocur))
              rhos=dmdparvec(ihalotag(irhosabg,ihalocur))
              alpha=dmdparvec(ihalotag(ialphaabg,ihalocur)) ! NOTE
              beta=dmdparvec(ihalotag(ibetaabg,ihalocur))
              gamma=dmdparvec(ihalotag(igammaabg,ihalocur))
              x=r/rs
              rho0=rhos*dsdmdprof_abg(x,alpha,beta,gamma)
              call dsdmdsetind(indexout,rho0)
              ihalotag(irho0,ihalocur)=indexout
              setrho=.true.
            endif  
          enddo
          if(.not.setrho) then
            write(*,*) 'DS: Error: In dsdmddriver_abg you are'
            write(*,*) 'DS: not setting neither rho0 nor rhos'
            write(*,*) 'DS: program stopped '
            stop
          endif
        else
c
c if this is not the Milky Way, rho0 is set to zero
c
          rho0=0.d0
          call dsdmdsetind(indexout,rho0)
          ihalotag(irho0,ihalocur)=indexout
        endif
c 
c if it is Milky Way and a non-temporarily-defined halo model load
c a label tag for CR green function tabulations:
c
        if(mwlab.and.(ihalocur.ne.ihalotmp)) then
          call dsdmdgetlabel(dsdmdsetlabel,labelout)
          halotabtag(ihalocur)=labelout
c
c NOTE: in the future you may need tabulations also for models
c different from the Milky Way, add them below  
c          
        else
          halotabtag(ihalocur)='none'
        endif  
        return
      endif
ccc
      if(iwin.eq.idmdparprint) then
        write(*,*)
        write(*,*) 'DS: the halo model with access label: '
     &    ,halotag(ihalocur)
        write(*,*) 'DS: has a corresponding authomatically generated'
        write(*,*) 'DS: label for tabulation files equal to: '
     &    ,halotabtag(ihalocur)
        write(*,*) 'DS: the list of parameters for this model is: '
        do ii=1,ndmdlabtot
           write(*,*) 'DS: par & value : ',chdmdlab(ii)
     &                ,dmdparvec(ihalotag(idmdlab(ii),ihalocur)) 
        enddo
        return
      endif
ccc
      if(iwin.eq.idmddensph) then
c        
c connect to spherical dark matter density profile with inner density
c cutoff and outer truncation radius:
c
        r=rein(idvrsph)
        radouttr=dmdparvec(ihalotag(iradouttr,ihalocur))
        if(r.gt.radouttr) then
          reout=0.d0
          return
        endif  
        radintr=dmdparvec(ihalotag(iradintr,ihalocur))
        r=max(r,radintr)
        rs=dmdparvec(ihalotag(irsabg,ihalocur))
        rhos=dmdparvec(ihalotag(irhosabg,ihalocur))
        alpha=dmdparvec(ihalotag(ialphaabg,ihalocur))  ! NOTE
        beta=dmdparvec(ihalotag(ibetaabg,ihalocur))
        gamma=dmdparvec(ihalotag(igammaabg,ihalocur))
        x=r/rs
        reout=rhos*dsdmdprof_abg(x,alpha,beta,gamma)  ! GeV cm^-3
        return
      endif
ccc
ccc you called this function with an iwin which is not set, an error
ccc message is written and the program is stopped      
ccc      
      write(*,*) 'DS: connecting to dsdmddriver_abg with iwin = ',iwin
      write(*,*) 'DS: out of the defined range of values'
      write(*,*) 'DS: program stopped'
      stop
      end

      
c_______________________________________________________________________
c  dark matter density profiles written in the form:
c
c    rho(r) = rhos * dsdmdprof(r/rs,alpha,beta,gamma)
c
c  namely dsdmdprof is the function g(x), e.g., in Eq. (10) of 
c  astro-ph/0207125.
c      
c  NOTE: version valid for a Zhao profile only      
c      
c  input:
c    x [1]
c    alpha [1]: exponents in abg profile
c    beta [1]
c    gamma [1]
c_______________________________________________________________________
      real*8 function dsdmdprof_abg(x,alpha,beta,gamma)
      implicit none
      real*8 x,alpha,beta,gamma
      if(x.lt.1.d-16) x=1.d-16
      dsdmdprof_abg=1.d0/(x**gamma
     &   *(1.d0+x**alpha)**((beta-gamma)/alpha))
      end
