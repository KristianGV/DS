      subroutine dsdmddriver(iwin,nrein,rein,nchin,chin,labin,reout)
c_______________________________________________________________________
c this is a sample driver to initialize, print or link functions
c corresponding to the halo model. 
c
c inputs:
c   - iwin: integer setting the action of the driver; it can be equal
c     to:
c       idmdinit -> intialization of a set of default halo model: you
c         can access via the label:
c          - 'mwnfwdef': the default Milky Way NFW profile
c          - 'mwburdef': the default Milky Way Burkert profile
c          - 'mweindef': the default Milky Way Einasto profile
c       idmdcreate -> a profile is created (started here and completed
c         in the sub-drivers called in this subroutine); unless another
c         model is loaded via idmdload below, this is assumed as current
c         halo model        
c       idmdload -> the profile corresponding to 'labin' is searched
c         for in the halo model database and loaded as current model
c       idmdparprint -> the current halo profile is printed
c         depending on the input labin, see below      
c       idmddensph -> the driver links to the spherical density profile
c       idmdsousph -> the driver links to the spherical dark matter
c         source profile, scaling with rho (dm decay) or rho**2 (dm pair
c         annihilation) [GeV cm^-3] or [GeV^2 cm^-6]; no substructure 
c         component in this sample implementation.
c       idmdmassph -> the driver links to the spherical mass profile
c     since in this implemetation all halo models are spherically 
c       symmetric calling with iwin idmddenaxi (axisymmetric coordinates)
c       or idmddentri (cartesian coordinates) is just turned to the case
c       iwin=idmddensph after proper coordinate transformation; the same
c       applies to idmdsouaxi and idmdsoutri
c   - labin: external label of the current halo model; the character
c       string MUST contain substrings for the driver initialization:
c       - one of the following is needed to specify which profile to
c         link:
c          a) 'nfw' a spherical NFW profile is assumed 
c          b) 'bur' a spherical Burkert profile is assumed
c          c) 'ein' a spherical Einasto profile is assumed
c          d) 'num' a spherical profile is given as a table interpolation
c       - if it contains 'hdb' the halo model is stored as one entry in
c         halo database; if it contains 'def' it is a default DS model,
c         also stored in the database; if it does not contain any of the
c         it is assumed to be a temporary model, one for which you cannot
c         store tabulations. Only one temporary model is available at any
c         given time. A 'num' models needs to be temporary.       
c   - rein(nrein) real*8 vector for different inputs according to
c      different iwin values       
c   - chin(nchin) character*10 vector for input parameter labels
c output:         
c   - reout: real*8 output for the different function linking as
c     specified by the input value iwin 
c_______________________________________________________________________
      implicit none
      include 'dsdmdcom.h'
      include 'dsdmddrvrcom.h'
      include 'dsdvcom.h'
      integer iwin
      character(*) labin
      integer nrein,nchin
      real*8 rein(nrein),reout
      character(*) chin(nchin)
ccc
      logical match,matchloc
      integer ii,ihatype
      logical first
      data first/.true./
      save first
ccc
      if(first) then
	nhalotag=0
        ndmdpar=0
        ndmdpartmp=0
	first=.false.
      endif
ccc
      if(iwin.eq.idmdinit) then
c
c include the possibility of a temporary halo model which can be 
c overloaded rather than uniquely added to the database
c
        nhalotag=nhalotag+1  ! add a halo model to the database
        ihalotmp=nhalotag
        ndmdpar=nhaloindmax  ! the first nhaloindmax entries reserved
                             ! for this model
c
c since this kind of models is temporary, an integer tracks the reloading
c
        dmdihalotagtmp=nhalotagmax ! to be incremented each time the
                             ! temporary halo is replaced by a new one
c
c add the set of default models        
c        
c   - default Milky Way NFW profile, a spherical profile:         
c
        nhalotag=nhalotag+1  ! add a halo model to the database
        ihalocur=nhalotag    ! update internal model linking
        halotag(ihalocur)=labnfwdef  ! default label 
        ihalotag(iwhichobject,ihalocur)=imwsuff ! Milky Way
        ihalotag(iwhichprof,ihalocur)=infwsuff ! NFW
        ihalotag(itypeprof,ihalocur)=itydvsph  ! spherical
        call dsdmddriver_nfw(idmdcreate,nnfwdef,renfwdef,nnfwdef,
     &     chnfwdef,labnfwdef,reout) ! associate the default parameters
             ! which are set via a parameter statement in the
             ! dsdmddrvrcom.h common block
c        
c   - default Milky Way Burkert profile, a spherical profile:         
c
        nhalotag=nhalotag+1  ! add a halo model to the database
        ihalocur=nhalotag    ! update internal model linking
        halotag(ihalocur)=labburdef  ! default label 
        ihalotag(iwhichobject,ihalocur)=imwsuff ! Milky Way
        ihalotag(iwhichprof,ihalocur)=ibursuff ! BUR
        ihalotag(itypeprof,ihalocur)=itydvsph  ! spherical
        call dsdmddriver_bur(idmdcreate,nburdef,reburdef,nburdef,
     &     chburdef,labburdef,reout) ! associate the default parameters
             ! which are set via a parameter statement in the
             ! dsdmddrvrcom.h common block
c        
c   - default Milky Way Einasto profile, a spherical profile:         
c
        nhalotag=nhalotag+1  ! add a halo model to the database
        ihalocur=nhalotag    ! update internal model linking
        halotag(ihalocur)=labeindef  ! default label 
        ihalotag(iwhichobject,ihalocur)=imwsuff ! Milky Way
        ihalotag(iwhichprof,ihalocur)=ieinsuff ! EIN
        ihalotag(itypeprof,ihalocur)=itydvsph  ! spherical
        call dsdmddriver_ein(idmdcreate,neindef,reeindef,neindef,
     &     cheindef,labeindef,reout) ! associate the default parameters
             ! which are set via a parameter statement in the
             ! dsdmddrvrcom.h common block
c
c set the NFW as current halo profile:
c        
        do ii=1,nhalotag
          if(halotag(ii).eq.labnfwdef) then
c
c CRs and l.o.s.i. interface quantities:       
c
            call dsdmddriver_if(ii)
          endif
        enddo
        return
      endif  
ccc
      if(iwin.eq.idmdcreate) then
c
c check whether you are creating a halo model to store in the database:
c        
        match=.false.
        do ii=1,nrepsuff
          call dslabcheck(nchhalotag,labin,nchrepsuff(ii),repsuff(ii),
     &      matchloc)
          if(matchloc) match=.true.
        enddo
c
c if not just re-initialize the temporary entry: 
c        
        if(.not.match) then
          ihalocur=ihalotmp ! set tmp as current profile
          ndmdpartmp=0  ! overwrite on parameters preaviusly stored
          dmdihalotagtmp=dmdihalotagtmp+1 ! unique integer labelling
                        ! this specific tmp halo model
          goto 100
        endif  
c
c a model given as a table interpolation needs to be a temporary entry: 
c        
        call dslabcheck(nchhalotag,labin,nchnumsuff,numsuff,match)
        if(match) then
          write(*,*) 'DS: call to dsdmddriver to create a halo model'
          write(*,*) 'DS: which you are trying to store although it is'
          write(*,*) 'DS: defined as a table interpolation; define it'
          write(*,*) 'DS: as temporary instead'
          write(*,*) 'DS: current label = ',labin
          write(*,*) 'DS: program stopped'
          stop
        endif
c
c otherwise:        
c check whether you are trying to create a halo with a tag that has
c already been assigned to previously set halo
c        
        do ii=1,nhalotag
          if(halotag(ii).eq.labin) then
            write(*,*) 'DS: call to dsdmddriver to create a halo model'
            write(*,*) 'DS: with a label that has already been assigned'
            write(*,*) 'DS: current label = ',labin
            write(*,*) 'DS: model number: ',ii,' with the same label = '
     &        ,halotag(ii)
            write(*,*) 'DS: program stopped'
            stop
          endif
        enddo  
c
c start with the creation of the dark matter density profile:
c
c add profile to the database:
        nhalotag=nhalotag+1 
c set it as current profile
        ihalocur=nhalotag
 100    continue
c store the user identification label: 
        halotag(ihalocur)=labin
c
c check whether you are creating a profile that applies to the Milky Way:
c to allow for another special cases which might be needed in the future 
c we use here a integer variable rather than a logical variable:
c        
        call dslabcheck(nchhalotag,labin,nchmwsuff,mwsuff,match)
        if(match) then
          ihalotag(iwhichobject,ihalocur)=imwsuff
        else  ! eventually add here other options:
          ihalotag(iwhichobject,ihalocur)=0
        endif
c
c check which kind of profile you want to link and complete the
c parameter assignment for this profile:        
c
        call dsdmddriver_choice(idmdcreate,nrein,rein,nchin,chin,labin,
     &      reout)
c        
c complete the initialization of the profile with          
c CRs and l.o.s.i. interface quantities:
c        
        call dsdmddriver_if(ihalocur)
        return
      endif
ccc
      if(iwin.eq.idmdload) then
c
c check whether you are loading a halo model which is supposed be to
c stored in the database:
c        
        match=.false.
        do ii=1,nrepsuff
          call dslabcheck(nchhalotag,labin,nchrepsuff(ii),repsuff(ii),
     &      matchloc)
          if(matchloc) match=.true.
        enddo
c
c if not just link to the last the temporary entry: 
c        
        if(.not.match) then
c
c CRs and l.o.s.i. interface quantities:       
c
          call dsdmddriver_if(ihalotmp)
          return
        endif
c
c otherwise find the model in the database:       
c        
        do ii=1,nhalotag
          if(halotag(ii).eq.labin) then
c
c CRs and l.o.s.i. interface quantities:       
c
            call dsdmddriver_if(ii)
            return
          endif
        enddo
        write(*,*) 'DS: attempt to load with dsdmddrive the label: '
     &     ,labin 
        write(*,*) 'DS: which has no match in the halo database'
        write(*,*) 'DS: program stopped'
        stop
      endif
ccc
      if(iwin.eq.idmdparprint) then
        if(labin.eq.printalltag) then
c
c force all models in the database to be printed and skip the
c check below:
c          
          match=.true.
          goto 200
        endif  
c
c check whether you are loading a halo model which is supposed be to
c stored in the database:
c        
        match=.false.
        do ii=1,nrepsuff
          call dslabcheck(nchhalotag,labin,nchrepsuff(ii),repsuff(ii),
     &      matchloc)
          if(matchloc) match=.true.
        enddo
        if(.not.match) then
          write(*,*)
          write(*,*) 'DS: the halo model with label = ',labin
          write(*,*) 'DS: was not stored in the database'
          write(*,*) 'DS: the latest temporary model is printed'
          write(*,*) 'DS: instead; the label priuted here for such'
          write(*,*) 'DS: model does not carry information'
          ihalocur=ihalotmp
          call dsdmddriver_choice(idmdparprint,nrein,rein,nchin,chin,
     &      labin,reout)
          return
        endif
c
c print the dark matter model corresponding to labin or all models:
c        
        match=.false.
 200    continue
        do ii=1,nhalotag
          if(halotag(ii).eq.labin.or.match) then
            ihalocur=ii
            call dsdmddriver_choice(idmdparprint,nrein,rein,nchin,chin,
     &      labin,reout)
          endif
        enddo
        return
      endif
ccc
      ihatype=ihalotag(itypeprof,ihalocur)
ccc
      if((iwin.eq.idmddensph.or.iwin.eq.idmdsousph)
     &   .and.ihatype.eq.itydvsph) then
c        
c link to spherical dark matter density profile:
c
        call dsdmddriver_choice(idmddensph,nrein,rein,nchin,chin,labin,
     &      reout)
        if(iwin.eq.idmdsousph) then
          if(ksoupow.eq.1) reout=reout ! decay, no substructures
          if(ksoupow.eq.2) reout=reout**2 ! annih., no substructures
        endif  
        return
      endif
ccc
      if((iwin.eq.idmddenaxi.or.iwin.eq.idmdsouaxi)
     &   .and.ihatype.eq.itydvsph) then
c        
c link to spherical dark matter density profile with axisymmetric
c coordinates, call a coordinate conversion first:
c
        call dsdvmatch(itydvaxi,nrein,rein,itydvsph,ntysph,dvsph)        
        call dsdmddriver_choice(idmddensph,ntysph,dvsph,nchin,chin,
     &       labin,reout)
        if(iwin.eq.idmdsouaxi) then
          if(ksoupow.eq.1) reout=reout ! decay, no substructures
          if(ksoupow.eq.2) reout=reout**2 ! annih., no substructures
        endif  
        return
      endif
ccc
      if((iwin.eq.idmddentri.or.iwin.eq.idmdsoutri)
     &   .and.ihatype.eq.itydvsph) then
c        
c link to spherical dark matter density profile with cartesian 
c coordinates, call a coordinate conversion first:
c
        call dsdvmatch(itydvtri,nrein,rein,itydvsph,ntysph,dvsph)        
        call dsdmddriver_choice(idmddensph,ntysph,dvsph,nchin,chin,
     &       labin,reout)
        if(iwin.eq.idmdsoutri) then
          if(ksoupow.eq.1) reout=reout ! decay, no substructures
          if(ksoupow.eq.2) reout=reout**2 ! annih., no substructures
        endif  
        return
      endif
ccc
      if(iwin.eq.idmdmassph.and.ihatype.eq.itydvsph) then
c        
c link to spherical dark matter mass profile:
c
        call dsdmddriver_choice(idmdmassph,nrein,rein,nchin,chin,labin,
     &      reout)
        return
      endif
ccc
ccc you called this function with an iwin for which we cannot find a 
ccc match, a message is written and the program is stopped      
ccc      
      write(*,*) 'DS: linking to dsdmddriver with iwin = ',iwin
      write(*,*) 'DS: out of the range of values which have been set'
      write(*,*) 'DS: program stopped'
      stop
      end



      
      subroutine dsdmddriver_if(ihaloif)
c
c CRs and l.o.s.i. interface quantities:       
c
      implicit none
      include 'dsdmdcom.h'
      include 'dsdmddrvrcom.h'
      integer ihaloif
      ihalocur=ihaloif
      dmdihalotag=ihalocur
      dmdtmp=.false.
      if(dmdihalotag.eq.ihalotmp) then
        dmdihalotag=dmdihalotagtmp
        dmdtmp=.true.
      endif  
      dmdmw=.false.
      if(ihalotag(iwhichobject,ihalocur).eq.imwsuff) dmdmw=.true. 
      dmdlabel=halotabtag(ihalocur)
      dmdobjdist=dmdparvec(ihalotag(iobjdist,ihalocur))
      dmdradintr=dmdparvec(ihalotag(iradintr,ihalocur))
      dmdradouttr=dmdparvec(ihalotag(iradouttr,ihalocur))
      dmdrho0=dmdparvec(ihalotag(irho0,ihalocur))
c
c cross check that MW profiles are consistently defined
c      
      if(dmdmw) then
        if(dmdobjdist.gt.dmdradouttr) then
          write(*,*) 'DS: inconsistent definition of MW halo model:'
          write(*,*) 'DS: for model with tag : ',dmdlabel
          write(*,*) 'DS: galactocentric distance = ',dmdobjdist
          write(*,*) 'DS: larger than halo size = ',dmdradouttr
          write(*,*) 'DS: program stopped'
          stop
        endif
      endif  
      return
      end
