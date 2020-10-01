*         -*- mode: fortran -*-

c
c internal variables within the dsdmddriver:      
c
      
ccc
ccc halo profile storing
ccc      
      integer ihalocur     ! index for the current halo tag
      integer ihalotmp     ! reserve one entry to a temporary halo model
                           ! which can be overloaded rather than stored
      integer nhalotag     ! number of entries currently stored
      integer nhalotagmax  ! maximum number of entries that can be stored
      parameter(nhalotagmax=100)
      integer nhaloindmax  ! maximum number of indices stored for each tag
      parameter(nhaloindmax=12)
      integer ihalotag(nhaloindmax,nhalotagmax) ! matrix of stored indices
      integer nchhalotag   ! maximum number of characters for halo tags
      parameter(nchhalotag=12) ! if you change this, change the declaration
        ! below for halotag to character*nhalotag. NOTE: you must keep
        ! the declaration for halotabtag fixed to character*10
      character*12 halotag(nhalotagmax) ! vector storing running halo tags
      character*10 halotabtag(nhalotagmax) ! vector storing running halo
        ! tabulation tag, automatically generated and not needed
        ! for all halotag 
      common/halodatabasecom/ihalotag,ihalocur,ihalotmp,nhalotag,
     &  halotag,halotabtag

ccc
ccc parameters storing:
ccc 
      integer ndmdpar ! number of parameters currently saved
      integer ndmdpartmp ! number of parameters currently saved for
                         ! temporary halo model (allowing overwriting)
      integer ndmdparmax ! maximum number of parameters that can be saved
      parameter(ndmdparmax=1000)
      real*8 dmdparvec(ndmdparmax) ! vector of parameters
      common/dsdmddrcomint/ndmdpar,ndmdpartmp
      common/dsdmddrcomre/dmdparvec
      
      save /halodatabasecom/,/dsdmddrcomint/,/dsdmddrcomre/
c
c label storing for parameters that should be set for any (Milky Way) halo
c the integer label selects the first entry in the ihalotag matrix above
c the character label selects the corresponding label for parameter loading
c
      integer iwhichobject   ! which kind of halo object: now only mw or not-mw
      parameter(iwhichobject=1)
      integer iwhichprof   ! specific halo driver to use
      parameter(iwhichprof=2)
      integer itypeprof  ! which kind of coordinates the halo driver assumes
      parameter(itypeprof=3)
      integer iobjdist ! distance from observer to center of the object
      character(*) chobjdist
      parameter(iobjdist=4,chobjdist='objdist')
      integer irho0 ! local halo density NOTE: Milky Way only
      character(*) chrho0
      parameter(irho0=5,chrho0='rho0')
      integer iradintr ! inner truncation radius for object profile
      character(*) chradintr
      parameter(iradintr=6,chradintr='radintr')
      integer iradouttr ! outer truncation radius for object profile
      character(*) chradouttr
      parameter(iradouttr=7,chradouttr='radouttr')

      integer ndmdlab ! number of tags for storing internal label tags
      integer ndmdlabtot ! total number of tags
      character*4 sufflab
      character*10 chdmdlab(nhaloindmax) ! vector storing tags for
                                         ! a given halo type
      integer idmdlab(nhaloindmax) ! vector storing tags for
         ! a given halo type (stored as integer to speed things up)
      common/dmdlabtagcom/idmdlab,ndmdlab,sufflab

      save /dmdlabtagcom/


ccc
ccc list of halo model suffixes to store a model into the database:
      integer nrepsuff
      parameter(nrepsuff=2)
      character*3 repsuff(nrepsuff)
      integer nchrepsuff(nrepsuff) ! number of characters for tmp suffix
      integer irepsuff
      parameter(repsuff=(/'def','hdb'/))
      parameter(nchrepsuff=(/3,3/))
      
ccc
ccc tags for objects requiring special care:
ccc   the Milky Way:      
      character(*) mwsuff
      integer nchmwsuff ! number of characters for MW suffix
      integer imwsuff
      parameter(mwsuff='mw',nchmwsuff=2,imwsuff=1)
c add below any other special case:


      
c specific for the different choice of halo profiles:
c
c spherical Navarro-Frenk-White profile, i.e.:
c
c     rho(r) = rhosnfw/((r/rsnfw)*(1+r/rsnfw)^2)
c      
      character(*) nfwsuff
      integer nchnfwsuff ! number of characters for nfw suffix
      integer infwsuff
      parameter(nfwsuff='nfw',nchnfwsuff=3,infwsuff=1)
      integer irhosnfw  ! refence density
      character(*) chrhosnfw
      parameter(irhosnfw=8,chrhosnfw='rhosnfw') ! NOTE: this is iradouttr+1
      integer irsnfw  ! scale length
      character(*) chrsnfw
      parameter(irsnfw=9,chrsnfw='rsnfw')

ccc default DS NFW:      
      integer nnfwdef
      parameter(nnfwdef=5)
      character*10 labnfwdef,chnfwdef(5)
      real*8 renfwdef(5)
      parameter(labnfwdef='mwnfwdef')
      parameter(chnfwdef=(/'rho0    ','radintr ','radouttr','rsnfw   ',
     & 'objdist '/))
      parameter(renfwdef=(/0.3d0,1.d-5,2.d2,2.d1,8.d0/)) 


c
c spherical Burkert profile, i.e.:
c
c     rho(r) = rhosbur/(1+(r/rsbur))/(1+(r/rsbur)^2)
c      
      character(*) bursuff
      integer nchbursuff ! number of characters for bur suffix
      integer ibursuff
      parameter(bursuff='bur',nchbursuff=3, ibursuff=2) ! NOTE ibursuff
                                  ! needs to be different from infwsuff
      integer irhosbur  ! refence density
      character(*) chrhosbur
      parameter(irhosbur=8,chrhosbur='rhosbur') ! NOTE: this is iradouttr+1
      integer irsbur  ! scale length
      character(*) chrsbur
      parameter(irsbur=9,chrsbur='rsbur')
      
ccc default DS Burkert:      
      integer nburdef
      parameter(nburdef=5)
      character*10 labburdef,chburdef(5)
      real*8 reburdef(5)
      parameter(labburdef='mwburdef')
      parameter(chburdef=(/'rho0    ','radintr ','radouttr','rsbur   ',
     & 'objdist '/))
      parameter(reburdef=(/0.3d0,1.d-5,2.d2,2.d1,8.d0/)) 

c
c spherical Einasto profile, i.e.:
c
c     rho(r) = rhosein*exp(-2/alphaein*((r/rsein)^alphaein-1))
c      
      character(*) einsuff
      integer ncheinsuff ! number of characters for ein suffix
      integer ieinsuff
      parameter(einsuff='ein',ncheinsuff=3, ieinsuff=3) ! NOTE ieinsuff
                     ! needs to be different from infwsuff and ibursuff
      integer irhosein  ! refence density
      character(*) chrhosein
      parameter(irhosein=8,chrhosein='rhosein') ! NOTE: this is iradouttr+1
      integer irsein  ! scale length
      character(*) chrsein
      parameter(irsein=9,chrsein='rsein')
      integer ialphaein  ! alpha parameter
      character(*) chalphaein
      parameter(ialphaein=10,chalphaein='alphaein')
      
ccc default DS Einasto:      
      integer neindef
      parameter(neindef=6)
      character*10 labeindef,cheindef(6)
      real*8 reeindef(6)
      parameter(labeindef='mweindef')
      parameter(cheindef=(/'rho0    ','radintr ','radouttr','rsein   ',
     & 'objdist ','alphaein'/))
      parameter(reeindef=(/0.3d0,1.d-5,2.d2,2.d1,8.d0,0.18d0/)) 
      
      
c
c spherical profile from a table of values (radius-density)
c      
      character(*) numsuff
      integer nchnumsuff ! number of characters for num suffix
      integer inumsuff
      parameter(numsuff='num',nchnumsuff=3, inumsuff=4) ! NOTE inumsuff
           ! needs to be different from infwsuff, ibursuff and ieinsuff
      integer itabnum  ! how one should tabulate the input values
      character(*) chtabnum
      parameter(itabnum=8,chtabnum='tabnum') ! NOTE: this is iradouttr+1
ccc
ccc some options for the tabulation:
ccc      
      integer tabnumspl ! radius & density spline tabulation
      integer tabnumlin ! radius & density linear tabulation
      integer tabnumlogspl ! log(radius) & log(density) spline tabulation
      integer tabnumloglin ! log(radius) & log(density) linear tabulation
      parameter(tabnumspl=1,tabnumlin=2,tabnumlogspl=3,tabnumloglin=4)
ccc      
ccc labels for number of entries of the table and entries themselves
ccc are not in a fixed number and are then generated in the subroutine
ccc dsdmdnumlab under dsdmddriver_num; fix here the maximum number of
ccc entries and the storing matrix:
ccc
      integer ntabdim ! maximum number of entries
      parameter(ntabdim=999) ! choose between 999, 9999 and 99999
      character*3 n3label ! to generate up to 999 labels
c      parameter(ntabdim=9999) ! choose between 999, 9999 and 99999
      character*4 n4label ! to generate up to 9999 labels per axis
c      parameter(ntabdim=99999) ! choose between 999, 9999 and 99999
      character*5 n5label ! to generate up to 99999 labels per axis
ccc common block to store the current tabulation
      integer nentries ! number of entris in current tabulation
      integer ntype ! number of entris in current tabulation
      integer tabnumtype ! tabulation option in current tabulation
      real*8 xtabn(ntabdim),ytabn(ntabdim),ytabn2(ntabdim)
      common/dsdmdintsphcom/xtabn,ytabn,ytabn2,nentries,tabnumtype
