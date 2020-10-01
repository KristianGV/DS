*         -*- mode: fortran -*-

c
c external interface to the dsdmddriver:      
c
ccc
ccc CRs and l.o.s.i. interface quantities:       
ccc
      logical dmdmw ! logical variable saying whether the model should 
        ! be used for Milky Way rates or not
      logical dmdtmp ! logical variable saying whether the model is a 
        ! temporary model or not
      integer dmdihalotag ! integer uniquely assigned to a model in the
        ! halo database
      integer dmdihalotagtmp    ! integer uniquely assigned to a temporary
        ! model, stored but to be overwritten in one entry of the
        ! database
      character*10 dmdlabel ! internal dmd label for tabulation files
      real*8 dmdobjdist ! distance from the center of dark matter object
        ! (the Sun galactocentric distance in case of the Milky Way)
      real*8 dmdradintr ! inner density truncation radius 
      real*8 dmdradouttr  ! outer radius at which the density is set to 0
      real*8 dmdrho0 ! if dmdmw=.true. this is the local halo density
      common/dmdinterlcom/dmdmw,dmdtmp
      common/dmdintericom/dmdihalotag,dmdihalotagtmp
      common/dmdinterrcom/dmdobjdist,dmdradintr,dmdradouttr,dmdrho0
      common/dmdinterccom/dmdlabel    
      save /dmdinterlcom/,/dmdintericom/,/dmdinterlcom/,/dmdinterccom/
      
ccc
ccc entries to link instructions/functions to the driver:      
ccc      
      integer idmdinit ! integer value which tells the driver to   
                 ! initialize a set of default dm density profiles
      integer idmdcreate ! integer value which tells the driver to   
                 ! initialize the dm density profile for given label tag 
      integer idmdload ! integer value which tells the driver to retrieve  
                 ! and store the integer code for a given label tag
      integer idmdparprint ! integer value which tells the driver to 
                 ! print a specific (set of) parameter (or all of them)      
      integer idmddensph ! integer value which tells the driver to link to 
        ! the spherical dark matter density profile [GeV cm^-3]
      integer idmddenaxi ! the same as above but for the axisymmetric dark
        ! matter density profile [GeV cm^-3]
      integer idmddentri ! the same as above but for the dark matter 
        ! density profile in cartesian coordinates [GeV cm^-3]
      integer idmdsousph ! integer value which tells the driver to link to 
        ! the spherical dark matter source profile, scaling with rho (dm 
        ! decay) or rho**2 (dm pair annihilation) + eventual substructure 
        ! effects, [GeV cm^-3] or [GeV^2 cm^-6]       
      integer idmdsouaxi ! the same as above but for the axisymmetric
        ! dark matter source profile, [GeV cm^-3] or [GeV^2 cm^-6]
      integer idmdsoutri ! the same as above but for the dark matter 
        ! source profile in cartesian coordinates, [GeV cm^-3] or 
        ! [GeV^2 cm^-6]      
      integer idmdmassph ! integer value which tells the driver to link to 
        ! the spherical dark matter mass profile [M_sun]
ccc
ccc internal arbitrary assignment:
ccc
      parameter(idmdcreate=0,idmdload=-1,idmdparprint=-2,idmdinit=-3,
     & idmddensph=1,idmddenaxi=2,idmddentri=3,idmdsousph=4,idmdsouaxi=5,
     & idmdsoutri=6,idmdmassph=7)

ccc
ccc power index for pair annihilations or decays
ccc      
      integer psannihi,psdecay
      parameter(psannihi=2,psdecay=1)
      integer ksoupow ! power selecting annihilations or decays 
      common/ksoupowcom/ksoupow      

      
ccc
ccc print tag
ccc      
      character*12 printalltag
      parameter(printalltag='printall')

ccc
ccc dummy entries for some of the calls to dsdmddriver:
ccc      
      integer ndummy
      parameter(ndummy=1)
      real*8 redummy(ndummy),redummy2
      character*1 chdummy(ndummy),labdummy
      common/dmddummy/ redummy,redummy2,chdummy,labdummy ! TB added...
