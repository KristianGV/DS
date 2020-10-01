**********************************************************************
*** program DMhalo_bypass shows how to provide your own hardcoded DM source,
*** bypassing the DM halo profile setting provided in DS.
***
*** WARNING: while simpler than DMhalo_new.f, this is still an example 
*** for 'expert users'! In particular, the drawback of this approach 
*** is that the system of automatic tabulation of quantities related 
*** to DM rates cannot be easily exploited
***
*** Concretely, we demonstrate below how the DS subroutine dsdmddriver
*** must be overwritten. The new hardcoded DM source is in this 
*** example provided by the function my_dm_sphsource further down 
*** (and for this specific demonstration example chosen to coincide 
*** with the NFW case). The example provided here assumes the hardcoded 
*** halo is of temporary kind and hence no tabulations are loaded.
***
*** For testing purposes, this program reads in a file DMhalo_bypass_sav.dat
*** that contains the same set of cosmic ray fluxes computed with the
*** NFW profile as provided in the standard DS setup (this file is
*** compute with DMhalo5.f)
**********************************************************************
      program DMhalo_bypass
      implicit none
ccc
      real*8 mwimp,sv,SI  ! generic wimp parameters
      integer pdgann
      logical selfconj      
ccc
      character*12 halotag
      real*8 egam,psi0,theta0,fluxgacdiff,dsgafluxsph,tpbar,pbflux,
     &  dspbdphidtaxi,tDbar,phiin,dsdbdphidtaxi,dbflux,eeplus,phiep,
     &  dsepdphidpaxi,mwrho0,dsdmdrho0mw
      integer istat
ccc
      integer testunit,testlevel
      character*200 savfile
      data savfile/'DMhalo_bypass_sav.dat'/
      data testunit/98/
ccc      
ccc Initialize DarkSUSY: in the simple example provided here also the
ccc initialization for the single hardcoded profile available is done
ccc through this call, containing an initialization call to dsdmddriver
ccc in the new version included below      
ccc
      call dsinit 
ccc
ccc compare against output computed in DMhalo5.f for the model that
ccc within the standard DS release is 'mwnfwdef'
ccc
c      testlevel=9
      testlevel=1
      if (testlevel.eq.9) then
        write(*,*) 'WRITING numerical values to DMtest savfile: ',
     &    savfile
	open(unit=testunit, file=savfile, status="unknown", err=3000)
      elseif (testlevel.eq.1.or.testlevel.eq.2) then
        write(*,*) 
     &    'COMPARING to numerical values from DMtest savfile: ',
     &    savfile
        open(unit=testunit, file=savfile, status="unknown", err=3000)
      else
        write(*,*) 'ERROR: unsupported value of testlevel=',testlevel 
      endif
ccc      
ccc sample generic WIMP model setup:
ccc      
      mwimp=100.d0              ! WIMP mass in GeV
      selfconj=.true.           ! self-conjugated WIMP
      sv=3.d-27                 ! WIMP pair annihilation cm^3/s
      pdgann=5                  ! b b-bar
      SI=1.d-5                  ! WIMP-proton scattering (pb)
ccc      
      call dsgivemodel_generic_wimp(mwimp,selfconj,sv,pdgann,SI)
ccc
      write(*,*) 'Sample results in case of a WIMP of mass = ',mwimp
      write(*,*) 'annihilating into the b-bbar channel'
ccc
ccc this example has only one DM halo model active, defined via the
ccc subroutine dsdmddriver given below, which bypasses the halo model
ccc structure in the DS release, replacing it by a single hardcoded
ccc model. in the specific example (and in any other unless you reset
ccc the procedure consistently) the halo model tag, while still
ccc formally needed to as an argument for the DS routines, is just a
ccc DUMMY variable!!!
ccc
      halotag='dummy'

      write(*,*)
      psi0=0.d0
      theta0=0.1d0 ! degree
      theta0=theta0*4.d0*datan(1.d0)/180.d0 ! rad
      egam=20.d0
      write(*,*) 'Calculating the gamma ray flux at E = ',egam
      fluxgacdiff=dsgafluxsph(egam,1,1.d0,halotag,psi0,theta0,istat)
      write(*,*) 'fluxgacdiff = ',fluxgacdiff,' ph/(cm^2 s GeV)'
      call testvalue(testunit,testlevel,'gamma flux',fluxgacdiff)

      write(*,*)
      tpbar=2.d0
      write(*,*) 'Calculating the antiproton flux at Tp = ',tpbar
      phiin=0.32d0 ! solar modulation parameter in GV
      pbflux=dspbdphidtaxi(tpbar,phiin,1,halotag) ! you cannot use
                               ! tables for the "confinement time"
      write(*,*) 'pbflux = ',pbflux,' GeV^-1 cm^-2 s^-1 sr^-1'
      call testvalue(testunit,testlevel,'pbar flux',pbflux)


      write(*,*)
      tDbar=1.d0
      write(*,*) 'Calculating the antideuteron flux at TD = ',tDbar
      phiin=0.32d0 ! solar modulation parameter in GV
      dbflux=dsdbdphidtaxi(tDbar,phiin,1,halotag) ! you cannot use
                               ! tables for the "confinement time"
      write(*,*) 'dbflux = ',dbflux,' GeV^-1 cm^-2 s^-1 sr^-1'
      call testvalue(testunit,testlevel,'Dbar flux',dbflux)

      write(*,*)
      eeplus=1.d0
      write(*,*) 'Calculating the positron flux at E = ',eeplus
      phiep=dsepdphidpaxi(eeplus,0.d0,4,halotag) ! use of tables
        ! necessary, reloaded whenever the halo model is changed
      write(*,*) 'phiep = ',phiep,' GeV^-1 cm^-2 s^-1 sr^-1'
      call testvalue(testunit,testlevel,'eplus flux',phiep)
      
      write(*,*)
      mwrho0=dsdmdrho0mw(halotag)
      write(*,*) 'Local halo density as imput for direct detection'
      write(*,*) 'rate routines, as well as neutrino flux from the'
      write(*,*) 'Sun and the Earth, rho0 = ',mwrho0
      call testvalue(testunit,testlevel,'MW rho0',mwrho0)
      
      write(*,*)
      write(*,*) 'The DMhalo_bypass test program ran successfully'
      call testvalue(testunit,testlevel,'final call',0.d0)
      stop
      
      
 3000 write(*,*) 'ERROR with file I/O: ',savfile     
      stop
      end ! main program

c_______________________________________________________________________
c_______________________ AUXILIARY ROUTINES ____________________________
c_______________________________________________________________________


      subroutine dsdmddriver(iwin,nrein,rein,nchin,chin,labin,reout)
c_______________________________________________________________________
c this is a sample driver to provide your own hardcoded DM source 
c function. This example works if you wish to use one single setting.
c
c inputs:
c   - iwin: integer setting the action of the driver; the minimum
c     setup MUST have at least the following entries:
c       idmdinit -> intialization called by dsinit; it can be even an
c         empty statement, in this example it contains initializations
c         of the global quantities associated to the DM source
c       idmdload -> loading of halo model from halo model database; 
c         not available here, just have a return statement        
c       idmdsousph -> spherical DM source (annihilation or decay)
c       idmdsouaxi -> axisymmetric DM source (annihilation or decay)
c         this is needed for MW dark matter sources only
c     in this minimal configuration you can you only DM rates without
c     tabulations; you cannot call, e.g., the DM density profile or any
c     of the input/output DM profile routines, such as, e.g.
c     dsdmdset_halomodel, dsdmdselect_halomodel or dsdmdprint_halomodel
c   - labin: dummy variable in the minimal setup
c   - rein(nrein): real*8 vector for different inputs according to
c      different iwin values; in this minimal setup you have only the
c      coordinates at which the DM source is computed
c   - chin(nchin): dummy character*10 vector in the minimal setup
c output:         
c   - reout: real*8 output for the different function linking as
c     specified by the input value iwin 
c_______________________________________________________________________
      implicit none
      include 'dsdmdcom.h' ! global variables for the DM source
      include 'dsdvcom.h' ! coordinates
      integer iwin
      character(*) labin
      integer nrein,nchin
      real*8 rein(nrein),reout
      character(*) chin(nchin)
ccc
      integer ihatype
      real*8 rsph,raxi,zaxi,my_dm_sphsource,my_dm_axisource
ccc
      if(iwin.eq.idmdinit) then
ccc
ccc all profiles will be temporary, no storing label, intialize
ccc dmdihalotag; 
ccc        
        dmdtmp=.true.
        dmdlabel='none'
        dmdihalotag=0 ! initialization, later incremented for different
                      ! dark matter source functions
ccc
ccc in this specific example in which one single hardcode profile is
ccc available you can set here those global variables that are needed by
ccc the DS library even before the first call to DM source function:
ccc the specific set of entries here match the default MW NFW in the
ccc DS release, i.e. the profile with tag mwnfwdef which is available
ccc in DS with the routine dsdmddriver which is being overwritten here
ccc NOTE: would you wish to implement a structure in which more than
ccc one model is available you need to reset model by model the full
ccc list from dmdihalotag to dmdrho0; note in particular that failing
ccc to change dmdihalotag would not force the reloading of some stored
ccc quantities and hence give wrong results!!!        
ccc
        dmdihalotag=dmdihalotag+1 ! identification number
        dmdmw=.true. ! this source CAN be used for MW CRs routines
c        dmdmw=.false. ! this source CANNOT be used for MW CRs routines
        dmdobjdist=8.d0 ! distance of the observer from the center
                        ! of the DM source in kpc
        dmdradintr=1.d-5 ! inner truncation radius in kpc, namely
                         ! source(r<dmdradintr) = source(r=dmdradintr) 
        dmdradouttr=200.d0  ! outer truncation radius in kpc, namely
                            ! source(r>dmdradouttr) = 0
        dmdrho0=0.3d0 ! local halo density for MW profile in GeV cm^-3
c        dmdrho0=0.d0 ! setting to zero if not MW profile
ccc
ccc a couple of cross checks that MW/non-MW profiles are consistently
ccc defined, remove them in case you know they are satisfied:
ccc      
        if(dmdmw) then
          if(dmdobjdist.gt.dmdradouttr) then
            write(*,*) 'DS: problem in dsdmddriver'
            write(*,*) 'DS: inconsistent definition of MW halo model:'
            write(*,*) 'DS: galactocentric distance = ',dmdobjdist
            write(*,*) 'DS: larger than halo size = ',dmdradouttr
            write(*,*) 'DS: program stopped'
            stop
          endif
        else
          if(dabs(dmdrho0).gt.1.d-16) then
            write(*,*) 'DS: inconsistent definition of dmdrho0 = '
     &        ,dmdrho0
            write(*,*) 'DS: for profile which has be declared as not MW'
            write(*,*) 'DS: program stopped'
            stop
          endif
        endif  
        return     
      endif
ccc
      if(iwin.eq.idmdload) then
ccc
ccc attempt to load a model form the halo database, which is not
ccc available here, just return
        return
      endif
ccc
ccc specify here if your model for the DM source is spherically 
ccc symmetric or axisymmetric (you could have also a homoeoid or
ccc a triaxial source, but in these cases no DM rates are available
ccc in the present DS release), this is not a global variable so you
ccc need to set it at every call to this routine:
ccc
      ihatype=itydvsph ! type: spherically symmetric -> r
c      ihatype=itydvaxi ! type: axisymmetric -> R,z 
ccc
ccc some generic linking, which you can eventually simplify, note
ccc however that in the current DS release, for MW profiles only,
ccc both idmdsousph and idmdsouaxi MUST be active:      
ccc
      if(iwin.eq.idmdsousph.and.ihatype.eq.itydvsph) then
ccc
ccc call to a spherically symmetric DM source function with spherical
ccc coordinates, extract the radius:
ccc
        rsph=rein(idvrsph) ! assumed to be in kpcc
        reout=my_dm_sphsource(rsph)
ccc
      elseif(iwin.eq.idmdsouaxi.and.ihatype.eq.itydvsph) then
ccc
ccc call to a spherically symmetric DM source function with axisymmetric
ccc coordinates, still extract the radius:
ccc
        rsph=dsqrt((rein(idvraxi))**2+(rein(idvzaxi))**2) ! kcp
        reout=my_dm_sphsource(rsph)
ccc
      elseif(iwin.eq.idmdsouaxi.and.ihatype.eq.itydvsph) then
ccc
ccc call to an axisymmetric DM source function with spherical
ccc coordinates, print an error and stop:
ccc
        write(*,*) 'DS: dsdmddriver with source of type itydvsph'
        write(*,*) 'DS: called with axisymmetric coordinates'
        write(*,*) 'DS: program stopped'
        stop
ccc
      elseif(iwin.eq.idmdsouaxi.and.ihatype.eq.itydvaxi) then
ccc
ccc call to an axisymmetric DM source function with axisymmetric
ccc coordinates, extract radial and vertical coordinate:
ccc
        raxi=rein(idvraxi)
        zaxi=rein(idvzaxi)
        reout=my_dm_axisource(raxi,zaxi)
ccc
ccc add here eventual other options you want to activate...
ccc        
      else
        reout=0.d0 ! dummy setting to avoid some compilers complaints
        write(*,*) 'DS: linking to dsdmddriver with iwin = ',iwin
        write(*,*) 'DS: out of the range of values which have been set'
        write(*,*) 'DS: program stopped'
        stop
      endif
      return
      end


ccc
ccc hardcode here your dark matter source; if you set ihatype=itydvsph
ccc activate the spherical function, if instead you set ihatype=itydvaxi
ccc activate the axisymmetric function
ccc      
      real*8 function my_dm_axisource(raxi,zaxi)
      implicit none
      real*8 raxi,zaxi
      my_dm_axisource=0.d0
      write(*,*) 'DS: linking to my_dm_axisource in a setup in which'
      write(*,*) 'DS: the source function has been defined as'
      write(*,*) 'DS: spherically symmetric'
      write(*,*) 'DS: program stopped'
      stop
      end
ccc
      real*8 function my_dm_sphsource(rin)
ccc
ccc this example is just the default MW NFW in the DS release
ccc      
      implicit none
      include 'dsdmdcom.h' ! passing some quantities initialized above
               ! as well as ksoupow, selecting annihilations or decays 
      real*8 rin
      real*8 rloc,rs,rho0,r0,x,x0
      rloc=rin
      if(rloc.lt.dmdradintr) rloc=dmdradintr ! inner truncation 
      if(rloc.gt.dmdradouttr) then   ! outer truncation
        my_dm_sphsource=0.d0
        return
      endif
      rho0=dmdrho0 ! GeV cm^{-3}
      r0=dmdobjdist ! kpc
ccc rs, the NFW scale radius, is the extra parameter still to be set
      rs=20.d0 ! kpc
      x=rloc/rs
      x0=r0/rs
      if(x.lt.1.d-16) x=1.d-16
      my_dm_sphsource=rho0*x0*(1.d0+x0)**2/x/(1.d0+x)**2 ! this is the
                                                  ! DM density profile
      if(ksoupow.eq.psdecay) return  ! for decaying DM, the DS  
                                     ! convention is that the DM source
                                     ! is just the DM density profile
      if(ksoupow.eq.psannihi) then
        my_dm_sphsource=my_dm_sphsource**2 ! for pair annihilating DM
          ! the DS convention is that the DM source is the square of
          ! the DM density profile; NOTE: no factor of 1/2 here!!!
          ! NOTE: in case of signals from a population of (unresolved)
          ! substructures, this has to be fixed appropriately
        return
      endif
ccc you shouldn't get here:      
      write(*,*) 'DS: problem in my_dm_sphsource'
      write(*,*) 'DS: ksoupow = ',ksoupow,' not properly initialized'
      write(*,*) 'DS: it should be either  = ',psannihi,' or = '
     &     ,psdecay
      write(*,*) 'DS: program stopped'
      stop
      end 
      
            
      subroutine testvalue(testunit,testlevel,msg,value)
      implicit none
      integer testunit,testlevel
      character*(*) msg
      real*8 value,tmp,eps,epscomp
      data eps/0.003/   !required accuracy for match with data file
ccc
      integer memory            ! to count total number of errors
      save memory
      data memory /0/
ccc
      if (msg.eq.'final call') then
        write(*,*) 'Total number of errors in DMhalo_bypass: ',memory
        return
      endif
ccc
      if (testlevel.eq.9) then
        write (testunit,*,err=30) value
      else
        read (testunit,*,err=30) tmp
        epscomp=abs((value-tmp)/tmp)
        if (epscomp.gt.eps) then
          memory=memory+1
          write(*,*)
          write(*,*) '====================================='
          write(*,*) 'ERROR in computation of ',msg,' !!!',
     &               '   (# ',memory,'in dstest)'
          write(*,*) 'calculated value:      ',value
          write(*,*) 'saved reference value: ',tmp
          write(*,*) '====================================='
          write(*,*)
        endif
      endif
      return
 30   write(*,*) 'testvalue: ERROR with file I/O!'     
      stop
      end !testvalue  
      
