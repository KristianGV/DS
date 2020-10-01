**********************************************************************
*** program DMhalo_bypass_prep.f generates the test file needed for 
*** DMhalo_bypass.f. Concretely, it computes the same DM rates as in 
*** DMhalo_bypass.f -- but for the DS default halo model 'mwnfwdef' 
*** provided with the DS release (rather than for the hardcoded DM 
*** source as provided in DMhalo_bypass.f).
*** NOTE: we need to write a separate main
*** file since the constructions are alternative one to the other.
**********************************************************************
      program DMhalo_bypass_prep
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
ccc Initialize DarkSUSY:
ccc
      call dsinit
ccc
ccc save outputs to check them in DMhalo_bypass.f
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
ccc subroutine dsdmddriver given below, which bypass the halo model
ccc structure in the DS release, replacing it by a single hardcoded
ccc model. in the specific example (and in any other unless you reset
ccc the procedure consistently) the halo model tag, while still
ccc formally needed to as an argument for the DS routines, is just a
ccc DUMMY variable!!!
ccc
      halotag='mwnfwdef'

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
      pbflux=dspbdphidtaxi(tpbar,phiin,4,halotag)
      write(*,*) 'pbflux = ',pbflux,' GeV^-1 cm^-2 s^-1 sr^-1'
      call testvalue(testunit,testlevel,'pbar flux',pbflux)


      write(*,*)
      tDbar=1.d0
      write(*,*) 'Calculating the antideuteron flux at TD = ',tDbar
      phiin=0.32d0 ! solar modulation parameter in GV
      dbflux=dsdbdphidtaxi(tDbar,phiin,4,halotag) ! you cannot use
                               ! tables for the "confinement time"
      write(*,*) 'dbflux = ',dbflux,' GeV^-1 cm^-2 s^-1 sr^-1'
      call testvalue(testunit,testlevel,'Dbar flux',dbflux)

      write(*,*)
      eeplus=1.d0
      write(*,*) 'Calculating the positron flux at E = ',eeplus
      phiep=dsepdphidpaxi(eeplus,0.d0,4,halotag)
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
        write(*,*) 'Total number of errors in dstest: ',memory
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
      
