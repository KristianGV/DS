      subroutine dslosidriver(iwin,nrein,rein,reout)
c_______________________________________________________________________
c this is a sample driver to initialize, tabulate or link functions
c needed for line of sign intergration of different dark matter
c emissivities, such as synchrotron, inverse Compton or Bremstrahlung. 
c
c this version is actually empty, since non of the above is included in
c the present release of the code.
c inputs:
c   - iwin: integer setting the action of the driver; ilosigasph is the
c        only one currently available
c   - rein(nrein) real*8 vector for different inputs according to
c      different iwin values       
c output:         
c   - reout: real*8 output for the different function linking as
c     specified by the input value iwin 
c_______________________________________________________________________
      implicit none
      include 'dslosidrvrcom.h'
      include 'dsdvcom.h'
      include 'dsdmdcom.h' ! for ksoupow
      include 'dsdmdintcom.h' ! for taureschoweq1,locilosi
      integer iwin
      integer nrein
      real*8 rein(nrein),reout
      real*8 radius,dsdmassph
ccc
ccc uncomment if options b) below are activated
ccc      
c      real*8 theta0,mydstausynchsph
ccc
c      real*8 xxmin,xxmax
c      integer nstep
c      common/dstauxxstorecom/xxmin,xxmax,nstep
c      real*8 xx,mydstausynchfunsph
ccc
      if(iwin.eq.ilosigasph.or.iwin.eq.ilosigasph+ilosiabssft) then
ccc
ccc dslosifunsph link for the gamma-ray flux in case of a spherical 
ccc source (or another source just scaling with the dark matter 
ccc emissivity, such as the neutrino flux)      
ccc      
        radius=rein(idvrsph)
        reout=dsdmassph(radius,ksoupow)
        return
      endif
ccc
      if(iwin.eq.ilosisysph.or.iwin.eq.ilosisysph+ilosiabssft) then
ccc
ccc dslosifunsph link for the synchrotron flux in case of a spherical
ccc source        
ccc NOTE: this option is not set in the current driver!!!
ccc      
        write(*,*) 'DS: linking dslosidriver with invalid option'
        write(*,*) 'DS: iwin = ',iwin
        write(*,*) 'DS: program stopped'
        stop
ccc e.g.: dsdmassph * some rescaling function
      endif
ccc
      if(iwin.eq.ilosisysph+ilosinsph) then
ccc
ccc dslosifunsph link for the function needed to compute synchrotron
ccc (self-)absorption in case of a spherical source
ccc NOTE: this option is not set in the current driver!!!
ccc      
        write(*,*) 'DS: linking dslosidriver with invalid option'
        write(*,*) 'DS: iwin = ',iwin
        write(*,*) 'DS: program stopped'
        stop
ccc e.g.: case a) below would be: dsdmassph * some rescaling function
      endif
ccc
      if(iwin.eq.ilosisysph+ilosiabssft+ilositautot) then
ccc
ccc link for the total optical depth for synchrotron flux in case of a   
ccc spherical source; you can implement two options: 
ccc   a) let it compute in dstausph via a line of sight integral of
ccc      dslosifunsph as linked via the integer lochow set here, 
ccc      including eventually a rescaling factor 
ccc   b) compute it here
ccc
ccc e.g.: a) would be:
ccc        
        locilosi=ilosisysph+ilosinsph
        taureschoweq1=1.d0      ! or whatever rescaling is needed
        return
ccc        
ccc b) would be something like:
ccc        
c        locilosi=ilosiabscomp
c        theta0=rein(idvrsph)
c        reout=mydstausynchsph(theta0) ! some function giving the optical
c                                      ! depth in direction theta0
c                                      ! to be included below          
      endif

ccc
      if(iwin.eq.ilosisysph+ilosiabssft+ilositauset) then
ccc
ccc link for the subroutine setting the function computing absoption
ccc along l.o.s.i. for synchrotron flux in case of a spherical source; 
ccc you can implement two options: 
ccc   a) let dstausetsph set  via a line of sight integral of
ccc      dslosifunsph as linked via the integer lochow set here, 
ccc      including eventually a rescaling factor 
ccc   b) save here quantities which are later needed for the computation
ccc
ccc e.g.: a) would be:
ccc        
        locilosi=ilosisysph+ilosinsph ! NOTE: same entry as for
                                      ! iwin=ilosisysph+ilosiabssft+ilositautot
        taureschoweq1=1.d0            ! or whatever rescaling is needed
                                      ! NOTE: also here the same entry as for
                                      ! iwin=ilosisysph+ilosiabssft+ilositautot
        return
ccc        
ccc b) would be:
ccc        
c        locilosi=ilosiabscomp
c        xxmin=rein(idvx)
c        xxmax=rein(idvy)
c        nstep=int(rein(idvz)
      endif

ccc
      if(iwin.eq.ilosisysph+ilosiabssft+ilositaufun) then
ccc
ccc link for the function computing absoption along l.o.s.i. for 
ccc synchrotron flux in case of a spherical source; 
ccc it uses quantities stored above with a call with
ccc     iwin=ilosisysph+ilosiabssft+ilositauset
ccc NOTE: this option is not set in the current driver!!!
ccc      
        if(locilosi.ne.ilosiabscomp) then
          write(*,*) 'DS: wrong sequence in this dslosidriver link'
          write(*,*) 'DS: with iwin = ',iwin
        endif
        write(*,*) 'DS: linking dslosidriver with invalid option'
        write(*,*) 'DS: iwin = ',iwin
        write(*,*) 'DS: program stopped'
        stop        
ccc it would be something like:
ccc        
c        xx=rein(idvrsph)
c        reout=mydstausynchfunsph(xx,xxmin,xxmax,nstep) ! some function
c                                      ! parametrizing absorption
c                                      ! to be included below          
      endif
      write(*,*) 'DS: linking dslosidriver with invalid option'
      write(*,*) 'DS: iwin = ',iwin
      write(*,*) 'DS: program stopped' 
      stop
      end
