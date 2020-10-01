      real*8 function dsdbtdaxi(R,z,td,powerin,isin,prec,nprec,labhalo
     &   ,loadlab)
************************************************************************
*** antideuteron "confinement time" at the position (R,z) 
***    (NOTE: at the moment only z=0 is implemented)
*** for an axisymmetric primary source assumed to extend over the 
*** whole diffusion region, as specified by the function dspbdmasaxi
***
*** input:
***     td - kinetic energy per nucleon (gev)
***
*** see the header of the function dspbtdaxi for details on other inputs
*** and outputs, since they are in perfect analogy to those here. 
*** NOTE: the setup with the inner cylinder at the Galactic center 
*** excluded from the spatial integration and treated with the point 
*** source green function to speed up convergence is NOT changed here, 
*** assuming that convergence issues, for a given source function and R 
*** are analogous
***
*** output in 10^15 s
***      
************************************************************************
      implicit none
      include 'dsmpconst.h'
      include 'dsdmdcom.h'  ! to interface dmd halo parameters
      real*8 R,z,td,prec
      integer powerin,isin,nprec
      character(*) labhalo
      logical loadlab
ccc      
      real*8 dspbsums
ccc
      integer pbtdaxipower
      common/pbtdaxipowercom/pbtdaxipower
ccc
      real*8 mnuc
      integer anuc,znuc
      common/pbnucleoncom/mnuc,anuc,znuc
ccc
      integer ibessel
      common/ibesselinicom/ibessel
ccc
      if(loadlab) then
        call dsdmdselect_halomodel(labhalo)      
ccc
ccc check whether you have called this function for a halo model which
ccc can be used for Milky Way rates or not
ccc      
        if(.not.dmdmw) then
          write(*,*) 'DS: call to dsdbtdaxi with the halo label: ',
     &       labhalo    
          write(*,*) 'DS: which has not been initialized as suitable'
          write(*,*) 'DS: for Milky Way rates. Program stopped'
          stop
        endif
      endif  
ccc
      if(ibessel.ne.1) then
c        write(*,*) 'initializing bessel stuff' 
        call dspbbesselini
        ibessel=1
      endif
ccc
      if(dabs(z).gt.1.d-15) then
        write(*,*) 'DS: call to dsdbtdaxi with z = ',z
        write(*,*) 'DS: so far, dsdbtdaxi is impemented only for z = 0'
        write(*,*) 'DS: program stopped'
        stop
      endif
      if(nprec.lt.20) then
        write(*,*) 'DS: call to dsdbtdaxi with too small nprec = ',nprec
        write(*,*) 'DS: program stopped'
        stop
      endif
ccc
ccc set the varaiable for annihilation or decay:
ccc      
      pbtdaxipower=powerin
ccc
ccc set the antideuteron mass and atomic number
ccc
      mnuc=m_d
      anuc=2
      znuc=1
ccc
      dsdbtdaxi=dspbsums(R,z,td,isin,prec,nprec)
      return
      end

