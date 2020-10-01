************************************************************************
*** electron/positron equilibrium number density per unit momentum, at 
*** the momentum pp and at the location (R,z), due to wimp pair 
*** annihilations in a static, axisymmetric dark matter halo and up to
*** a normalization factor, see below.
*** set of further assumptions: 
***   - axisymmetric diffusion region with radial boundary condition 
***     neglected (do not use this function with R too close to the 
***     radial size of the diffusion region);
***   - spatial diffusion coefficient and the energy loss rate
***     independent of R and z; 
***   - reacceleration and convection neglected.
*** inputs:
***     R - radial coordinate in kpc
***     z - vertical coordinate in kpc
***     pp - positron momentum in GeV
***     power - integer selecting annihilations (=2) or decays (=1)
***     iv = 1 - the matching between  pp and the variable v is done 
***              assuming fully general momentum loss rate and spatial
***              diffusion coefficient, with a numerical integral + a
***              tabulation involved
***        = 2 - the matching between  pp and the variable v is done 
***              assuming the momentum loss rate scales with pp^2,
***              and the diffusion coefficient has a form which 
***              allows for an analitical inversion between v and pp,
***              see Eq.~(6) & (7) in Baltz & Edsjo (1998)
***     how = 2 - a table for green function is created on first  call, 
***               and then interpolated
***           3 - as 2, but also write the table to disk at the 
***               first call
***           4 - read table from disk on first call, and use the 
***               subsequent calls. If the file does not exist, it 
***               will be created (as in 3). (use as default)
***     labhalo = halo model access label, within the set of models 
***             defined in the halo model repository  
***     loadlab = logical variable: in case it is set to true the halo
***       labhalo is selected within this function, otherwise it is
***       assumed that it has been loaded before linking to this 
***       function such as in the function dsepdphidpaxi
*** output: cm^-3 GeV^-1
***
*** author: Piero Ullio
*** modified: Torsten Bringmann, 12/06/2015 
***           (made replaceable & linked to dscrsource_line)
************************************************************************
      real*8 function dsepdndpaxi(R,z,pp,power,iv,how,labhalo,loadlab)
      implicit none
      include 'dscraxicom.h'
      include 'dsdmdcom.h'  ! to interface dmd halo parameters
      real*8 R,z,pp
      integer power,iv,how
      character(*) labhalo
      logical loadlab
ccc      
      real*8 vvar,pinf,pup,par,pinfloc,puploc,eps,prec,result,
     &  deltaV,dsepgreenaxitab,dseppdotm,dsepdndpaxi_int,
     &  dsepdmasaxi
      real*8 dsmwimp, dscrsource_line            ! function
      real*8 eline,widthline, linecontrib
      integer i, pdg2line, istat
      logical getout
      external dsepdndpaxi_int
ccc
      real*8 Rint,zint,vvint,norm
      integer powerint,dhow,iivv
      common/epdndpaxiintcom/Rint,zint,vvint,norm,powerint,dhow,iivv
      character*12 labhaloloc 
      logical loadlabloc
      common/epdndpaxiint2com/labhaloloc,loadlabloc
      save /epdndpaxiintcom/,/epdndpaxiint2com/
ccc
      if(loadlab) then
         call dsdmdselect_halomodel(labhalo)

ccc
ccc check whether you have called this function for a halo model which
ccc can be used for Milky Way rates or not
ccc      
        if(.not.dmdmw) then
          write(*,*) 'DS: call to dsepdndpaxi with the halo label: ',
     &       labhalo    
          write(*,*) 'DS: which has not been initialized as suitable'
          write(*,*) 'DS: for Milky Way rates. Program stopped'
          stop
        endif
      endif
      loadlabloc=.false.
      labhaloloc=labhalo
ccc
ccc  transfer flags to the common block
ccc
      Rint=R
      zint=z
      powerint=power
      dhow=how
      iivv=iv
ccc
ccc if iv=1 check whether tabulations need to be reloaded
ccc
      if(iv.eq.1) then
        call vvarnumsetup
        call uvarnumsetup
      endif
      vvint=vvar(pp,iv)
ccc
      pinf=pp
      pup=dsmwimp()
      getout=.false.
      par=0.d0
      do i=1,100
        pinfloc=pinf*10.d0**(i-1)
        puploc=pinf*10.d0**i
        if(puploc.gt.pup) then
          puploc=pup
          getout=.true.
        endif
        norm=1.d0
        norm=0.5d0*(dsepdndpaxi_int(pinfloc)+dsepdndpaxi_int(puploc))
        if(norm.lt.1.d-100) goto 200
        eps=1.d-7
        prec=1.d-4
c        write(*,*) i,pinfloc,puploc,norm
        call dsfun_int(dsepdndpaxi_int,pinfloc,puploc,eps,prec,result)

ccc   back with 1.d-30: note that the previous fix missed the 1.d30
ccc in the lines !!! 
        par=par+result*norm*1.d-30
 200    continue
        if(getout) goto 100
      enddo
      write(*,*) 'DS: trouble in dsepdndpaxi, program stopped'
      stop
 100  continue
      dsepdndpaxi=par ! dimensionless
ccc sum the line contribution
c... FIXME: consider taking out the dependence on pspower in dscrsource
      linecontrib=
     &     dscrsource_line(-11,1,power,0.0d0,eline,widthline,pdg2line,
     &        istat)
      if (pdg2line.ne.11) linecontrib=0.0d0 ! TB FIXME: contribution from other
                                            !           lines so far missing
      deltaV=vvint-vvar(dsmwimp(),iv)
      dsepdndpaxi=dsepdndpaxi
c     &     +dssigmav0(11,-11)/dssigmav0tot()
     &     +linecontrib
     &     *dsepgreenaxitab(R,z,DeltaV,power,how,labhalo,loadlabloc) 
      dsepdndpaxi=dsepdndpaxi/dseppdotm(pp,iv)*1.d16 
      dsepdndpaxi=dsepdndpaxi*dsepdmasaxi(R,z,power) ! cm^-3 GeV^-1
      return
      end
ccc
ccc
ccc
      real*8 function dsepdndpaxi_int(pp0)
      implicit none
      real*8 pp0
      real*8 dndp,dsepdndpi,deltav,vvar,dsepgreenaxitab
ccc
      real*8 Rint,zint,vvint,norm
      integer powerint,dhow,iivv
      common/epdndpaxiintcom/Rint,zint,vvint,norm,powerint,dhow,iivv
      character*12 labhaloloc 
      logical loadlabloc
      common/epdndpaxiint2com/labhaloloc,loadlabloc
      save /epdndpaxiint2com/
ccc
c... which do I need this 1.d30 ????      
      dndp=1.d30*dsepdndpi(pp0,powerint) ! GeV^-1
ccc   
ccc check whether the source function is zero:
ccc      
      if(dndp.lt.1.d-30) then
        dsepdndpaxi_int=0.d0
        return
      endif
      deltav=vvint-vvar(pp0,iivv) ! kpc^2
      dsepdndpaxi_int=dndp
     &     *dsepgreenaxitab(Rint,zint,deltav,powerint,dhow,labhaloloc
     &        ,loadlabloc)/norm
      return
      end



