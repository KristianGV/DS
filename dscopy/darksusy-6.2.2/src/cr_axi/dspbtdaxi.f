      real*8 function dspbtdaxi(R,z,tp,powerin,isin,prec,nprec,labhalo
     &   ,loadlab)
************************************************************************
*** antiproton "confinement time" at the position (R,z) 
***    (NOTE: at the moment only z=0 is implemented)
*** for an axisymmetric primary source assumed to extend over the 
*** whole diffusion region, as specified by the function dspbdmasaxi
***
*** output in 10^15 s
***
*** the antiproton flux from wimp pair annihilation is obtained by
*** multiplying the result of this function by the factor:
***   unit_fact * 1/(4 pi) * vp * norm * (sigma v)/(2*mwimp**2) * dN/dtp
*** where:
***   vp = vp(tp) = antiproton velocity for given antiproton kinetic
***      energy tp in unit of c 
***   norm = normalization of dspbdmasaxi = (with standard definition)
***      (rho(R,z))**2 = halo density at (R,z) squared in GeV^2 cm^-6
***   sigma v = dssigmav0tot in cm^3 s^-1
***   mwimp = WIMP mass in GeV
***   dN/dtp = dN/dtp(tp) = number of antiproton within the interval 
***      (tp,tp+dtp) of kinetic energy tp (the latter in GeV) in GeV^-1
***   unit_fact = factor taking into account units =
***    1.d15 s sr^-1 2.99792d10 cm s^-1 GeV^2 cm^-6 cm^3 s^-1 GeV^-3
***      = 2.99792d25 cm^-2 s^-1 GeV^-1 sr^-1
*** 
*** for dark matter particle decays, the antiproton flux is obtained by 
*** multiplying the result of this function by the factor:
***   unit_fact * 1/(4 pi) * vp * norm * (dec-rate)/mwimp * dN/dtp
*** having defined (as opposed to the previous case):
***   norm = normalization of dspbdmasaxi = (with standard definition)
***      rho(R,z) = halo density at (R,z) in GeV cm^-3
***   dec-rate = dsdecratewimp = dark matter particle decay rate in s^-1
***   mwimp = decaying particle mass in GeV
***
*** input:
***   R = radial coordinate for the position at which the flux is
***      measured in kpc (cylindrical coordinate system)
***   z = vertical coordinate for the position at which the flux is
***      measured in kpc (cylindrical coordinate system)
***      NB: only z=0 works at the moment
***   tp = antiproton kinetic energy tp in GeV
***   powerin = integer selecting annihilations (=2) or decays (=1)
***   isin = propagation model option : isin=1 -> model with delta 
***      function approximation for the gas disc; isin=2 -> two-zone
***      model with finite thickness for the gas disk
***   prec = value setting the truncation of the sum in the following 
***      way: the sum oscillates around a mean value; a series of 
***      estimated mean values is computed with a method that is a 
***      slight variant of the method proposed in the PhD Thesis of 
***      Torsten Bringmann, pag. 72, Eq. 5.33. 
***      Educated guess: set, e.g., prec=1.d-3 
***   nprec = number of zeros within which the convergence in the sum
***      is required; if it does not happen, the code stops. This is a
***      flag preventining the code to spend too much time within this
***      routine for propagation or halo models for which the form in 
***      which the analytic solution is not converging. nprec cannot be 
***      lower than 20. educated guess: set, e.g., nprec=80 for prec=1.d-3 
***   labhalo = halo model access label, within the set of models 
***             defined in the halo model repository  
***   loadlab = logical variable: in case it is set to true the halo
***     labhalo is selected within this function, otherwise it is
***     assumed that it has been loaded before linking to this function
***     such as in the function dspbdphidtaxi or dspbtdaxitab      
***
*** other inputs from common blocks:
***   parameters for the propagation model
***   setup for dspbdmasaxi
***
***************************************************************************
      implicit none
      include 'dsmpconst.h'
      include 'dsdmdcom.h'      ! to interface dmd halo parameters
      real*8 R,z,tp,prec
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
          write(*,*) 'DS: call to dspbtdaxi with the halo label: ',
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
        write(*,*) 'DS: call to dspbtdaxi with z = ',z
        write(*,*) 'DS: so far, dspbtdaxi is impemented only for z = 0'
        write(*,*) 'DS: program stopped'
        stop
      endif
      if(nprec.lt.20) then
        write(*,*) 'DS: call to dspbtdaxi with too small nprec = ',nprec
        write(*,*) 'DS: program stopped'
        stop
      endif
ccc
ccc set the varaiable for annihilation or decay:
ccc      
      pbtdaxipower=powerin
ccc
ccc set the antiproton mass and atomic number
ccc
      mnuc=m_p
      anuc=1
      znuc=1
ccc
      dspbtdaxi=dspbsums(R,z,tp,isin,prec,nprec)
      return
      end


      real*8 function dspbsums(R,z,tp,isin,prec,nprec)
************************************************************************
*** sum over zeros of the bessel function J0, plus some checks
*** [1.d15 s]
************************************************************************
      include 'dscraxicom.h'
      real*8 R,z,tp,prec
      integer isin,nprec
      integer s,ii,smax(2)
      integer nn,nnout,iimax,deltas,jj,sload
      real*8 partial,dspbintrs,sum,res(0:100,1000),resave(100),diff
      real*8 ck1,ck2,ckcut,dspbdmasaxi
ccc
      integer pbtdaxipower
      common/pbtdaxipowercom/pbtdaxipower
ccc
      real*8 diffstore
      common/pbcintrloccom/diffstore
ccc
      real*8 Rcint,zcint
      common/pbcintrcom/Rcint,zcint
ccc
      integer nrstep,nzstep
      common/pbintstepcom/nrstep,nzstep
ccc
      real*8 pbnpairsnorm
      common/pbintnormcom/pbnpairsnorm
ccc
ccc fix the normalization of dspbdmasaxi
ccc
      pbnpairsnorm=dspbdmasaxi(R,z,pbtdaxipower)
ccc
ccc decide whether to do the integrations in 1 shot or to split it
ccc      
      ck1=dspbdmasaxi(max(Rcint,0.d0),0.d0,pbtdaxipower)      
      ck2=dspbdmasaxi(diffRh,0.d0,pbtdaxipower) 
      ckcut=6.d0  ! what is the best value here ???
      if(dlog10(ck1/ck2).gt.ckcut) then
        nrstep=2
        nzstep=2
      else
        nrstep=1
        nzstep=1
      endif
ccc
      s=0
      diffstore=1.d10
      smax(1)=0
      smax(2)=0
      nmax=0
      ii=0
      sum=0.d0
      diff=0.d0
 100  s=s+1
      partial=dspbintrs(R,z,tp,s,isin)
      sum=sum+partial
      res(0,s)=sum
c      write(*,*) s,res(0,s)
ccc
ccc  improving conversion implementing averages of cycles of oscillations, 
ccc  with a method that is a slight variant of the method proposed in the 
ccc  PhD Thesis of Torsten Bringmann, pag. 72, Eq. 5.33;
ccc  try to find the oscillation length between 4th and 5th maximum:
ccc
      if(s.gt.2.and.res(0,s-1).gt.res(0,s).and.res(0,s-1).gt.res(0,s-2)
     &   .and.smax(2).eq.0) then
        if(nmax.lt.4) then
          nmax=nmax+1
        endif
        if(smax(1).eq.0.and.nmax.ge.4) then ! fourth max
          smax(1)=s-1
        elseif(smax(1).gt.0.and.nmax.ge.4) then ! fith max
          smax(2)=s-1
          deltas=smax(2)-smax(1)
c          write(*,*) 'deltas = ',deltas
ccc
ccc fill in the first table of the sums
ccc
          nnmax=int(dble(s)/dble(deltas-1))
          do nn=1,nnmax
            iimax=s-(deltas-1)*nn
            do ii=1,iimax
              partial=0.d0
              do jj=ii,ii+deltas-1
                partial=partial+res(nn-1,jj)
              enddo
              res(nn,ii)=partial/dble(deltas)
            enddo
          enddo
          sload=s
        endif
      endif
      if(s.le.sload+deltas.and.smax(2).ne.0) then
ccc
ccc fill in more the table with one extra nn
ccc
        nnmax=nnmax+1
        do nn=1,nnmax
          ii=s-(deltas-1)*nn
          partial=0.d0
          do jj=ii,ii+deltas-1
            partial=partial+res(nn-1,jj)
          enddo
          res(nn,ii)=partial/dble(deltas)
          sload=s
        enddo
ccc
ccc check whether there is convergence        
ccc       
        if(nnmax.ge.3) then
          do nn=1,nnmax
            iimax=s-(deltas-1)*nn
            if(iimax.ge.deltas) then
              ii=iimax
              resave(nn)=0.d0
              do ii=iimax-deltas+1,iimax
                resave(nn)=resave(nn)+res(nn,ii)
              enddo
              resave(nn)=resave(nn)/dble(deltas)
              if(nn.gt.3) then
                diff=dabs(resave(nn)-resave(nn-1))/dabs(resave(nn))
                if(diff.lt.prec) then
                  ii=iimax
c                  write(*,*) nnmax,nn,iimax,resave(nn)
c                  do ii=iimax-deltas+1,iimax
c                    write(*,*) ii,res(nn,ii)
c                  enddo
                  nnout=nn
                  goto 200
                endif
              endif
            endif
          enddo
        endif
      endif
      if(s.ge.nprec.and.diff.gt.prec) then
        write(*,*) 'DS: in dspbsums the tool for faster convergence'
        write(*,*) 'DS: is not working; you need to debug this case'
        write(*,*) 'DS: R,z,tp,isin,nprec : ',R,z,tp,isin,nprec
        write(*,*) 'DS: diff / prec : ',diff,prec
        write(*,*) 'DS: program stopped'
        stop
      endif
      goto 100
 200  continue
      dspbsums=resave(nnout)
      return
      end


      real*8 function dspbintrs(R,z,tp,s,isin)
************************************************************************
*** R integration [1.d15 s]
************************************************************************
      implicit none
      include 'dscraxicom.h'
ccc
      real*8 R,z,tp
      integer s,isin
ccc
      integer si
      real*8 Ri,tpi
      common/pbintrcom/Ri,tpi,si
ccc
      integer is
      common/pboptioncom/is
ccc
      real*8 rcint,zcint
      common/pbcintrcom/rcint,zcint
ccc
      integer nrstep,nzstep
      common/pbintstepcom/nrstep,nzstep
ccc
      real*8 rmin,rmax,prec,eps,res,par,dspbint_intgR2
      external dspbint_intgR2

ccc
      real*8 tollf2,tollab
      integer nmax,ihow
ccc
      Ri=R
      tpi=tp
      si=s
      is=isin
ccc
      if(nrstep.eq.1) then
        rmin=0.d0
        rmax=diffRh
        eps=1.d-10
        prec=1.d-3
        call dsfun_intb(dspbint_intgR2,rmin,rmax,eps,prec,res)
        dspbintrs=res                          ! 1.d15 s
      else
        rmin=0.d0
        rmax=diffrcpb
        if(rmax.lt.1.d-10) rmax=1.d-10
        if(rmax.lt.rcint) rmax=rcint
        prec=1.d-3
        eps=min(1.d-10,dabs(dspbint_intgR2(Rmax))*Rmax**2*prec*1.d-2)
        call dsfun_intb(dspbint_intgR2,rmin,rmax,eps,prec,res)
        par=res
        rmin=rmax
        rmax=diffRh
c
        tollf2=1.d5
        tollab=(dlog(rmax)-dlog(rmin))/1000.d0
        nmax=1000
        prec=1.d-3
        ihow=2
        call dsfun_intparb(dspbint_intgR2,dspbint_intgR2,rmin,rmax,
     &      tollf2,tollab,nmax,prec,ihow,res)
        par=par+res                          ! 1.d15 s
        dspbintrs=par
      endif
      return
      end



      real*8 function dspbint_intgR2(R0)
************************************************************************
*** integrand for R integration [1.d15 s kpc^-1]
************************************************************************
      implicit none
      include 'dscraxicom.h'
      real*8 R0,dspbint_intgRis1,dspbint_intgRis2
ccc
      integer si
      real*8 Ri,tpi
      common/pbintrcom/Ri,tpi,si
ccc
      integer is
      common/pboptioncom/is
ccc
      dspbint_intgr2=0.d0
      if(is.eq.1) then
        dspbint_intgR2=dspbint_intgRis1(si,tpi,R0,Ri)
      elseif(is.eq.2) then
        dspbint_intgR2=dspbint_intgRis2(si,tpi,R0,Ri)
      endif
      dspbint_intgR2=2.d0*R0*dspbint_intgR2/diffRh**2 
                                           ! 1.d15 s kpc^-1
      return
      end
ccc
ccc
      real*8 function dspbint_intgRis1(s,tp,R0,R)
************************************************************************
*** integrand for R integration [1.d15 s]
************************************************************************
      implicit none
      include 'dscraxicom.h'
      include 'dsmpconst.h'
      integer s
      real*8 dbesj0
      real*8 tp,R0,R
      real*8 diff,sum,pp,ee,rig,dskdiff,beta,axsec,dspbsigmavpbar,
     &  dsdbsigmavdbar
      real*8 dspbbesselzeroj0,dspbbesselfactorj0,dspbint_intz
ccc
      real*8 tpstore,nusk,factorj
      integer anucstore,sstore,isstore
      common/pbintrlocis12com/tpstore,nusk,factorj,anucstore,sstore,
     & isstore
ccc
      real*8 dpbarh,ratiodpbar,vbar,pebarg,pebarh,lambdag,lambdah
     & ,ratiolambda
      common/pbintglobcom/dpbarh,ratiodpbar,vbar,pebarg,pebarh,lambdag
     & ,lambdah,ratiolambda 
ccc
      real*8 mnuc
      integer anuc,znuc
      common/pbnucleoncom/mnuc,anuc,znuc
ccc
      diff=dabs(tp-tpstore)
      sum=dabs(tp+tpstore)
      if((diff.gt.1.d-8*sum).or.(anuc.ne.anucstore).or.(isstore.ne.1)) 
     &  then
ccc
ccc set tp dependent variables:
ccc
        ee=anuc*tp+mnuc
        pp=dsqrt(dabs(ee**2-mnuc**2))
        rig=pp/dble(znuc)
        beta=pp/ee
        dpbarh=dskdiff(rig,1,beta)
        vbar=diffcvel/dpbarh   ! 10^5 cm/s * 10^-27 s/cm^-2 = 
                               ! 10^-22 cm^-1 
        vbar=vbar*(kpc/10.d0)  ! kpc^-1
        if(anuc.eq.1) then
          axsec=dspbsigmavpbar(ee)  ! mb * 1.d10 cm*s^-1 
        elseif(anuc.eq.2) then
          axsec=dsdbsigmavdbar(ee)  ! mb * 1.d10 cm*s^-1
        else
          write(*,*) 'DS: call to dspbint_intgRis1 without proper'
          write(*,*) 'DS: intialization of mnuc, anuc : ',mnuc,anuc
          write(*,*) 'DS: program stopped'
          stop
        endif
        pebarg=2.d0*diffhg*diffng*axsec/dpbarh  
            ! kpc * cm^-3 * 1.d-27 cm^2 * 1.d10 cm*s^-1 * 10^-27 s/cm^-2
            ! kpc* 10^-44 * cm^-2
        pebarg=pebarg*(kpc/10.d0)**2 ! kpc * kpc^-2 = kpc^-1
ccc
ccc you neeed to reset s dependent variables as well
ccc
        sstore=-1
        tpstore=tp
        isstore=1
      endif
      if(s.ne.sstore) then
ccc
ccc set s dependent variables:
ccc
        nusk=dspbbesselzeroj0(s)
        factorj=dspbbesselfactorj0(s,R)
        if(vbar.gt.1.d-16) then
          lambdah=dsqrt((nusk/diffRh)**2+vbar**2/4.d0) ! kpc^-1
        else
          lambdah=(nusk/diffRh) ! kpc^-1
        endif
        sstore=s
      endif
      dspbint_intgRis1=dspbint_intz(R0)     ! kpc^2
      dspbint_intgRis1=dspbint_intgRis1*dbesj0(nusk*R0/diffRh)
     &                 *factorj/dpbarh  ! kpc^2 * 10^-27 s/cm^-2
      dspbint_intgRis1=dspbint_intgRis1*kpc**2 ! 1.d15 s
      return
      end
ccc
ccc
      real*8 function dspbint_intgRis2(s,tp,R0,R)
************************************************************************
*** integrand for R integration [1.d15 s]
************************************************************************
      implicit none
      include 'dscraxicom.h'
      include 'dsmpconst.h'
      real*8 dbesj0
      integer s
      real*8 tp,R0,R
      real*8 diff,sum,pp,rig,ee,dpbarg,dskdiff,beta,axsec,
     &  dspbsigmavpbar,dsdbsigmavdbar
      real*8 dspbbesselzeroj0,dspbbesselfactorj0,dspbinte_intz
ccc
      real*8 tpstore,nusk,factorj
      integer anucstore,sstore,isstore
      common/pbintrlocis12com/tpstore,nusk,factorj,anucstore,sstore,
     & isstore
ccc
      real*8 dpbarh,ratiodpbar,vbar,pebarg,pebarh,lambdag,lambdah
     & ,ratiolambda
      common/pbintglobcom/dpbarh,ratiodpbar,vbar,pebarg,pebarh,lambdag
     & ,lambdah,ratiolambda 
ccc
      real*8 mnuc
      integer anuc,znuc
      common/pbnucleoncom/mnuc,anuc,znuc
ccc
      diff=dabs(tp-tpstore)
      sum=dabs(tp+tpstore)
      if((diff.gt.1.d-8*sum).or.(anuc.ne.anucstore).or.(isstore.ne.2)) 
     &  then
ccc
ccc set tp dependent variables:
ccc
        ee=anuc*tp+mnuc
        pp=dsqrt(dabs(ee**2-mnuc**2))
        rig=pp/dble(znuc)
        beta=pp/ee
        dpbarh=dskdiff(rig,1,beta)
        vbar=diffcvel/dpbarh   ! 10^5 cm/s * 10^-27 s/cm^-2 
                               ! = 10^-22 cm^-1 
        vbar=vbar*(kpc/10.d0)  ! kpc^-1
        if(anuc.eq.1) then
          axsec=dspbsigmavpbar(ee)  ! mb * 1.d10 cm*s^-1 
        elseif(anuc.eq.2) then
          axsec=dsdbsigmavdbar(ee)  ! mb * 1.d10 cm*s^-1
        else
          write(*,*) 'DS: call to dspbint_intgRis2 without proper'
          write(*,*) 'DS: intialization of mnuc, anuc : ',mnuc,anuc
          write(*,*) 'DS: program stopped'
          stop
        endif
        pebarh=diffnh*axsec/dpbarh  
            ! cm^-3 * 1.d-27 cm^2 * 1.d10 cm*s^-1 * 10^-27 s/cm^-2
            ! 10^-44 * cm^-2
        pebarh=pebarh*(kpc/10.d0)**2 ! kpc^-2 = kpc^-2
        dpbarg=dskdiff(rig,2,beta)
        pebarg=diffng*axsec/dpbarg  
            ! cm^-3 * 1.d-27 cm^2 * 1.d10 cm*s^-1 * 10^-27 s/cm^-2
            ! 10^-44 * cm^-2
        pebarg=pebarg*(kpc/10.d0)**2 ! kpc^-2 = kpc^-2
        ratiodpbar=dpbarg/dpbarh
ccc
ccc you neeed to reset s dependent variables as well
ccc
        sstore=-1
        tpstore=tp
        isstore=2
      endif
      if(s.ne.sstore) then
ccc
ccc set s dependent variables:
ccc
        nusk=dspbbesselzeroj0(s)
        factorj=dspbbesselfactorj0(s,R)
        if(vbar.gt.1.d-16) then
          lambdah=dsqrt((nusk/diffRh)**2+vbar**2/4.d0) ! kpc^-1
        else
          lambdah=(nusk/diffRh) ! kpc^-1
        endif
        if(pebarg.gt.1.d-16) then
          lambdag=dsqrt((nusk/diffRh)**2+pebarg) ! kpc^-1
        else
          lambdag=(nusk/diffRh) ! kpc^-1
        endif
        ratiolambda=lambdah/lambdag
        sstore=s
      endif
      dspbint_intgRis2=dspbinte_intz(R0)    ! kpc^2
      dspbint_intgRis2=dspbint_intgRis2*dbesj0(nusk*R0/diffRh)
     &                 *factorj/dpbarh ! kpc^2 * 10^-27 s/cm^-2
      dspbint_intgRis2=dspbint_intgRis2*kpc**2 ! 1.d15 s
      return
      end



      real*8 function dspbint_intz(R0)
************************************************************************
*** z integration, is=1 case [kpc^2]
************************************************************************
      implicit none
      include 'dscraxicom.h'
      real*8 R0
ccc
      real*8 rcint,zcint
      common/pbcintrcom/rcint,zcint
ccc
      integer nrstep,nzstep
      common/pbintstepcom/nrstep,nzstep
ccc
      real*8 dpbarh,ratiodpbar,vbar,pebarg,pebarh,lambdag,lambdah
     & ,ratiolambda
      common/pbintglobcom/dpbarh,ratiodpbar,vbar,pebarg,pebarh,lambdag
     & ,lambdah,ratiolambda
ccc
      real*8 R0i
      common/pbintzloccom/R0i
ccc
      real*8 zmin,zmax,tollf2,tollab,prec,res,partial,eps,dspbint_intgz
      integer nmax,ihow
      external dspbint_intgz
ccc
      R0i=R0
      if(nzstep.eq.1) then
        zmin=0.d0
ccc
ccc eventually skip the inner cylinder
ccc
        if(R0.lt.rcint.and.zmin*diffhh.lt.zcint) zmin=zcint/diffhh
        zmax=1.d0
        eps=1.d-10
        prec=1.d-3
        call dsfun_int(dspbint_intgz,zmin,zmax,eps,prec,res)
        partial=res*diffhh             ! kpc 
      else
        zmin=0.d0
        zmax=(max(diffrcpb,1.d-10))/diffhh
ccc
ccc eventually skip the inner cylinder
ccc
        if(R0.lt.rcint.and.zmax*diffhh.lt.zcint) then
          zmax=zcint/diffhh
          partial=0.d0
          goto 333
        endif
        prec=1.d-3
        eps=dabs(dspbint_intgz(zmax))*zmax*prec*1.d-2
        call dsfun_int(dspbint_intgz,zmin,zmax,eps,prec,res)
        partial=res
c        write(*,*) s,zmin,zmax,res
 333    zmin=zmax
        zmax=1.d0
        tollf2=1.d5
        tollab=(dlog(zmax)-dlog(zmin))/1000.d0
        nmax=1000
        prec=1.d-3
        ihow=2
        call dsfun_intpar(dspbint_intgz,dspbint_intgz,zmin,zmax,
     &      tollf2,tollab,nmax,prec,ihow,res)
        partial=partial+res              ! dimensionless
        partial=partial*diffhh             ! kpc
      endif
      partial=partial
     &          /((pebarg+vbar)/2.d0+lambdah/dtanh(lambdah*diffhh))
      dspbint_intz=partial            ! kpc^2
      return
      end
ccc
ccc
ccc
      real*8 function dspbint_intgz(z0i)
************************************************************************
*** integrand for z integration, is=1 case dimensionless
************************************************************************
      implicit none
      include 'dscraxicom.h'
      real*8 z0i
      real*8 partial,arg,dspbdmasaxi
ccc
      real*8 dpbarh,ratiodpbar,vbar,pebarg,pebarh,lambdag,lambdah
     & ,ratiolambda
      common/pbintglobcom/dpbarh,ratiodpbar,vbar,pebarg,pebarh,lambdag
     & ,lambdah,ratiolambda
ccc
      real*8 R0i
      common/pbintzloccom/R0i
ccc
      real*8 pbnpairsnorm
      common/pbintnormcom/pbnpairsnorm
ccc
      integer pbtdaxipower
      common/pbtdaxipowercom/pbtdaxipower
ccc
      arg=lambdah*diffhh*(1.d0-dabs(z0i))
      if(arg.lt.37.d0) then
        if(dabs(z0i).gt.0.1d0) then
          partial=dsinh(arg)/dsinh(lambdah*diffhh)
        else
          partial=dcosh(lambdah*diffhh*dabs(z0i))
     &            -dsinh(lambdah*diffhh*dabs(z0i))/dtanh(lambdah*diffhh)
        endif
      else
        arg=-lambdah*diffhh*dabs(z0i)
        partial=dexp(arg)
      endif
      dspbint_intgz=partial*dspbdmasaxi(R0i,diffhh*z0i,pbtdaxipower)
     &  /pbnpairsnorm
      return
      end


      real*8 function dspbinte_intz(R0)
************************************************************************
*** z integration, is=2 case [kpc^2]
************************************************************************
      implicit none
      include 'dscraxicom.h'
      real*8 R0
ccc
      real*8 rcint,zcint
      common/pbcintrcom/rcint,zcint
ccc
      real*8 dpbarh,ratiodpbar,vbar,pebarg,pebarh,lambdag,lambdah
     & ,ratiolambda
      common/pbintglobcom/dpbarh,ratiodpbar,vbar,pebarg,pebarh,lambdag
     & ,lambdah,ratiolambda
ccc
      real*8 R0i
      common/pbintzloccom/R0i
ccc
      integer nrstep,nzstep
      common/pbintstepcom/nrstep,nzstep
ccc
      real*8 zmin,zmax,tollf2,tollab,prec,eps,res,partial,dspbinte_intgz
      integer nmax,ihow
      external dspbinte_intgz
ccc
      R0i=R0
      if(nzstep.eq.1) then
        zmin=0.d0
ccc
ccc eventually skip the inner cylinder
ccc
        if(R0.lt.rcint.and.zmin*diffhh.lt.zcint) zmin=zcint/diffhh
        zmax=1.d0
        eps=1.d-10
        prec=1.d-3
        call dsfun_int(dspbinte_intgz,zmin,zmax,eps,prec,res)
        partial=res*diffhh             ! kpc
      else
        zmin=0.d0
        zmax=(max(diffrcpb,1.d-10))/diffhh
ccc
ccc eventually skip the inner cylinder
ccc
        if(R0.lt.rcint.and.zmax*diffhh.lt.zcint) then
          zmax=zcint/diffhh
          partial=0.d0
          goto 333
        endif
        prec=1.d-3
        eps=dabs(dspbinte_intgz(zmax))*zmax*prec*1.d-2
        call dsfun_int(dspbinte_intgz,zmin,zmax,eps,prec,res)
        partial=res
c        write(*,*) s,zmin,zmax,res
 333    zmin=zmax
        zmax=1.d0
        tollf2=1.d5
        tollab=(dlog(zmax)-dlog(zmin))/1000.d0
        nmax=1000
        prec=1.d-3
        ihow=2
        call dsfun_intpar(dspbinte_intgz,dspbinte_intgz,zmin,zmax,
     &      tollf2,tollab,nmax,prec,ihow,res)
        partial=partial+res             ! dimensionless
        partial=partial*diffhh             ! kpc
      endif
      partial=partial/(vbar/2.d0+lambdah/dtanh(lambdah*(diffhh-diffhg))
     &                 +ratiodpbar*lambdag*dtanh(lambdag*diffhg))
      dspbinte_intz=partial            ! kpc^2
      return
      end
ccc
ccc
ccc
      real*8 function dspbinte_intgz(z0i)
************************************************************************
*** integrand for z integration, is=2 case [dimensionless]
************************************************************************
      implicit none
      include 'dscraxicom.h'
      real*8 z0i
      real*8 partial,arg,dspbdmasaxi
ccc
      real*8 dpbarh,ratiodpbar,vbar,pebarg,pebarh,lambdag,lambdah
     & ,ratiolambda
      common/pbintglobcom/dpbarh,ratiodpbar,vbar,pebarg,pebarh,lambdag
     & ,lambdah,ratiolambda
ccc
      real*8 R0i
      common/pbintzloccom/R0i
ccc
      real*8 pbnpairsnorm
      common/pbintnormcom/pbnpairsnorm
ccc
      integer pbtdaxipower
      common/pbtdaxipowercom/pbtdaxipower
ccc
      if(dabs(z0i)*diffhh.ge.diffhg) then
        arg=lambdah*diffhh*(1.d0-dabs(z0i))
        if(arg.lt.37.d0) then
          partial=dsinh(arg)
     &             /dsinh(lambdah*(diffhh-diffhg))/dcosh(lambdag*diffhg)
        else
          arg=-lambdah*diffhh*dabs(z0i)+diffhg*(lambdah-lambdag)
          partial=2.d0*dexp(arg)
        endif
      else
        arg=lambdag*(diffhg-diffhh*dabs(z0i))
        if(arg.lt.37.d0) then
          partial=(dcosh(arg)/dcosh(lambdag*diffhg)+
     &                     dsinh(arg)/dcosh(lambdag*diffhg)/ratiodpbar*
     &     (vbar/2.d0/lambdag
     &      +ratiolambda/dtanh(lambdah*(diffhh-diffhg))))
        else
          arg=-lambdag*diffhh*dabs(z0i)
          partial=dexp(arg)*(1.d0+1.d0/ratiodpbar*
     &     (vbar/2.d0/lambdag
     &      +ratiolambda/dtanh(lambdah*(diffhh-diffhg))))
        endif
      endif
      if(dabs(z0i)*diffhh.ge.diffhg) 
     &  partial=partial*dexp(-vbar*(diffhh*dabs(z0i)-diffhg)/2.d0)
      dspbinte_intgz=partial*dspbdmasaxi(R0i,diffhh*z0i,pbtdaxipower)
     &  /pbnpairsnorm
      return
      end
