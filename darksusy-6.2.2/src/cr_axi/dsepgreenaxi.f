      real*8 function dsepgreenaxi(R,z,DeltaV,powerin,labhalo,loadlab)
************************************************************************
*** green function to derive the electron/positron equilibrium number 
*** density due to wimp pair annihilations in a static, axisymmetric 
*** dark matter halo. see dsepdndpaxi for a complete list of underlying
*** assumptions.
*** inputs:
***     R - radial position of the observer in kpc
***     z - vertical position of the observer in kpc
***     DeltaV - increment in the variable vvar in kpc^2       
***     4*DeltaVin = lambda**2, with lambda the diffusion length in kpc 
***     labhalo = halo model access label, within the set of models 
***             defined in the halo model repository  
***     loadlab = logical variable: in case it is set to true the halo
***       labhalo is selected within this function, otherwise it is
***       assumed that it has been loaded before linking to this 
***       function such as in the function dsepdphidpaxi, dsepdndpaxi or
***       dsepgreenaxitab    
***
*** output: dimensionless
************************************************************************
      implicit none
      include 'dscraxicom.h'
      include 'dsmpconst.h' 
      include 'dsdmdcom.h'  ! to interface dmd halo parameters
      real*8 R,z,DeltaV
      integer powerin
      character(*) labhalo
      logical loadlab
ccc
      real*8 rmin,rmax,eps,prec,res,dsepgreenzaxi,par,tollf2,tollab,ck1,
     &  ck2,dsepdmasaxi,ckcut,rstore,resguess,rck
      integer nmax,ihow
      external dsepgreenzaxi
ccc
      integer epgraxipower
      common/epgraxipowercom/epgraxipower
ccc
      real*8 dv,Ltilde,Rtilde,zhat,R0tildei,npairs_norm
      common/epintegrcom/dv,Ltilde,Rtilde,zhat,R0tildei,npairs_norm
ccc
      integer nrstep,nzstep
      common/epintstepcom/nrstep,nzstep
ccc
      if(DeltaV.lt.1.d-16) then
        write(*,*) 'DS: calling dsepgreenaxi with invalid DeltaV = ',
     &    DeltaV
        dsepgreenaxi=0.d0
        return
      endif
ccc
      if(loadlab) then
        call dsdmdselect_halomodel(labhalo)
ccc
ccc check whether you have called this function for a halo model which
ccc can be used for Milky Way rates or not
ccc      
        if(.not.dmdmw) then
          write(*,*) 'DS: call to dsepgreenaxi with the halo label: ',
     &       labhalo    
          write(*,*) 'DS: which has not been initialized as suitable'
          write(*,*) 'DS: for Milky Way rates. Program stopped'
          stop
        endif  
      endif
ccc
ccc set the varaiable for annihilation or decay:
ccc      
      epgraxipower=powerin
ccc
      dv=DeltaV
      Ltilde=diffhh/dsqrt(dv)
      Rtilde=R/dsqrt(dv)
      zhat=z/diffhh
      npairs_norm=dsepdmasaxi(R,z,epgraxipower)
ccc
ccc decide whether to do the integrations in 1 shot or to split it
ccc
      rck=diffrcep      
      ck1=dsepdmasaxi(rck,0.d0,epgraxipower)
     &     *dexp(-(R-rck)**2/dv/4.d0)*rck**2 
c      write(*,*) rck,ck1
      rck=max(R*0.9d0,diffrcep)
      ck2=dsepdmasaxi(rck,0.d0,epgraxipower)
     &     *dexp(-(R-rck)**2/dv/4.d0)*rck**2 
      ckcut=5.d0  ! is this the best value ???
c      write(*,*) rck,ck2,dlog10(ck1/ck2)
      if(dlog10(ck1/ck2).gt.ckcut) then
        nrstep=2
        nzstep=2
      else
        nrstep=1
        nzstep=1
      endif
ccc
      if(nrstep.eq.1) then
        rmin=0.d0
        rmax=Rtilde
        eps=1.d-10
        prec=1.d-3
        call dsfun_intb(dsepgreenzaxi,rmin,rmax,eps,prec,res)
        par=res
c        write(*,*) 'dv,rmin,rmax,res : ',dv,rmin,rmax,res
        rmin=rmax
        rmax=diffRh/dsqrt(dv)
c        write(*,*) diffRh,dsqrt(dv),rmax
        eps=1.d-10
        prec=1.d-3
        call dsfun_intb(dsepgreenzaxi,rmin,rmax,eps,prec,res)
c        write(*,*) 'dv,rmin,rmax,res : ',dv,rmin,rmax,res
        dsepgreenaxi=par+res    
      else
        rmin=0.d0
        rmax=diffrcep
        if(rmax.lt.1.d-10) rmax=1.d-10
        rmax=rmax/dsqrt(dv)
        rstore=rmax
ccc
ccc add it later only if you estimate that it is needed
ccc
        rmin=rmax
        rmax=Rtilde
        tollf2=1.d5
        tollab=(dlog(rmax)-dlog(rmin))/1000.d0
        nmax=1000
        prec=1.d-3
        ihow=2
        call dsfun_intparb(dsepgreenzaxi,dsepgreenzaxi,rmin,rmax,
     &      tollf2,tollab,nmax,prec,ihow,res)
        par=res
c        write(*,*) '2: ',par,res
        rmin=rmax
        rmax=diffRh/dsqrt(dv)
c        write(*,*) diffRh,dsqrt(dv),rmax
        eps=1.d-10
        prec=1.d-3
        call dsfun_intb(dsepgreenzaxi,rmin,rmax,eps,prec,res)
        par=par+res
c        write(*,*) '3: ',par,res
        rmin=0.d0
        rmax=rstore
        resguess=dsepgreenzaxi(Rmax)*Rmax/2.d0
        if(resguess.gt.par*prec*1.d-3) then
          prec=1.d-3
          eps=resguess*prec*1.d-2
          call dsfun_intb(dsepgreenzaxi,rmin,rmax,eps,prec,res)
          par=par+res
c          write(*,*) '1: ',par,res,resguess
        endif
        dsepgreenaxi=par
      endif
      dsepgreenaxi=dsepgreenaxi/dsqrt(pi)/4.d0
      return
      end
ccc
ccc
ccc
      real*8 function dsepgreenzaxi(R0tilde)
      implicit none
      include 'dscraxicom.h'
      real*8 R0tilde
ccc
      real*8 zhatmin,zhatmax,eps,prec,res,dsepgreenzaxi2,par,tollf2,
     &  tollab,dsbessei0
      integer nmax,ihow
      external dsepgreenzaxi2
ccc
      real*8 dv,Ltilde,Rtilde,zhat,R0tildei,npairs_norm
      common/epintegrcom/dv,Ltilde,Rtilde,zhat,R0tildei,npairs_norm
ccc
      integer nrstep,nzstep
      common/epintstepcom/nrstep,nzstep
ccc
      R0tildei=R0tilde
ccc
      if(dabs(zhat)*diffhh.lt.1.d-15) then
ccc
      if(nzstep.eq.1) then
        zhatmin=0.d0
        zhatmax=1.d0
        eps=1.d-10
        prec=1.d-3
        call dsfun_int(dsepgreenzaxi2,zhatmin,zhatmax,eps,prec,res)
        dsepgreenzaxi=2.d0*res   
      else
        zhatmin=0.d0
        zhatmax=diffrcep
        if(zhatmax.lt.1.d-10) zhatmax=1.d-10
        zhatmax=zhatmax/diffhh
        prec=1.d-3
        eps=min(1.d-10,dabs(dsepgreenzaxi2(zhatmax))*zhatmax*prec*1.d-2)
        call dsfun_int(dsepgreenzaxi2,zhatmin,zhatmax,eps,prec,res)
        par=res
        zhatmin=zhatmax
        zhatmax=1.d0
        tollf2=1.d5
        tollab=(dlog(zhatmax)-dlog(zhatmin))/1000.d0
        nmax=1000
        prec=1.d-3
        ihow=2
        call dsfun_intpar(dsepgreenzaxi2,dsepgreenzaxi2,zhatmin,zhatmax,
     &      tollf2,tollab,nmax,prec,ihow,res)
        par=par+res                     
        dsepgreenzaxi=2.d0*par
      endif
      else
        zhatmin=-1.d0
        zhatmax=zhat
        eps=1.d-10
        prec=1.d-3
        call dsfun_int(dsepgreenzaxi2,zhatmin,zhatmax,eps,prec,res)
        par=res
        zhatmin=zhatmax
        zhatmax=1.d0
        eps=1.d-10
        prec=1.d-3
        call dsfun_int(dsepgreenzaxi2,zhatmin,zhatmax,eps,prec,res)
        dsepgreenzaxi=par+res
      endif
c      write(*,*) R0tilde,dsepgreenzaxi
      dsepgreenzaxi=dsepgreenzaxi*R0tilde
     &       *dsbessei0(R0tilde*Rtilde/2.d0)
     &       *dexp(-(R0tilde-Rtilde)**2/4.d0)
c      write(*,*) R0tilde,dsepgreenzaxi,dsbessei0(R0tilde*Rtilde/2.d0)
      return
      end
ccc      
ccc
ccc
      real*8 function dsepgreenzaxi2(z0hat)
      implicit none
      include 'dscraxicom.h'
ccc
      integer epgraxipower
      common/epgraxipowercom/epgraxipower
ccc
      real*8 dv,Ltilde,Rtilde,zhat,R0tildei,npairs_norm
      common/epintegrcom/dv,Ltilde,Rtilde,zhat,R0tildei,npairs_norm
ccc
      real*8 z0hat
      real*8 sign,znp,znn,add,dsepgreenzaxi2par,sum,dsepdmasaxi
      integer n
ccc
      sum=dsepgreenzaxi2par(zhat,z0hat,Ltilde)
      n=1
      sign=-1.d0
 100  znp=sign*zhat+2.d0*n
      znn=sign*zhat-2.d0*n
      add=sign*(dsepgreenzaxi2par(znp,z0hat,Ltilde)
     &          +dsepgreenzaxi2par(znn,z0hat,Ltilde))
      sum=sum+add ! dimensionless
c      write(*,*) z0hat,Rtilde,n,sum,add,sign
      n=n+1
      sign=-sign
      if(n.lt.5.or.dabs(add).gt.1.d-10*abs(sum)) goto 100
      dsepgreenzaxi2=sum*dsepdmasaxi(dsqrt(dv)*R0tildei,diffhh*z0hat,
     &                 epgraxipower)/npairs_norm
c      write(*,*) z0hat,Rtilde,dsepgreenzaxi2
      return
      end
ccc
ccc
ccc
      real*8 function dsepgreenzaxi2par(zhat,z0hat,Ltilde)
ccc dimensionless 
      implicit none
      real*8 zhat,z0hat,Ltilde,d2z
      d2z=(zhat-z0hat)**2*Ltilde**2
      if(d2z.lt.1.d3) then
        dsepgreenzaxi2par=Ltilde*dexp(-d2z/4.d0)
      else
        dsepgreenzaxi2par=0.d0
      endif
      return
      end
