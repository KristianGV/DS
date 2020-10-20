c_______________________________________________________________________
c  Function dsfi2to2int is the integrand to be integrated over s=pcm^2.
c 
c     input:
c       x=1/s - variable to integrate over
c       T - bath temperature
c       xi = mi/T 
c       etai=\pm exp^(mu_i/T)  - parameter in phase-space distribution
c          pluss (minus) for bosons (fermions)    
c       stat - determines if quantum corrections is to be considered 
c         (stat=0 for Maxwell-Boltzmann statistics) (stat\neq 0 for FD/BE)
c
c
c  author: Kristian Gjestad Vangsnes (kristgva@uio.no)   2020-09-07
c=======================================================================
      real*8 function dsfi2to2int(x,T,x1,x2,xX,eta1,eta2,etaX)
      implicit none
      include 'dsmpconst.h'
      include 'dsficom.h'

      real*8 T,x1,x2,xX,eta1,eta2,etaX,zero,k1,k1X,k1bess,s,x,p_ij,sigma,
     &dsbessek1,m1,m2,sstar,p,dssigmavpartial
      zero=0
      s=1/x
      ! sstar=1/x
      ! s=sstar*4*mdm**2

      call dsfik1(sqrt(s)/T,x1,x2,zero,eta1,eta2,k1)
      call dsfik1(sqrt(s)/T,x1,xX,zero,zero,etaX,k1X)
      call dsfik1(sqrt(s)/T,x1,xX,zero,zero,zero,k1bess)

      m1=x1*T;m2=x2*T

      if(stat.eq.0) then
    !     dsfi2to2int=sigma(s)*p_ij(mdm,mdm,s)**2*sqrt(s)*
    !  &dsbessek1(sqrt(s)/T)/exp(sqrt(s)/T)*s**2

        p=sqrt(s/4-mdm**2)
        if(isnan(dssigmavpartial(ichannel_22,p))) then
          dsfi2to2int=dssigmavpartial(ichannel_22,p)
          dsfi2to2int=0
    !       write(*,*) "dssigmavpartial is NaN at sqrts=", sqrt(s),
    !  &"for channel=",ichannel_22 
        else
        dsfi2to2int=dssigmavpartial(ichannel_22,p)*p_ij(mdm,mdm,s)*
     &(s-2*mdm**2)*dsbessek1(sqrt(s)/T)/exp(sqrt(s)/T)*s**2/gev2cm3s/2
        end if

      else
        if(s.lt.4*mdm**2) then
          dsfi2to2int=0
        else if(etaX.eq.0) then
          dsfi2to2int= s**2*p_ij(mdm,mdm,s)**2*sigma(s)*sqrt(s)*k1
        else
          dsfi2to2int=s**2*p_ij(mdm,mdm,s)**2*sigma(s)*sqrt(s)*k1*k1X/k1bess
        end if
      end if
      return
      end 

c     Function that calculated p_ij that is used in the annihilation rate      
      real*8 function p_ij(mi,mj,s)
      implicit none
      real*8 s,mi,mj

      p_ij=sqrt(s-(mi+mj)**2)*sqrt(s-(mi-mj)**2)/(2*sqrt(s))
      
      return
      end

c     dsfi2to2int made into a function of one dimension used to integrate over    
      real*8 function dsfi2to2int_simp(lgx)
      implicit none

      include 'dsficom.h'

      real*8 k1,s,x,dsfi2to2int,lgx
      x=exp(lgx)

      dsfi2to2int_simp=dsfi2to2int(x,T_22,x1_22,x2_22,xX_22,eta1_22
     &,eta2_22,etaX_22)*x
      return
      end


