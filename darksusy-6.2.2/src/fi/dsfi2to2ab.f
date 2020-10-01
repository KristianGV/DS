c_______________________________________________________________________
c  Function dsfi2to2ab calculates abundance of FIMP's Y_dm
c   
c    input:
c       TR - reheating temperature
c       Tmin - minimum temperature
c       etai=\pm exp^(mu_i/T)  - parameter in phase-space distribution
c          pluss (minus) for bosons (fermions)    
c        gi - internal degrees of freedom for particle i
c        c12 - constant = 1/2 if particle 1= particle 2, else c12=1.
c
c
c  author: Kristian Gjestad Vangsnes (kristgva@uio.no)   2020-09-07
c=======================================================================
      real*8 function dsfi2to2ab(TR,Tmin,m1,m2,eta1,eta2,etaX,g1,g2,c12)
      implicit none

      include 'dsmpconst.h'
      include 'dsficom.h'

      real*8 TR,Tmin,eps,a,b,dsf_int,sum,g1,g2,c12,m1,m2,eta1,eta2,etaX ,
     &dsf_int2,abserr,resabs,resasc,res,prec
      real*8, external :: dsfi2to2rhs

      eta1_22=eta1;eta2_22=eta2;etaX_22=etaX;g1_22=g1;g2_22=g2
      c12_22=c12;m1_22=m1;m2_22=m2
c     Integration limits      
c      a=log(Tmin);b=log(TR);eps=1E-5
      ! a=Tmin;b=TR;eps=1.D-5
      a=log(Tmin);b=log(TR);eps=1.d-5
c     Integrate using adaptive Gaussian integration routine      
      call dgadap(a,b,dsfi2to2rhs,eps,sum)
      ! sum=dsf_int2(dsfi2to2rhs,a,b,eps)
      ! call dqk21(dsfi2to2rhs,a,b,sum,eps,resabs,resasc)

      ! prec=1.d-2
      ! call dsfun_int(dsfi2to2rhs,a,b,eps,prec,res)
      if(stat.eq.0) then
            dsfi2to2ab=sum*g1_22*g2_22*c12_22/(8*pi**4)
      else
            dsfi2to2ab=sum*g1_22*g2_22*abs(eta1_22*eta2_22)
     &*c12_22/(8*pi**4)
      end if
      return
      end