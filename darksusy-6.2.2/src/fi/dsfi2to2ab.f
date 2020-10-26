c_______________________________________________________________________
c  Function dsfi2to2ab calculates abundance of FIMP's Y_dm
c   
c    input:
c       Tmin - minimum temperature
c       TR - reheating temperature
c       mi - mass of particle i
c       etai=\pm exp^(mu_i/T)  - parameter in phase-space distribution
c          pluss (minus) for bosons (fermions)    
c       gi - internal degrees of freedom for particle i
c       c12 - constant = 1/2 if particle 1= particle 2, else c12=1.
c
c
c  author: Kristian Gjestad Vangsnes (kristgva@uio.no)   2020-09-07
c=======================================================================
      real*8 function dsfi2to2ab(Tmin,TR,m1,m2,eta1,eta2,etaX,g1,g2,c12)
      implicit none

      include 'dsmpconst.h'
      include 'dsficom.h'

      real*8 TR,Tmin,eps,a,b,dsf_int,sum,g1,g2,c12,m1,m2,eta1,eta2,etaX ,
     &dsf_int2,abserr,resabs,resasc,res,prec
      real*8, external :: dsfi2to2rhs

c     Common block variables are set
      eta1_22=eta1;eta2_22=eta2;etaX_22=etaX;g1_22=g1;g2_22=g2
      c12_22=c12;m1_22=m1;m2_22=m2

c     Integration limits      
      ! a=Tmin;b=TR;eps=1.d-5
      a=log(Tmin);b=log(TR);eps=1.d-5

c     Integrate using integration routine    
      ! call dqk21(dsfi2to2rhs,a,b,sum,abserr,resabs,resasc)  
      ! call dgadap(a,b,dsfi2to2rhs,eps,sum)
      sum=dsf_int2(dsfi2to2rhs,a,b,eps)
      ! write(*,*) "abserr=",abserr
      ! prec=1.d-2
      ! call dsfun_int(dsfi2to2rhs,a,b,eps,prec,res)
      if(stat.eq.0) then
c       USES MY OWN ROUTINE      
        ! dsfi2to2ab=2*sum*g1_22*g2_22*c12_22/(8*pi**4)
c       USES <sv>
        dsfi2to2ab=2*sum*g1_22*g2_22*c12_22*mdm**4/(4*pi**4)
      else
        dsfi2to2ab=2*sum*g1_22*g2_22*abs(eta1_22*eta2_22)
     &*c12_22/(8*pi**4)
      end if
      return
      end