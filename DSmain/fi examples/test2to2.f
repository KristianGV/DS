
c_______________________________________________________________________
c                   TEST PROGRAM
c
c
c_______________________________________________________________________
      program test
        implicit none
        real*8 zero,T,x1,x2,xX,eta1,eta2,g1,g2,c12,dsfi2to2ab,x,
     &TR,Tmin,dsfi2to2rhs,dsf_int,a,b,err,dsfi2to2oh2,etaX
        real*8, external :: dsfi2to2int_simp


        include 'dsficom.h'

        zero=0; stat=1

        x=1E-3;T=1E4;x1=1;x2=1;eta1=1;eta2=1;g1=1;g2=1;c12=1
        TR=1E9; Tmin=10;a=1E-4;b=1;err=1E-3;mdm=1;xX=1;etaX=0
        
        write(*,*) "Oh2=", dsfi2to2oh2(TR,Tmin,eta1,eta2,etaX,g1,g2,c12,c12)

      end program