
c_______________________________________________________________________
c                   TEST PROGRAM
c
c
c_______________________________________________________________________
      program test
        implicit none

        include 'dsficom.h'

        real*8 dsfistat,P,x1,x2,x3,eta1,eta2,eta3,k1,dsfidecint_simp
        real*8 tmp, m, eta, dsfidecint,dsfidecab, TR, w, g,dsfidecoh2
        real*8 dsfi2to2int_simp,x1_2,x2_2,eta1_2,eta2_2,g1,g2,c12,T,mdm
        real*8 dsfi2to2int,x,zero,dsfi2to2rhs,Tmin,dsfi2to2ab,dsbessek1
        integer  i


        zero=0
c       Note M>m1+m2 => x1>x2+x3
        x=1;P=1;x1=3;x2=1; x3=0;eta1=0; eta2=1; eta3=1; stat=0
c       General variables
        Tmin=100; TR=1.D10

c       1->2 variables
        tmp=1; M_dec=1E2; eta_dec=1; w=1; g=1; mdm=10
    

        write(*,*) "S=",
     & dsfistat(P,x1,x2,x3,eta1,eta2)
        call  dsfik1(x,x2,x3,eta1,eta2,eta3,k1)
        write(*,*) "K1 =", k1
      

        write(*,*) "INT", dsfidecint(log(tmp),M_dec,eta_dec)
        write(*,*) "Y=", dsfidecab(TR,w,g)
        write(*,*) "Oh^2=", dsfidecoh2(TR,w,g,mdm)
      end program