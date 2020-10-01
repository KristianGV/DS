c_______________________________________________________________________
c                   TEST PROGRAM
c
c
c_______________________________________________________________________
      program test
        implicit none
        include 'dsmpconst.h'

        real*8 dsfistat,P,x1,x2,x3,eta1,eta2,eta3,k1,dsfidecint_simp
        real*8 tmp, m, eta, dsfidecint,dsfidecab, TR, w, g,dsfidecoh2
     &,temp2
        integer stat, i, temp
        character*80 filename
        common /dsfidecvar/ M,eta
        common /dsfik1var/ x1,x2,x3,eta1,eta2,eta3
        common /global_fi/ stat

        P=4; x1=3;x2=1; x3=1;eta1=1; eta2=1; eta3=1; stat=0
        tmp=1.  ; M=1E4; eta=1
        TR=1.D9; w=1.E-18; g=1

        filename='dydt.dat'

        open(unit=100,file=filename,status='unknown',form='formatted')
        write(100,'(A)') '   T [GeV]                  dY/dlog(T)'
        do temp=0,7000,1
          temp2=real(temp)/1000
          write(100,*) 10**temp2,w*g*M*M*dsfidecint_simp(log(10**temp2))
     &*tmp*log(10._8)/(pi*pi)
        end do
        close(100)
      end program
