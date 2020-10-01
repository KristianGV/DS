c_______________________________________________________________________
c                   TEST PROGRAM
c
c
c_______________________________________________________________________
      program two
        implicit none
        include 'dsmpconst.h'


        character*80 filename

        real*8 T,x1,x2,xX,eta1,eta2,etaX,g1,g2,c12,zero,k1,s,x,
     &step,a,b,f,dsfi2to2int_simp,y,dsfi2to2rhs,TR,Tmin
        integer stat,i,length

        common /dsfi2to2_var/ T,x1,x2,xX,eta1,eta2,etaX,g1,g2,c12
        common /global_fi/ stat

        T=1.;x1=1;x2=1;eta1=0;eta2=0;g1=1;g2=1;c12=1;stat=1;zero=0
        etaX=0;xX=1;y=1E-3;TR=1E10;Tmin=10
        filename='2to2rhs.dat'
        a=log(Tmin);b=log(TR)
c        a=1E-9; b=1/((x1+x2)*T)**2
        open(unit=100,file=filename,status='unknown',form='formatted')
        write(100,'(A)') '   x                        int'
        length=1000 
        do i=0,length
          step=real(i)/real(length)
          x=f(step,a,b)
          write(100,*) x, dsfi2to2rhs(x)
c          write(100,*) x, dsfi2to2int_simp(x)
        end do
        close(100)
      end program

      real*8 function f(x,a,b)
      implicit none
      real*8 x,a,b
      f=(b-a)*x+a
      return
      end