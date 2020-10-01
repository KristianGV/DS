c_______________________________________________________________________
c                   GENERIC FIMP DECAY
c
c
c_______________________________________________________________________
      program two
        implicit none
        include 'dsmpconst.h'
        include 'dsficom.h'


        character*80 filename

        real*8 T,x1,x2,xX,eta1,eta2,etaX,g1,g2,c12,zero,k1,s,x,
     &step,a,b,f,dsfi2to2int_simp,y,dsfidecoh2,TR,Tmin,w,g,mdm,g_dm
     &,width
        integer i,length


        TR=1.D5; mdm=10; stat=0; g_dm=1.D-12;eta_dec=1;stat=0;g=1

        filename='generic_FIMP.dat'

        a=10;b=1E3

        open(unit=100,file=filename,status='unknown',form='formatted')
        write(100,'(A)') '   M_dec                        oh2'
        length=100
        do i=0,length
          step=real(i)/real(length)
          M_dec=f(step,a,b)
          w=width(M_dec,mdm,g_dm)
          write(100,*) M_dec, dsfidecoh2(TR,w,g,mdm)
        end do
        close(100)
      end program

      real*8 function f(x,a,b)
      implicit none
      real*8 x,a,b
      f=(b-a)*x+a
      return
      end

      real*8 function width(M_dec,mdm,g_dm)
      implicit none
      include 'dsmpconst.h'
      real*8 M_dec,mdm,g_dm,p,lambda
      lambda=1E-3*sqrt(M_dec/mdm)  
      p=sqrt(M_dec**2/4-mdm**2)
      width=lambda**2*M_dec/(8*pi)
      ! g_dm**2*p**3/(8*pi*M_dec**2)
      return
      end
