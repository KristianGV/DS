c_______________________________________________________________________
c                   SCALAR FIMP DECAY
c
c
c_______________________________________________________________________
      program scalar_FIMP

        implicit none

        include 'dsmpconst.h'
        include 'dsparticles.h'
        include 'dssm.h'
        include 'dsidtag.h'
        include 'dsficom.h'

        real*8 T,x1,x2,xX,eta,eta1,eta2,etaX,g1,g2,c12,zero,k1,s,x,
     &step,mmin,mmax,f,dsfi2to2int_simp,dsfidecoh2,TR,Tmin,M,w,g
     &,inputlambda,dsgammahpartial,dsmwimp,lgmdm,dsfi2to2oh2,dsfi2to2rhs
     &,rddec,rd2to2,sigma,sig,rd2to2_tmp,sum,tmp,dsfi2to2ab,lgtmp,
     &Y,Y_dec,Y_22,dsfidecab,dsfidecint,lgs
      
      real*8 lb,ub,k1_ds,k1_my,dsbessek1
        character*80 filename,filename1,filename_debug,filename2,
     &filename3

        integer i,length,ierr,iwarn,ichannel,j

          
        call dsinit
        write (*,*) '-------------------------------------------------------'
        write (*,*) 

        TR=1.D6; ; inputlambda=1.D-11;eta=1;stat=0;g=1;
        eta_dec=1; zero=0; Tmin=0

        
        eta1_22=1;eta2_22=1;etaX_22=0;g1_22=1;g2_22=1;c12_22=1.d0
        selfcon=1


        filename='scalar_FIMP.dat'
        filename1='scalar_FIMP_Y.dat'
        filename2='scalar_FIMP_decay.dat'
        filename3='scalar_FIMP_2to2.dat'
        filename_debug='debug.dat'

        M_dec=mass(khsm)
        m1_22=M_dec;m2_22=M_dec

        length=100
        
c ... scale setup
        mmin=1.d-3;mmax=1.d3

        open(unit=100,file=filename,status='unknown',form='formatted')
        write(100,'(A)') '   mdm                       oh2'

        
        open(unit=200,file=filename1,status='unknown',form='formatted')
        write(200,'(A)') '   tmp                       Y'

        open(unit=400,file=filename2,status='unknown',form='formatted')
        write(400,'(A)') '   mdm                       oh2'

        open(unit=500,file=filename3,status='unknown',form='formatted')
        write(500,'(A)') '   mdm                       oh2'

        open(unit=300,file=filename_debug,status='unknown',
     &form='formatted')
        write(300,*) '%%%%%%%%%%%%%%    DEBUG FILE    %%%%%%%%%%%%%%'


        do i=0,length
          if(i.ne.0) then
            do j=1,10
              if(real(i)/real(length).eq.real(j)/real(10)) then
                write(*,*) j*10,'%'
              end if
            end do
            
          end if
          step=real(i)/real(length)

c       RELIC ABUNDACE TOTAL        

    !       lgmdm=f(step,log(mmin),log(mmax))
    !       mdm=exp(lgmdm)
    !       Tmin=1.d-3
    !       call dsgivemodel_silveira_zee(inputlambda,mdm)
    !       ichannel=19
    !       w=dsgammahpartial(ichannel,zero) 
                 

    !       rddec=dsfidecoh2(TR,Tmin,w,M_dec,
    !  &g,eta_dec)


    !       rd2to2=0.d0
    !       do ichannel_22=1,18
    !         rd2to2_tmp=dsfi2to2oh2(TR,Tmin,m1_22,m2_22
    !  &,eta1_22,eta2_22,etaX_22,g1_22,g2_22,c12_22)
    !         rd2to2=rd2to2+rd2to2_tmp
    !       end do


    !       write(100,*) mdm, (rddec+rd2to2)


c       RELIC ABUNDACE 2to2        

          lgmdm=f(step,log(mmin),log(mmax))
          mdm=exp(lgmdm)
          Tmin=1.d-3
          call dsgivemodel_silveira_zee(inputlambda,mdm)

          rd2to2=0
    !       dsfi2to2oh2(TR,Tmin,m1_22,m2_22
    !  &,eta1_22,eta2_22,etaX_22,g1_22,g2_22,c12_22)
          do ichannel_22=1,18
            rd2to2_tmp=dsfi2to2oh2(TR,Tmin,m1_22,m2_22
     &,eta1_22,eta2_22,etaX_22,g1_22,g2_22,c12_22)
            rd2to2=rd2to2+rd2to2_tmp
          end do

          write(500,*) mdm, rd2to2

c     RELIC ABUNDANCE DECAY

    !       lgmdm=f(step,log(mmin),log(mmax))
    !       mdm=exp(lgmdm)
    !       Tmin=1.d-3
    !       ichannel=19
    !       call dsgivemodel_silveira_zee(inputlambda,mdm)
    !       w=dsgammahpartial(ichannel,zero) 
                 

    !       rddec=dsfidecoh2(TR,Tmin,w,M_dec,
    !  &g,eta_dec)
    !       write(400,*) mdm, rddec



c     CALC Y

    !       Tmin=1
    !       lgtmp=f(step,log(Tmin),log(TR))
    !       tmp=exp(lgtmp)
    !       mdm=1
    !       call dsgivemodel_silveira_zee(inputlambda,mdm)
    !       w=dsgammahpartial(ichannel,zero) 
    !       Y_dec=dsfidecab(TR,tmp,w,M_dec,g,eta_dec)
    !       Y_22=0
    !       do ichannel_22=1,18
    !         Y_22=Y_22+dsfi2to2ab(TR,tmp,m1_22,m2_22
    !  &,eta1_22,eta2_22,etaX_22,g1_22,g2_22,c12_22)
    !       end do
    !       Y=Y_dec+Y_22
    !       write(200,*) tmp, Y


c       TESTING dsfi2to2rhs

          ! mdm=25
          ! Tmin=20000.d0

          ! tmp=f(step,Tmin,TR)
          ! x1_22=m1_22/tmp
          ! x2_22=x1_22
          ! sum=0

          ! do ichannel_22=1,18
          !   sum=sum+dsfi2to2rhs(tmp)
          ! end do

          ! write(300,*) tmp, sum

c       TESTING dsfi2to2int_simp

          ! T_22=1.d4
          ! x1_22=m1_22/T_22
          ! x2_22=x1_22

          ! mdm=1.d2
          ! sum=0.d0
          ! call dsgivemodel_silveira_zee(inputlambda,mdm)

          ! x=f(step,4*mdm**2+1.d-30,4*mdm**2+1.d10)
          ! do ichannel_22=1,18
          !   sum=sum+dsfi2to2int_simp(x)   
          ! end do  

          ! write(300,*) x, sum

c       TESTING dsfidecint


          ! mdm=20
          ! Tmin=mdm/100

          ! tmp=f(step,Tmin,TR)          
          ! write(300,*) tmp, dsfidecint(tmp,M_dec,eta_dec)


c     TESTING CROSS SECTION VALUES          

        ! mdm=1.d-3
        ! call dsgivemodel_silveira_zee(inputlambda,mdm)
        ! !  lgs=f(step,log((2*mdm)**2),log(1.D10))
        ! !  s=exp(lgs)  

        ! s=f(step,m1_22**2-2,m1_22**2+2)

        ! sig=0.0
        ! do ichannel_22=1,18
        !   sig=sig+sigma(s)
        ! end do
        ! write(300,*) s, sig

c     Testing Bessel function

      ! mdm=1.d-2; lb=2*mdm; ub=50.d0; stat=1
      ! x=f(step,lb,ub)
      ! k1_ds=dsbessek1(x)/exp(x)
      ! call dsfik1(x,zero,zero,zero,zero,zero,k1)
      ! k1_my=k1
      ! write(300,*) x,k1_ds, k1_my, k1_ds-k1_my

      end do  
      close(100)
      close(200)
      close(300)
      close(400)
      close(500)
      write(*,*) 'Done!'
      end program

      real*8 function f(x,mmin,mmax)
      implicit none
      real*8 x,mmin,mmax
      f=(mmax-mmin)*x+mmin
      return
      end

