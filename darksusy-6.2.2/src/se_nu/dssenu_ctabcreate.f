      subroutine dssenu_ctabcreate(wh,i)

***********************************************************************
*** Creates tabulated capture rates (apart from cross section)
*** Input: wh = 'su' or 'ea' for sun or earth
***        i = table number to store the results in
*** Author: Joakim Edsjo
*** Date: 2003-11-27
***********************************************************************

      implicit none

      include 'dssecap.h'
      include 'dshmcom.h'
      real*8 dssenu_capearthnum,dssenu_capsunnum,
     &  sigsi,sigsd,mx,tmp1,tmp2,rho
      integer i,j
     
      character*2 wh

      sigsi=0.0d0
      sigsd=0.0d0
      rho=0.3d0
      do j=0,nc
        mx=10**(dble(j)*5.0d0/dble(nc))
        write(*,80) '  Taking care of mass number ',j,'/',nc,
     &       ' (',mx,' GeV)'
 80     format(A,I4,A,I4,A,F11.4,A)
        
        if (wh.eq.'su'.or.wh.eq.'SU') then
          sigsi=1.0d-40
          sigsd=0.0d0
          tmp1=dssenu_capsunnum(mx,rho,sigsi,sigsd)
          ctabsusi(j,i)=tmp1
          sigsi=0.0d0
          sigsd=1.0d-40
          tmp2=dssenu_capsunnum(mx,rho,sigsi,sigsd)
          ctabsusd(j,i)=tmp2
          write(*,*) '    result: ',tmp1,tmp2
        else
          sigsi=1.0d-40
          tmp1=dssenu_capearthnum(mx,rho,sigsi)
          ctabea(j,i)=tmp1
          write(*,*) '    result: ',tmp1
        endif
      enddo

      return

      end


      

