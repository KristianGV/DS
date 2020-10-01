      subroutine dssenu_ctabffcreate(wh,i)

***********************************************************************
*** Creates tabulated capture rates (apart from couplings)
*** Input: wh = 'su' or 'ea' for sun or earth
***        i = table number to store the results in
*** Author: Joakim Edsjo
*** Date: 2015-06-12
***********************************************************************

      implicit none

      include 'dssecap.h'
      include 'dshmcom.h'
      real*8 dssenu_capsunnumff,
     &  gps,gns,gpa,gna,mx,rho,tmp1,tmp2,tmp3
      integer i,j
     
      character*2 wh

      gps=0.0d0
      gns=0.0d0
      gpa=0.0d0
      gna=0.0d0
      rho=0.3d0 ! Default local halo density
      do j=0,ncff
        mx=10**(dble(j)*5.0d0/dble(ncff))
        write(*,80) '  Taking care of mass number ',j,'/',ncff,
     &       ' (',mx,' GeV)'
 80     format(A,I4,A,I4,A,F11.4,A)
        
        if (wh.eq.'su'.or.wh.eq.'SU') then
           gpa=0.d0
           gna=0.d0

           gps=1.d0
           gns=0.d0
           tmp1=dssenu_capsunnumff(mx,rho,gps,gns,gpa,gna)
           ctabffsu(j,i,1)=tmp1
           gps=0.0d0
           gns=1.d0
           tmp2=dssenu_capsunnumff(mx,rho,gps,gns,gpa,gna)
           ctabffsu(j,i,2)=tmp2
           gps=1.0d0 ! cross term
           gns=1.d0
           tmp3=dssenu_capsunnumff(mx,rho,gps,gns,gpa,gna)
           ctabffsu(j,i,3)=tmp3-tmp1-tmp2
           write(*,90) '    result (s): ',tmp1,tmp2,tmp3-tmp1-tmp2

           gps=0.d0
           gns=0.d0

           gpa=1.d0
           gna=0.d0
           tmp1=dssenu_capsunnumff(mx,rho,gps,gns,gpa,gna)
           ctabffsu(j,i,4)=tmp1
           gpa=0.0d0
           gna=1.d0
           tmp2=dssenu_capsunnumff(mx,rho,gps,gns,gpa,gna)
           ctabffsu(j,i,5)=tmp2
           gpa=1.0d0 ! cross term
           gna=1.d0
           tmp3=dssenu_capsunnumff(mx,rho,gps,gns,gpa,gna)
           ctabffsu(j,i,6)=tmp3-tmp1-tmp2
           write(*,90) '    result (a): ',tmp1,tmp2,tmp3-tmp1-tmp2
 90        format(A,3(1x,E14.6))
        else
c...Earth not implemented
c          sigsi=1.0d-40
c          tmp1=dssenu_capearthnum(mx,sigsi)
c          ctabea(j,i)=tmp1
c          write(*,*) '    result: ',tmp1
        endif
      enddo

      return

      end


      

