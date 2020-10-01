      subroutine dsddg2sigma(mwimp,swimp,v,er,a,z,gg,sigij,ierr)
c_______________________________________________________________________
c  Calculate the unpolarized *equivalent* scattering cross section
c      sigma = E_{max} d\sigma/dE_R
c  The calculation depends on the definition of the couplings G_a, i.e.,
c    the dd scheme,
c      sigma = \sum_{ij} sig_{ij}
c    where sig_{ij} = \sum_{NN'} G_i^{N*} G_j^{N'} P_{ij}^{NN'}.
c  input:
c    mwimp : real*8: WIMP mass in GeV
c    swimp : real*8: WIMP spin in hbar
c    v   : real*8  : WIMP-nucleus relative velocity in km/s
c    er  : real*8  : nucleus recoil energy in keV
c    a,z : integer : nucleus mass number and atomic number, respectively
c    gg  : real*8 array gg(27,2)  : list of couplings G_i^N
c  output:
c    sigij : real*8  : sig_{ij} in cm^2
c    ierr  : integer : error flag (0=no error)
c            1-10  : general error
c            11-20 : error in SI part      
c            21-30 : error in SD part
c            31-40 : error in SI and SD part
c     Note: if there is an error in any part, thos cross sections are
c     put to zero.
c
c
c  20161113 Paolo Gondolo   first version
c  20161119 Paolo Gondolo   second version
c=======================================================================
      implicit none

      include 'dsddcom.h'
      include 'dsnuclides.h'
      include 'dsmpconst.h'
      include 'dsio.h'

      real*8 mwimp,swimp,v,er
      integer a,z,ierr
      complex*16 gg(ddng,2)
      real*8 sigij(ddng,ddng)
      integer p,n
      parameter (p=1)
      parameter (n=2)
      integer i,j,jerr
      integer ii,dsnucldindx
      real*8 mx,mni,muxi,q,kin,vperpsq

      real*8 aux
      complex*16 pp(2,2)

c      do i=1,ddng
c        write (*,'(a,i2,a,"(",g12.5,",",g12.5,")",a,i2,a,"(",g12.5,",",g12.5,")")') 'PG> gg(',i,',1) = ',gg(i,1),'     gg(',i,',2) = ',gg(i,2)
c      enddo

      ierr=0
      do i=1,ddng
        do j=1,ddng
          sigij(i,j)=0.d0
        enddo
      enddo
      
      mx=mwimp
      ii=dsnucldindx(a,z)
      if (ii.eq.0) then
         if (prtlevel.gt.0) write(*,*)
     &       'DS Warning in dsddg2sigma: non-existent nuclide (',
     &       a,',',z,')'
        ierr=-1
        return
      endif
      mni = nucldm(ii)
      muxi = mni*mx/(mni+mx)
      q = dsqrt(2.d0*mni*er*1.d-6)
      kin=muxi**2/pi*gev2cm2
      vperpsq=(v/c_light)**2-(0.5d0*q/muxi)**2
c      write (*,*) 'PG> vperpsq=',vperpsq 
      if (vperpsq.lt.0.d0) return

      do i=1,ddng
        do j=1,ddng
           aux=abs(conjg(gg(i,p))*gg(j,p))+
     &         abs(conjg(gg(i,p))*gg(j,n))+
     &         abs(conjg(gg(i,n))*gg(j,p))+
     &         abs(conjg(gg(i,n))*gg(j,n))
          if (aux.ne.0.d0) then            
            call dsddpp(swimp,vperpsq,q,a,z,i,j,pp,jerr)
            if (jerr.ne.0) then
              ierr=jerr
            else
              sigij(i,j)=kin*real(
     &             conjg(gg(i,p))*gg(j,p)*pp(p,p)+
     &             conjg(gg(i,p))*gg(j,n)*pp(p,n)+
     &             conjg(gg(i,n))*gg(j,p)*pp(n,p)+
     &             conjg(gg(i,n))*gg(j,n)*pp(n,n))
            endif
          endif
        enddo
      enddo

c... TB debug
c      write(*,*) 'gg(1,p)),gg(1,n) = ',gg(1,p),gg(1,n)
c      write(*,*) 'gg(4,p)),gg(4,n) = ',gg(4,p),gg(4,n)
c      write(*,*) 'pp(p,p),pp(p,n),pp(n,p),pp(n,n) = ', pp(p,p),pp(p,n),pp(n,p),pp(n,n)
        
      return
      end
