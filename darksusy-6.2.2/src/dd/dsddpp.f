      subroutine dsddpp(swimp,vperpsq,q,a,z,i,j,pp,ierr)
c_______________________________________________________________________
c  Structure functions for direct detection.
c  input:
c    swimp   : real*8: WIMP spin in hbar
c    vperpsq : real*8  : vperp^2 in c^2
c    q       : real*8  : momentum transfer in GeV ( q=sqrt(2*M*E) )
c    a       : integer : mass number
c    z       : integer : atomic number
c    i       : integer : index of P_{ij}^{NN'}
c    j       : integer : index of P_{ij}^{NN'}
c  output:
c    pp : complex*16 : array pp(2,2) containing P_{ij}^{NN'}
c    ierr : integer : error flag (0=no error)
c
c  20161119 Paolo Gondolo   second version
c=======================================================================
      implicit none
      include 'dsddcom.h'
      include 'dsio.h'
      real*8 swimp,vperpsq,q
      integer a,z
      complex*16 pp(2,2)
      integer ierr

      integer p,n
      parameter (p=1)
      parameter (n=2)
      integer i,j,jerr
      real*8 ffsi,ffsdpp,ffsdpn,ffsdnn
      real*8 cspin
      
      ierr=0
      pp(p,p)=cmplx(0.d0,0.d0)
      pp(p,n)=cmplx(0.d0,0.d0)
      pp(n,p)=cmplx(0.d0,0.d0)
      pp(n,n)=cmplx(0.d0,0.d0)

      if (i.eq.1.and.j.eq.1) then
        if (
     &       ddsf(1)(2:).eq.'best'
     &       .or.ddsf(1)(2:).eq.'FB'
     &       .or.ddsf(1)(2:).eq.'SOG'
     &       .or.ddsf(1)(2:).eq.'Fermi'
     &       .or.ddsf(1)(2:).eq.'L-S'
     &       .or.ddsf(1)(2:).eq.'ds4.1'
     &       .or.ddsf(1)(2:).eq.'gauss'
     &       .or.ddsf(1)(2:).eq.'gould'
     &       ) then
          call dsddffsi(q,a,z,ffsi,jerr)
          if (jerr.ne.0) then
            write (*,*) 'Warning dsddpp : error in F_m structure functions ',jerr
            ierr=100+jerr
            return
          endif
          pp(p,p)=cmplx(z*z*ffsi,0.d0)
          pp(p,n)=cmplx(z*(a-z)*ffsi,0.d0)
          pp(n,p)=pp(p,n)
          pp(n,n)=cmplx((a-z)*(a-z)*ffsi,0.d0)
        else if (
     &         ddsf(1)(2:).eq.'haxton'
     &         ) then
          write (*,'(a,a,a)') 'Warning dsddpp: ''sf_m''=''',trim(ddsf(1)(2:)),''' not available'
        else
          write (*,'(a,a,a)') 'Warning dsddpp: ''sf_m''=''',trim(ddsf(1)(2:)),''' not available'
        endif
      endif
        
      if (i.eq.4.and.j.eq.4) then
        if (
     &       ddsf(2)(2:).eq.'best'
     &       .or.ddsf(2)(2:).eq.'ISM'
     &       .or.ddsf(2)(2:).eq.'OddG'
     &       .or.ddsf(2)(2:).eq.'SPSM'
     &       .or.ddsf(2)(2:).eq.'SIMS'
     &       .or.ddsf(2)(2:).eq.'ISMR'
     &       ) then
          call dsddffsd(q,a,z,ffsdpp,ffsdnn,ffsdpn,jerr)
          if (jerr.ne.0) then
            write (*,*) 'Warning dsddpp : error in F_sigma structure functions ',jerr
            ierr=200+jerr
            return
          endif
          cspin=4.d0*swimp*(swimp+1.d0)/3.d0
          pp(p,p)=cmplx(cspin*4.d0*ffsdpp,0.d0)
          pp(p,n)=cmplx(cspin*2.d0*ffsdpn,0.d0)
          pp(n,p)=pp(p,n)
          pp(n,n)=cmplx(cspin*4.d0*ffsdnn,0.d0)
        else if (
     &         ddsf(1)(2:).eq.'haxton'
     &         ) then
          write (*,'(a,a,a)') 'Warning dsddpp: ''sf_sigma''=''',trim(ddsf(2)(2:)),''' not available'
        else
          write (*,'(a,a,a)') 'Warning dsddpp: ''sf_sigma''=''',trim(ddsf(2)(2:)),''' not available'
        endif
      endif

      return
      end
