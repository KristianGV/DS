      subroutine dsddsigmanucleon(v,e,sigsip,sigsin,sigsdp,sigsdn,ierr)
c_______________________________________________________________________
c  dark matter - nucleon cross section.
c
c  type : commonly used
c  desc : Calculate nuclear cross sections
c
c  common:
c    'dssusy.h' - file with susy common blocks
c  output:
c    sigsip, sigsin : proton and neutron spin-independent cross sections
c    sigsdp, sigsdn : proton and neutron spin-dependent cross sections
c    units: cm^2
c  author: paolo gondolo (gondolo@lpthe.jussieu.fr) 1994,1995,2002
c     13-sep-94 pg no drees-nojiri twist-2 terms
c     22-apr-95 pg important bug corrected [ft -> ft mp/mq]
c     06-apr-02 pg drees-nojiri treatment added
c     06-feb-16 pg renamed dsddsigmanucleon (was dsddneunuc)
c     16-mar-19 tb bug fixes (catch errors independently)
c=======================================================================
      implicit none

      real*8 v,e,sigsip,sigsdp,sigsin,sigsdn
      real*8 sigij(27,27)
      integer ierr, ierr2

c...Initialize
      sigsin=0.d0
      sigsdn=0.d0
      sigsip=0.d0
      sigsdp=0.d0
      ierr=0

      call dsddsigma(v,e,1,0,sigij,ierr)
      if (ierr.eq.0) then
        sigsin=sigij(1,1)
        sigsdn=sigij(4,4)
      endif  
      call dsddsigma(v,e,1,1,sigij,ierr2)
      if(ierr.eq.0) then
        sigsip=sigij(1,1)
        sigsdp=sigij(4,4)
      endif
      ierr = ierr + 100*ierr2
      return
      end
