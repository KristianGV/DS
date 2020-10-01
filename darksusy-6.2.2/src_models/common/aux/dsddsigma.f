      subroutine dsddsigma(v,e,a,z,sigij,ierr)
c_______________________________________________________________________
c
c  type : INTERFACE
c  desc : UNpolarized *equivalent* WIMP nucleus cross section including form factors
c      sigma = E_{max} d\sigma/dE_R
c  The calculation depends on the definition of the couplings G_a, i.e.,
c    the dd scheme,
c      sigma = \sum_{ij} sig_{ij}
c    where sig_{ij} = \sum_{NN'} G_i^{N*} G_j^{N'} P_{ij}^{NN'}.
c  input:
c    v : real*8    : WIMP-nucleus relative velocity in km/s
c    e : real*8    : nucleus recoil energy in keV
c    a : integer   : nucleus mass number
c    z : integer   : nucleus atomic numbers
c  output:
c    sigij(27,27) : real*8  : partial cross section array in cm^2
c                      note that the usual spin-independent scattering
c                      cross section is sigij(1,1) and the usual 
c                      spin-dependent one is sigij(4,4)
c    ierr  : integer : error code (0=no error)
c
c  Note: If you call this routine with arguments
c    0.0d0,0.0d0,1,1 etc you get the usual scattering cross sections
c       on protons      
c    0.0d0,0.0d0,1,0 etc you get the usual scattering cross sections
c       on neutrons
c
c  Note: This particular version of dsddsigma can be used (as an interface
c        function) by all particle modules that provide the function dsddgpgn.
c        Alternatively, a particle module can provide its own version of
c        dsddsigma, in which case this version will be overriden.
c
c  author: paolo gondolo (paolo.gondolo@utah.edu) 2016
c  mod: torsten bringmann (generalized and placed in /common, and hence
c                          usable by all particle modules) 2018
c=======================================================================
      implicit none
      include 'dsddcom.h'
      include 'dsio.h'
      real*8 v,e,sigij(ddng,ddng)
      integer a,z,ierr
      integer i,j
      real*8 mx,sx, dsmwimp, dsdmspin
      complex*16 gg(ddng,2)
      
      ierr=0
      mx = dsmwimp()
      sx = dsdmspin()
 
      call dsddgpgn(gg,ierr)
c      write (*,*) 'TB debug : ', gg(1,1), gg(1,2), gg(4,1), gg(4,2)
      if (ierr.ne.0) then
        if (prtlevel.gt.0) 
     &     write(*,*)'WARNING from dsddsigma: error in nucleon couplings!'
        ierr=0 ! we still continue here... 
      endif
      call dsddg2sigma(mx,sx,v,e,a,z,gg,sigij,ierr)
      if (ierr.ne.0) then
        do i=1,27
           do j=1,27
              sigij(i,j)=0.d0
           enddo
        enddo
        return
      endif

      return
      end
