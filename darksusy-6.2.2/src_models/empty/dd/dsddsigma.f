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
c  author: torsten bringmann, 2018
c=======================================================================
      implicit none
      real*8 v,e,sigij(27,27)
      integer a,z,ierr,i,j

c... This internal consistency check makes sure that the correct particle module
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('empty','dsddsigma')

      
c... In the empty model, we simply set everything to zero. By using the routine
c... dsddsigma in src_models/common, we would achieve the same -- because
c... dsddgpgn returns zero in the empty module -- but this shows an explicit
c... way of bypassing/overriding that routine.

      do i=1,27
         do j=1,27
            sigij(i,j)=0.d0
         enddo
      enddo
      ierr=0
      return
      end
