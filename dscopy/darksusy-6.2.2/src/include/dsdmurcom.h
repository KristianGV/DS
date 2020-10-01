*         -*- mode: fortran -*-

c
c codes for dm user-replaceable functions:
c
c type of coordinates used
c      
      integer tysph  ! type: spherically symmetric
      integer tyaxi  ! type: axisymmetric
      integer tytri  ! type: triaxial
      parameter(tysph=1,tyaxi=2,tytri=3)
c      
c for each type specify entries for coordinate vector
c
      integer irsph
      parameter(irsph=1)
      integer iraxi,izaxi
      parameter(iraxi=1,izaxi=2)
      integer ix,iy,iz
      parameter(ix=1,iy=2,iz=3)
c
c which kind of function linked in set of functions which need to be
c consistently defined:      
c      
      integer krho        ! density profile rho [GeV cm^-3]
      integer ksoupow1    ! source function with power = 1, i.e. rho +
                          ! eventual substructure effects [GeV cm^-3]
      integer ksoupow2    ! source function with power = 2, i.e. rho**2 +
                          ! eventual substructure effects [GeV^2 cm^-6]
      integer kmass       ! mass (volume integral of rho) [M_sun] 
      integer kdrhodr     ! first derivative of rho w.r.t. radius
                          ! [GeV cm^-3 kpc^-1]
      integer kd2rhodr2   ! second derivative of rho w.r.t. radius
                          ! [GeV cm^-3 kpc^-2]
      integer kphi        ! potential energy phi [c^2]
      integer kdphidr     ! derivative of phi w.r.t. radius [c^2 kpc^-1]
      integer kd2phidr2   ! second derivative of phi w.r.t. radius
                          ! [c^2 kpc^-2]
      parameter(krho=1,ksoupow1=2,ksoupow2=3,kmass=4,kdrhodr=5,
     &  kd2rhodr2=5,kphi=6,kdphidr=7,kd2phidr2=8)      

      
