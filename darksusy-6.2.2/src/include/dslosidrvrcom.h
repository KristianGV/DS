*         -*- mode: fortran -*-

c
c  options for a valid call to line of sight integration routines, as linked
c  within the driver dslosidriver 
c
      integer ilosigasph  ! link for gamma-ray flux spherical source
      integer ilosisysph  ! link for synchrotron flux spherical source
      integer ilosiICsph  ! link for Inverse Compton flux spherical source
      integer ilosibrsph  ! link for Bremstrahlung flux spherical source
      integer ilosinsph   ! maximum number of links 
      parameter(ilosigasph=1,ilosisysph=2,ilosiICsph=3,ilosibrsph=4,
     &  ilosinsph=4)

      integer ilosiabssft  ! how much you shift the previous indices
                           ! in case (self-)absorption is included
      parameter(ilosiabssft=100)

      integer ilositautot  ! how much you shift the previous indices
                           ! in case you link dslosidriver to compute or
                           ! or set the total optical depth
      parameter(ilositautot=200)

      integer ilositauset  ! how much you shift the previous indices
                           ! in case you link dslosidriver to set the
                           ! function computing absoption along l.o.s.i.
      parameter(ilositauset=300)

      integer ilositaufun  ! how much you shift the previous indices
                           ! in case you link dslosidriver to link the
                           ! function computing absoption along l.o.s.i.
      parameter(ilositaufun=400)
      
      integer ilosiabscomp      ! absorption computed in dslosidriver
      parameter(ilosiabscomp=-1)
