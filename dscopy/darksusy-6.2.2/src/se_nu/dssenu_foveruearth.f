***********************************************************************
*** input: velocity relative to Earth [ km s^-1 ]
*** output: f(u) / u [ cm^-3 (cm/s)^(-2) ]
*** Date: 2004-01-28
***********************************************************************

      real*8 function dssenu_foveruearth(u)
      implicit none

      include 'dssecom.h'
      real*8 u,dshmuDFearth

      dssenu_foveruearth=dshmuDFearth(u) ! use user-defined halo profile
     &  *1.0d-10                      ! (km/s)^(-2) -> (cm/s)^(-2)
     &  *(serho/semx)                  ! normalize to local number density

      return
      end
