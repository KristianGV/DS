      subroutine dshmucmhset(mucmh, rucmh, rcore, dist)
****************************************************************
*** subroutine dshmucmhset:                                  ***
*** Perform some initialisation required when the UCMH mass, ***
*** radii or distance have been set/reset.                   ***
***                                                          ***
*** author: pat scott (patscott@physics.mcgill.ca)           ***
*** date: 2010-10-17                                         ***
****************************************************************
      implicit none
      include 'dsucmh.h'

      real*8 mucmh, rucmh, rcore, dist

      Mh_ucmh = mucmh
      Rh_ucmh = rucmh
      rhcut = rcore
      r_0 = dist

      end
