***********************************************************************
*** This function returns the WIMP mass                             
***
*** author: torsten.bringmann@fys.uio.no, 2014-05-09
***********************************************************************

      real*8 function dsmwimp()
      implicit none

      include 'dsparticles.h'

      dsmwimp=mass(kdm)

      return
      end
