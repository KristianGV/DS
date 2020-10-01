
c_______________________________________________________________________
c  Function dsrdxi returns the temperature ratio for the temperature of
c  the heat bath that keeps DM in thermal equilibirum, Tdark, and the
c  photon temperature, Tphoton.
c
c  input:
c    x - mass/Tphoton
c
c  output:
c    xi = Tdark/Tphoton
c
c  NB: this particular version resides in src/, and simply returns
c      a constant 1.0d0. A particle module can replace this function
c      to implement a constant ratio different from 1, or a more
c      complicated behaviour. In this case, all routines in src/rd/ and
c      arc/kd/ will assume that the temperature of the heat bath 
c      in contact with DM does *not* coincide with the photon temperature.
c
c  author: Torsten Bringmann (torsten.bringmann@fys.uio.no), 2018-05-28
c=======================================================================
      real*8 function dsrdxi(x)
      implicit none
      real*8 x

      if (x.le.0.0d0) then
        write(*,*) 'ERROR in dsrdxi: called with x = ',x
        stop
      endif

      dsrdxi = 1.0d0
      
      return
      end
