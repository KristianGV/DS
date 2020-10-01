
c_______________________________________________________________________
c  Function dsrddeltaneff returns the additional number of effective
c  relativistic degrres of freedom, in units of the contribution by
c  one neutrino species.
c
c  input:
c    T - photon temperature
c
c  NB: This functions requires the functiom dsrddofDS that returns the
c      relativistic degrees of freedom in the dark sector (typically
c      provided as an interface function by the particle module, see
c      e.g. src_models/vdSIDM/rd/) 
c
c  author: Torsten Bringmann (torsten.bringmann@fys.uio.no), 2018-05-31
c=======================================================================
      real*8 function dsrddeltaneff(T)
      implicit none
      real*8 T
  
      real*8 Tnu, Tnutmp, Td, x, sqrtgstar,heff, heffsav, heffep

c functions
      real*8 dsrddofDS, dsrdxi, dsmwimp

      x   = dsmwimp()/T
      Td  = T*dsrdxi(x)
      Tnu = T

      if (T.lt. 2.d-3) then
        Tnu = 0.71379*T  ! (4./11.)**0.3333*T
        if (T.gt.3.d-5) then ! improve estimate for low-T result
          heffep = 5.5d0 ! fotons and e^+-
          call dsrddof(T,sqrtgstar,heffsav)
  10      Tnutmp = Tnu
          heff = heffsav - 3.*7./4.*(Tnu/T)**3 !subtract neutrino contribution
c          write(*,*) 'heff: ',heffsav, heff
          Tnu = (heff/heffep)**0.3333*T
          if (abs(Tnutmp/Tnu-1.).gt.0.01) goto 10
        endif  
c        write(*,*) 'Tnu/T = ', Tnu/T
      endif 
 
      dsrddeltaneff = 4/7.*dsrddofDS(Td)*(TD/Tnu)**4
      
      return
      end
