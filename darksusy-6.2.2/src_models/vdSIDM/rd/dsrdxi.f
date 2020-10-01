
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
c  NB: this particular version assumes that x=1 at very high temperatures, 
c      when the DS decouples. In the vdSIDM module, this decoupling temperature
c      is set to 'infinity' in dsinit_module.
c      Further evolution after decoupling is only determined by how the 
c      effective number of d.o.f. changes in the visible and dark sector, 
c      respectively.
c
c  author: Torsten Bringmann (torsten.bringmann@fys.uio.no), 2018-06-01
c  mod: 2018-08-17 (TB, added decoupling temperature as free parameter)
c=======================================================================
      real*8 function dsrdxi(x)
      implicit none
      include 'dsvdSIDM.h'
      real*8 x
      real*8 xp, logdx, xmin, xmax, dsrdxinotab 
      
      integer kmax, dk, k
      parameter (kmax=300) 
      real*8 xdat(kmax),xidat(kmax)
      integer klow, khigh
      save xdat, xidat, klow, khigh

c to make sure that we only tabulate once per model
      integer dsidnumber
      integer idold
      data idold/-123456789/
      save idold

c...TB debug     
c      dsrdxi=1.0      
c      return

      if (idold.ne.dsidnumber()) then ! new model: first tabulate xi
        xmin = min(0.05d0, mass(kdm)/1.d4) ! DM is relativistic, but at least 1 TeV
        xmax = mass(kdm)/1.d-8             ! we tabulate down to 10 eV 
        do k = 1, kmax
          logdx = (k-1.)/(1.0d0*kmax-1.)
          xp = xmin*(xmax/xmin)**logdx
          xdat(k) = xp
          xidat(k) = dsrdxinotab(xp)
        enddo
        klow = 1
        khigh = kmax
        idold=dsidnumber()
      endif

      if (x.le.xdat(1)) then
         klow = 1
         dsrdxi = xidat(1)
      elseif (x.ge.xdat(kmax)) then
         khigh = kmax
         dsrdxi = xidat(kmax)
      else
c... find an interval [xdat(klow),xdat(khigh)] around t, 
c... starting with values from previous call
         dk=1
 100     if (xdat(khigh).lt.x) then
            khigh=khigh+dk
            dk=2*dk
            if (khigh.lt.kmax) goto 100
            khigh=kmax
         endif
 110     if (xdat(klow).gt.x) then
            klow=klow-dk
            dk=2*dk
            if (klow.gt.1) goto 110
            klow=1
         endif

c... and then narrow it down
 120     if (khigh-klow.gt.1) then
            k=(khigh+klow)/2
            if (xdat(k).lt.x) then
               klow=k
            else
               khigh=k
            endif
            goto 120
         endif
         
c... interpolate between closest tabulated points
         dsrdxi = xidat(klow)+(xidat(khigh)-xidat(klow))
     &            *(x-xdat(klow))/(xdat(khigh)-xdat(klow))

      endif

      return
      end


c_______________________________________________________________________
c Actual calculation of xi
c_______________________________________________________________________
      real*8 function dsrdxinotab(x)
      implicit none
      include 'dsvdSIDM.h'
      real*8 x
      
      real*8 gDShigh, gDS, GSMhigh, GSM, sqrtgstar ! NB: these g should actually be
                                                   ! entropy d.o.f.      
      real*8 T, Td, xiold, xi
      integer steps
      data xiold /1.0d0/
      save xiold
      
c... functions
      real*8 dsrddofDS
      
      T = mass(kdm)/x
      if (T.ge.DSTdec) then ! prior to decoupling
        dsrdxinotab=1.0d0
        return
      endif

      if (DSTdec.gt.1.d3) then
        gSMhigh = 106.75           ! all SM d.o.f.
      else  
c        call dskdgeff(DSTdec,gSMhigh) ! this would be energy d.o.f.
        call dsrddof(DSTdec,sqrtgstar,gSMhigh)
      endif
      gSM=GSMhigh
c      if (T.lt.1.d3) call dskdgeff(T,gSM) ! this would be energy d.o.f.
      if (T.lt.1.d3) call dsrddof(T,sqrtgstar,gSM)
      gDShigh = dsrddofDS(DSTdec) ! this is typically when even DM is relativistic

      steps = 0 
 10   steps = steps + 1
      Td = T*xiold      
      gDS = dsrddofDS(Td)
      xi = (gSM/gSMhigh*gDShigh/gDS)**0.33333333
      xi = xiold**0.4*xi**0.6
c      write(*,*) 'xi, xiold =', xi, xiold 
      if (abs(T*xi/Td-1.0d0).gt.0.001) then
        if (steps.lt.5000) then
          xiold=xi
          goto 10
        endif  
        xi = (xi+xiold)/2.
      endif

      dsrdxinotab = xi
            
      return
      end
