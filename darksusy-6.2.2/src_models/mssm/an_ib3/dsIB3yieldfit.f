*****************************************************************************
*** subroutine dsIB3yieldfit returns the fit parameters necessary to approximate
*** the full 3-body spectrum as a linear combination of the most extreme 
*** spectra that are possible, following the procedure described in 1510.02473.
***
***   input:   qch    -- quark channel (7-up to 12-bottom)
***            yieldk -- as in dsanyield_ch. For the moment, only gamma
***                      rays and antiprotons are implemented.
***   output:  mix      -- 1 (0) for large (small) squark mixing
***            y        -- interpolation parameter, c.f. Eq.(3,4) in 1510.02473.
***
*** Author: torsten.bringmann@fys.uio.no
*** Date: 2016-04-02
*****************************************************************************

      subroutine dsIB3yieldfit(qch,yieldk,mix,y)
      implicit none
      include 'dsmssm.h'
      include 'dsIB3com.h'


c------------------------ variables ------------------------------------

      real*8 y
      integer qch, yieldk, mix

      real*8 rvib, rheavy, rtot, gmix(2), tmp, mdm
      real*8 dx, x, xpeak, dNdxpeak, dNdx, dNdxlast
      integer spectype, i, ksq, iq
      integer khi,klo,km, dk
      data khi/nmass/
      data klo/1/
      save khi,klo

      integer dsidnumber
      integer idold
      data idold/-123456789/
      save idold

      real*8 rtrue(7:12) 
      data rtrue/6*-1.2345d50/
      save rtrue


c------------------------ functions ------------------------------------

      real*8 dsIByieldone, dsIB3yieldtab

c-----------------------------------------------------------------------

      mix=0
      y=1.0d0 ! default: pure VIB (if fit does not exist)
      mdm=mass(kn(1))

c... determine rtrue
      if (idold.ne.dsidnumber()) then
        do iq=7,12 

c very rough, but fast,  way to determine xmax 
c (improved method wrt fit --> need eventually to redo fit!?)
          xpeak = 1.d-10
          dNdxpeak = 0.0d0
          dNdxlast = 0.0d0
          x=0.8
          dx = 0.09

          if (mdm.gt.mass(iq)) then 
            do while (x.gt.1.d-4.and.x.lt.1d0.and.1.d-3.lt.abs(dx))
 		
              dNdx = mdm*dsIByieldone(mdm*x,iq,152,i)
            
              if (abs(dNdxpeak).lt.abs(dNdx)) then ! Check differential rate against maximum.
                dNdxpeak = dNdx 
                xpeak = x
              endif

		      ! Check if amplitude increasing with steps, if decreasing turn around and reduce step size.
              if (abs(dNdxlast).gt.abs(dNdx).or.(x+dx).ge.1.d0.or.(x+dx).le.1.d-4) 
c     &          dx = -dx*(abs(dx)**0.2)
     &          dx = -(abs(dx)**1.2) ! unfortunately, the fit has been performed with a wrong peak finding
                                     ! algorithm; for consistency we thus have to keep it instead of the
                                     ! commented version above...
                x = x + dx;
                dNdxlast = dNdx;
                
            enddo

            rtrue(iq) = dNdxpeak/dsIByieldone(0.001*mdm,iq,52,i) ! normalize spectrum to unity

          endif
        enddo
        idold=dsidnumber()
        
      endif


c... determine whether squark mixing is relevant
      do i=1,2
        ksq = 2*qch + i + 16
        gmix(i) = dble(gr(ksq,qch,kn(1))*(gl(ksq,qch,kn(1))))/
     &            (abs(gl(ksq,qch,kn(1)))**2 + abs(gr(ksq,qch,kn(1)))**2)
      enddo
      if (gmix(1).ge.1.d-4.or.gmix(2).ge.1.d-4) mix=1
      mix=1

c... call relevant dsIB3yieldtab, just to make sure that tables are loaded 
      tmp = dsIB3yieldtab(0.1*mdm,mdm,qch,2-mix,yieldk)
      tmp = dsIB3yieldtab(0.1*mdm,mdm,qch,4-mix,yieldk)

c... determine spectype 
      spectype = 2-mix ! this is for VIB (assume number of simulated masses for
                         ! heavy squarks are the same
      if (yieldk/100.eq.0) then    !  integrated yields
        spectype = spectype + 0
      elseif (yieldk/100.eq.1) then ! differential yields
        spectype = spectype + 4 
      endif
      if (mod(yieldk,100).eq.52) then    ! photons
        spectype = spectype + 0
      elseif (mod(yieldk,100).eq.54) then ! antiprotons
        spectype = spectype + 8 
      endif

      if (.not.fitexists(qch,spectype)) return


c... find an interval [fitmass(klo),fitmass(khi)] around mdm, 
c... starting with values from previous call
      dk=1
 200  if (fitmass(qch,spectype,khi).lt.mdm) then
        khi = khi + dk
        dk=2*dk
        if (khi.ge.nmass) then
          khi = nmass
        elseif (fitmass(qch,spectype,khi).
     &        lt.fitmass(qch,spectype,nmass)) then
          goto 200
        endif
      endif
      dk=1
 210  if (fitmass(qch,spectype,klo).gt.mdm) then
         klo = klo - dk
         dk=2*dk
         if (klo.le.1) then
           klo = 1
         else
           goto 210
         endif
      endif
c... and then narrow it down
 220  if (khi-klo.gt.1) then
         km=(khi+klo)/2
         if (fitmass(qch,spectype,km).lt.mdm) then
            klo=km
         else
            khi=km
         endif
         goto 220
      endif

c... log-interpolate of Rprime between closest tabulated points
      if (((fitmass(qch,spectype,khi)-fitmass(qch,spectype,klo))/
     &    fitmass(qch,spectype,khi)).gt.1.d-6) then
          rvib = exp(dlog(Rprime(qch,spectype,klo)) +
     &                  (dlog(Rprime(qch,spectype,khi)) -
     &                   dlog(Rprime(qch,spectype,klo))) * 
     &                  (dlog(mdm)-dlog(fitmass(qch,spectype,klo)))/
     &                   (dlog(fitmass(qch,spectype,khi))-
     &                    dlog(fitmass(qch,spectype,klo))))
      else
          rvib = Rprime(qch,spectype,khi)
      endif

c... same again, for heavy squark spectrum
      if (((fitmass(qch,spectype+2,khi)-fitmass(qch,spectype+2,klo))/
     &    fitmass(qch,spectype+2,khi)).gt.1.d-6) then
          rheavy = exp(dlog(Rprime(qch,spectype+2,klo)) +
     &                  (dlog(Rprime(qch,spectype+2,khi)) -
     &                   dlog(Rprime(qch,spectype+2,klo))) * 
     &                  (dlog(mdm)-dlog(fitmass(qch,spectype+2,klo)))/
     &                   (dlog(fitmass(qch,spectype+2,khi))-
     &                    dlog(fitmass(qch,spectype+2,klo))))
      else
          rheavy = Rprime(qch,spectype+2,khi)
      endif

      rtot = abs((abs(rtrue(qch)) - abs(rheavy)) / (abs(rvib) - abs(rheavy)))
      
      y = 0.0d0
      do i=1,3
        y = y + ci(i,qch,spectype)*rtot**ni(i,qch,spectype)
      enddo
      y = rtot*10.**y

      if (y.gt.1.05) then
c        write(*,*) 'WARNING in dsIB3yieldfit: y = ',y, qch, rtot, rtrue(qch), rvib, rheavy, spectype
        y = 1.d0
      endif  
 
      return

      end


