**********************************************************************
*** returns tabulated yield resulting from the subtracted 3-body qqg 
*** spectrum
***
***   input:   TGeV   -- Kinetic energy in GeV
***            mdm    -- dark matter mass in GeV
***            qch    -- quark channel (7-up to 12-bottom)
***            threebodytype -- 1: VIB, maximal mixing
***                             2: VIB, minimal mixing
***                             3: heavy squark limit, maximal mixing
***                             4: heavy squark limit, minimal mixing
***                                (for definitions, see arxiv: 1510.02473)
***            yieldk -- as in dsanyield_ch. For the moment, only gamma
***                      rays and antiprotons are implemented.
***
*** the units are (annihilation into IBch)**-1
*** for the differential yields, the units are the same times gev**-1.
***
*** author: torsten.bringmann@fys.uio.no, 2015-08-07
***********************************************************************

      real*8 function dsIB3yieldtab(TGev,mdm,qch,threebodytype,yieldk)

      implicit none
      real*8 TGeV, x, mdm, resultm1, resultm2
      integer qch, spectype, threebodytype, yieldk, k,dk
      integer ehitmp, elotmp, mhitmp,mlotmp

      include 'dsIB3com.h'

      integer ntot
      parameter (ntot=6*2*ntype) ! 6 quark types
      logical initialized(7:18,ntype) 
      data initialized / ntot*.false. /
      save initialized

      dsIB3yieldtab=0.0d0
      x=TGeV

c... check that qch is in the allowed range
      if (qch.ge.13.and.qch.le.18) then
c        write(*,*) 'WARNING from dsIB3yieldtab: Assessing 2-body yields',
c     &             ' rather than 3-body yields!'        
      elseif(qch.lt.7.or.qch.gt.12) then
        write(*,*) 'ERROR:  dsIB3yieldtab called with unsupported qch = ',qch
        stop
        return        
      endif

c... combine threebodytype and yieldk to spectype
      if (threebodytype.ge.1.and.threebodytype.le.4) then
        spectype = threebodytype
      else
        goto 1000
      endif
      if (yieldk/100.eq.0) then    !  integrated yields
        spectype = spectype + 0
      elseif (yieldk/100.eq.1) then ! differential yields
        spectype = spectype + 4 
      else
        goto 1000
      endif
      if (mod(yieldk,100).eq.52) then    ! photons
        spectype = spectype + 0
        x=x/mdm                      ! for gamma (only!), we use dimensionless 
                                     ! units for the energy in tabulated data
      elseif (mod(yieldk,100).eq.54) then ! antiprotons
        spectype = spectype + 8 
      else
        goto 1000
      endif

      if (.not.initialized(qch,spectype)) then
        call dsIB3importdata(qch,spectype) 
        initialized(qch,spectype) = .true.
      endif

 
c... this is introduced just to make the following code a bit simpler to read
      ehitmp = ehi(qch,spectype)
      elotmp = elo(qch,spectype)
      mhitmp = mhi(qch,spectype)
      mlotmp = mlo(qch,spectype)
 
c... find an interval [xdat(klo),xdat(khi)] around TGev, 
c... starting with values from previous call
         dk=1
 100     if (xdat(qch,spectype,ehitmp).lt.x) then
            ehitmp = ehitmp + dk
            dk=2*dk
            if (ehitmp.ge.nenergy) then
              ehitmp = nenergy
            elseif (xdat(qch,spectype,ehitmp).
     &            lt.xdat(qch,spectype,nenergy)) then
              goto 100
            endif
         endif
         dk=1
 110     if (xdat(qch,spectype,elotmp).gt.x) then
            elotmp = elotmp - dk
            dk=2*dk
            if (elotmp.le.1) then
              elotmp = 1
            else
              goto 110
            endif
         endif

c... and then narrow it down
 120     if (ehitmp-elotmp.gt.1) then
            k=(ehitmp+elotmp)/2
            if (xdat(qch,spectype,k).lt.x) then
               elotmp=k
            else
               ehitmp=k
            endif
            goto 120
         endif


c... now do the same for the limiting masses
         dk=1
 200     if (fitmass(qch,spectype,mhitmp).lt.mdm) then
            mhitmp = mhitmp + dk
            dk=2*dk
            if (mhitmp.ge.nmass) then
              mhitmp = nmass
            elseif (fitmass(qch,spectype,mhitmp).
     &            lt.fitmass(qch,spectype,nmass)) then
              goto 200
            endif
         endif
         dk=1
 210     if (fitmass(qch,spectype,mlotmp).gt.mdm) then
            mlotmp = mlotmp - dk
            dk=2*dk
            if (mlotmp.le.1) then
              mlotmp = 1
            else
              goto 210
            endif
         endif
 220     if (mhitmp-mlotmp.gt.1) then
            k=(mhitmp+mlotmp)/2
            if (fitmass(qch,spectype,k).lt.mdm) then
               mlotmp=k
            else
               mhitmp=k
            endif
            goto 220
         endif

c... store limiting tags in common block to speed up point identification for next call
      ehi(qch,spectype) = ehitmp
      elo(qch,spectype) = elotmp
      mhi(qch,spectype) = mhitmp
      mlo(qch,spectype) = mlotmp

c... log-interpolate of dNdE between closest tabulated points, for both masses
      if (((xdat(qch,spectype,ehitmp)-xdat(qch,spectype,elotmp))/
     &    xdat(qch,spectype,ehitmp)).gt.1.d-6) then
          resultm1 = exp(dlog(dNdEdat(qch,spectype,mlotmp,elotmp)) +
     &                  (dlog(dNdEdat(qch,spectype,mlotmp,ehitmp)) -
     &                   dlog(dNdEdat(qch,spectype,mlotmp,elotmp))) * 
     &                  (dlog(x)-dlog(xdat(qch,spectype,elotmp)))/
     &                   (dlog(xdat(qch,spectype,ehitmp))-
     &                    dlog(xdat(qch,spectype,elotmp))))
          resultm2 = exp(dlog(dNdEdat(qch,spectype,mhitmp,elotmp)) +
     &                  (dlog(dNdEdat(qch,spectype,mhitmp,ehitmp)) -
     &                   dlog(dNdEdat(qch,spectype,mhitmp,elotmp))) * 
     &                  (dlog(x)-dlog(xdat(qch,spectype,elotmp)))/
     &                   (dlog(xdat(qch,spectype,ehitmp))-
     &                    dlog(xdat(qch,spectype,elotmp))))
      else
          resultm1 = dNdEdat(qch,spectype,mlotmp,ehitmp)
          resultm2 = dNdEdat(qch,spectype,mhitmp,ehitmp)
      endif
  
c... re-introduce dimensions for differential photon spectrum!
      if (((spectype-1)/8).eq.0.and.((spectype-1)/4).eq.1) then
        resultm1=resultm1/mdm
        resultm2=resultm2/mdm
      endif

c... now log-interpolate between the masses

      if ((abs(resultm2-resultm1)/resultm2).gt.1.d-6) then
        dsIB3yieldtab = exp(dlog(resultm1) + (dlog(resultm2) - dlog(resultm1)) * 
     &       (dlog(mdm)-dlog(fitmass(qch,spectype,mlotmp)))/
     &       (dlog(fitmass(qch,spectype,mhitmp))-
     &        dlog(fitmass(qch,spectype,mlotmp))))       
      else
          dsIB3yieldtab = resultm1
      endif

c      write(*,*) 'masses = ',fitmass(qch,spectype,mlotmp),fitmass(qch,spectype,mhitmp)
c      write(*,*)'E= ',xdat(qch,spectype,elotmp), xdat(qch,spectype,ehitmp)
c      write(*,*) 'dNdEdat = ', dNdEdat(qch,spectype,mlotmp,elotmp),
c     &                         dNdEdat(qch,spectype,mlotmp,ehitmp),
c     &                         dNdEdat(qch,spectype,mhitmp,elotmp),
c     &                         dNdEdat(qch,spectype,mhitmp,ehitmp)
c      write(*,*) 'resultm1,2 = ', resultm1, resultm2

      return
 
 1000 write(*,*) 'ERROR in dNdE! Requested combination of threebodytype and ',
     &           'yieldk is not implemented: ',threebodytype,yieldk
      stop
      return 

      end
