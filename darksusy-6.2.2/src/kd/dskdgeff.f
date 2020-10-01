**********************************************************************
*** returns effective degrees of freedom during kinetic decoupling
***
***   input:   T    -- Temperature in GeV
***   output:  geff -- rel. energy d.o.f [\rho = geff*(pi^2/30)*T^4 ]
***
*** author: Torsten Bringmann (torsten.bringmann@fys.uio.no), 2017-05-12
***         (modified version of Paolo Gondolos dsrddof)
***********************************************************************

      subroutine dskdgeff(T,geff)

      implicit none
      real*8 T,geff
      integer k,dk

      include 'dsrdcom.h'
      include 'dskdcom.h'


      logical initialized
      data initialized / .false. /
      save initialized

c... NB: this assumes that tgev(1) [tgev(nf)] is the highest [lowest] tabulated T !!!

      if (.not.initialized) then
        kklo = 1
        kkhi = nf
        initialized = .true.
      endif

      if (t.ge.tgev(1)) then
         geff = fe(1)
         kklo = 1
      elseif (t.le.tgev(nf)) then
         geff = fe(nf)
         kkhi = nf
      else

c... This check should no longer be necessary with above
c... initialization check unless kklo or kkhigh are modified in
c... another routine
c        if (kklo.lt.1.or.kklo.gt.nf) kklo=1
c        if (kkhi.gt.nf.or.kkhi.lt.1) kkhi=nf

c... find an interval [tgev(kklo),tgev(kkhi)] around t, 
c... starting with values from previous call
         dk=1
 100     if (tgev(kkhi).gt.t) then
            kkhi=kkhi+dk
            dk=2*dk
            if (kkhi.lt.nf) goto 100
            kkhi=nf
         endif
         dk=1
 110     if (tgev(kklo).lt.t) then
            kklo=kklo-dk
            dk=2*dk
            if (kklo.gt.1) goto 110
            kklo=1
         endif

c... and then narrow it down
 120     if (kkhi-kklo.gt.1) then
            k=(kkhi+kklo)/2
            if (tgev(k).gt.t) then
               kklo=k
            else
               kkhi=k
            endif
            goto 120
         endif
         
c... interpolate between closest tabulated points
         geff = fe(kklo)+(fe(kkhi)-fe(kklo))
     &          *(t-tgev(kklo))/(tgev(kkhi)-tgev(kklo))

      endif

      return
      end
