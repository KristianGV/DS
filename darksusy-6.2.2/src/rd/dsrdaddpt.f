      subroutine dsrdaddpt(wrate,pres,deltap,deltapminr)
c_______________________________________________________________________
c  add a point in rdrate table
c  input:
c    wrate - invariant annihilation rate (real, external)
c    pres - momentum of the point to add
c    deltap - scaling factor used in dsrdtab
c    deltapminr - minimum delta p between existing point and new point
c                 (in units of WIMP mass) for it to be added      
c    pmax - maximum p used in dsrdtab (from common block)
c  common:
c    'dsrdcom.h' - included common blocks
c  used by dsrdtab
c  author: joakim edsjo (edsjo@fysik.su.se)
c  modified: 01-01-31 paolo gondolo (paolo@mamma-mia.phys.cwru.edu)
c  modified: 15-12-08 Joakim Edsjo, to allow for wrate=0      
c=======================================================================
      implicit none
      include 'dsrdcom.h'
      real*8 wrate,dsrdlny
      external wrate
      integer i,j
      real*8 pres,deltap
      real*8 deltapminr
      logical pexist
c-----------------------------------------------------------------------

c...Check if we have exceeded time constraint
      call CPU_TIME(rdt_end)
      if ((rdt_end-rdt_start).gt.rdt_max) then
         rderr=ibset(rderr,10)
         write(*,*) 'DS ERROR in dsrdaddpt: time limit exceeded.'
         write(*,*) '  Time spent: ',rdt_end-rdt_start,' seconds'
         write(*,*) '  Maximum allowed time: ',rdt_max,' seconds'
         write(*,*) '  omega calculation stopping.'
         return
      endif

      if (pres.le.0.0d0) return
c...check if p already exists [serious bug corrected pg 01-01-31]
c...Floating Point exception now working 2008-09-09 JE
      pexist=.false.
      do i=1,nr
         if (pp(i).eq.0.d0.and.pres.eq.0.d0) then
            pexist=.true.
         else if
c     &     (abs((pp(i)-pres/deltap)/max(pp(i),1.d-9)).lt.1.d-9) then ! org
     &     (abs((pp(i)-pres/deltap)/max(pp(i),deltapminr))
     &        .lt.deltapminr) then ! new
c     &     (abs((pp(i)-pres/deltap)/max(pp(i),dpminr)).lt.dpminr) then ! alt2
c     &           (abs(pp(i)-pres/deltap).lt.deltapminr) then ! JE alternative
            pexist=.true.
         endif
      enddo
      if (pexist) return

      if (pres.gt.pmax) then
        if (rdprt.gt.0) write(*,*)
     &    'dsrdaddpt: pmax raised from ',pmax,' to ',pres
        pmax=pres
        if (nr+1.ge.nrmax) then
          write (rdluerr,*) 
     &      'error in dsrdaddpt: array capacity exceeded'
          write (rdluerr,*) '  for model ',rdtag
          write(rdluerr,*) 'Omega calculation stopping.'
          rderr=ibset(rderr,0)
          return
        endif
        nr=nr+1
        pp(nr)=pres/deltap
        yy(nr)=dsrdlny(pres,wrate) ! JE Correction 2015-12-08
        indx(nr)=nr
      else
        if (nr+1.ge.nrmax) then
          write (rdluerr,*) 
     &      'error in dsrdaddpt: array capacity exceeded'
          write (rdluerr,*) '  for model ',rdtag
          write(rdluerr,*) 'Omega calculation stopping.'
          rderr=ibset(rderr,0)
          return
        endif
        nr=nr+1
        pp(nr)=pres/deltap
        yy(nr)=dsrdlny(pres,wrate) ! JE Correction 2015-12-08
        do i=1,nr-2
          if (pp(indx(i)).lt.pp(nr).and.pp(indx(i+1)).gt.pp(nr)) then
            do j=nr-1,i+1,-1
              indx(j+1)=indx(j)
            enddo
            indx(i+1)=nr
          endif
        enddo
      endif

      return
      end


