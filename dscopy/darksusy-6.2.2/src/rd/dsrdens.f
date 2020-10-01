      subroutine dsrdens(wrate,oh2,xf,fast,ierr,iwar)
c_______________________________________________________________________
c  present density in units of the critical density times the
c    hubble constant squared. Called from dsrdomega.
c  input:
c    wrate - invariant annihilation rate (real, external)
c    fast - see meaning of options in dsrdomega header
c  output:
c    oh2   - relic density parameter times h**2 (real*8)
c    xf    - m_WIMP / T_f where T_f is the photon temperature at freeze-out
c    ierr  - error code (integer)
c      bit 0 (1) = array capacity exceeded. increase nrmax in dsrdcom.h
c          1 (2) = a zero vector is given to dsrdnormlz.
c          2 (4) = step size underflow in dsrdeqn
c          3 (8) = stepsize smaller than minimum hmin in dsrdeqn
c          4 (16) = too many steps in dsrdeqn
c          5 (32) = step size underflow in dsrdqrkck
c          6 (64) = step size smaller than miminum in dsrdqrkck
c          7 (128) = too many steps in dsrdqrkck
c          8 (256) = gpindp integration failed in dsrdthav
c          9 (512) = threshold array too small. increase tharsi in dsrdcom.h
c         10 (1024) = calculation took longer than max time in rdt_max,
c                     calculation aborted      
c    iwar  - warning code (integer)
c      bit 0 (1) = a difference of >5waccd in the ratio of w_spline
c                  and w_linear is obtained due to delta_p<dpmin.
c          1 (2) = a difference of >10waccd in the ratio of w_spline
c                  and w_linear is obtained due to delta_p<dpmin.
c          2 (4) = a difference of >15waccd in the ratio of w_spline
c                  and w_linear is obtained due to delta_p<dpmin.
c          3 (8) = wimp too heavy, d.o.f. table needs to be
c                  extended to higher temperatures. now the solution
c                  is started at a higher x than xinit (=2).
c          4 (16) = spline interpolated value too high (overflow) during
c                   check of interpolation accuracty (dsrdwintpch)
c  common:
c    'dsrdcom.h' - included common blocks
c  uses dsrdtab, dsrdeqn.
c  authors: Paolo Gondolo (gondolo@lpthe.jussieu.fr) 1994-1996 and
c           Joakim Edsjo (edsjo@fysik.su.se) 30-april-98
c  modified:
c     2013-10-04 paolo gondolo: totally decouple from mssm
c     2018-04-26 torsten bringmann: version without tabulation
c=======================================================================
      implicit none
      include 'dsrdcom.h'
      include 'dsio.h'
      integer nfcn,ierr,iwar,fast
      real*8 wrate,dsrdwintp,oh2,xstart,xend,yend,xf
      external wrate,dsrdwintp
      integer k,i
      real*8 tstart
c-----------------------------------------------------------------------

      call CPU_TIME(rdt_start) ! save start time of calculation

      rderr=0
      rdwar=0

      nrd=0                     ! to have dsrdwx initialized

c...make sure that d.o.f. tables are read in
      if (rdinit.ne.1234) then
         write(*,*) 'DS ERROR: dsrdinit (called from dsinit) is',
     &        ' not called prior to dsrdomega.'
         stop
      endif
      
      if (rderr.ne.0) goto 999

c-------------------------------------------------------------- initial x

      xstart=max(xinit,1.0001d0*mco(1)/tgev(1))  ! je corr 97-03-14
      if (xstart.ne.xinit) rdwar=ibset(rdwar,3)

c--------------------------------------------------- locate tstart in dof

      tstart=mco(1)/xstart
      khi=nf
      klo=1
  100 if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if (tgev(k).lt.tstart) then
          khi=k
        else
          klo=k
        endif
        goto 100
      endif

c--------- setup resonance flags
      do i=1,nres
         resinc(i)=.true.
         resfit(i)=.false.
      enddo
      
c------------------------------------------- tabulate the invariant rate

      if (fast.eq.20) then          ! new method, only do resonance optimization
         call dsrdbw_setup(wrate,xstart)
      elseif (fast.ne.99) then ! fast.ne.20, or fast.ne.99 tabulate first
         call dsrdtab(wrate,xstart,fast)
         if (rderr.ne.0) goto 999
      endif

c------------------------------------------ determine integration limits

      call dsrdthlim

c---------------------------------------------------------- call dsrdeqn

      if (fast.eq.20) then ! don't use tabulated result
        call dsrdeqn(wrate,xstart,xend,yend,xf,nfcn)
      elseif (fast.eq.99) then ! don't use tabulated result
        call dsrdeqn(wrate,xstart,xend,yend,xf,nfcn)
      else
        call dsrdeqn(dsrdwintp,xstart,xend,yend,xf,nfcn)
      endif  

c----------------------------------------- normalize to critical density

 999  if (rderr.eq.0) then
        oh2=0.70365d8*fh(nf)*mco(1)*yend  ! je 970404, t_0=2.726 k
      else
        oh2=0.0d0     ! something went wrong, changed to 0.0 020816, je
        xf=-1.0d0
      endif
      ierr=rderr
      iwar=rdwar
c      write (*,*) 'PG-DEBUG dsrdens: oh2,yend,mco(1),fh(nf)',
c     &     oh2,yend,mco(1),fh(nf)

      call CPU_TIME(rdt_end)

c...Write out for possible tests or debugging of fast=20 option
c      if (fast.eq.20) then
c      open(unit=55,file=idtag(1:12)//'-dsrdwx.dat',status='unknown',
c     &  form='formatted')
c      do k=1,nrd 
c         write(55,*) k,ppp(k),rdwx(k)
c      enddo
c      close(55)

c      dpminr=dpminr*10.d0 ! to avoid new tabulations
c      open(unit=55,file=idtag(1:12)//'-dsrdwx-det.dat',status='unknown',
c     &  form='formatted')
c      do k=0,3000
c         ptmp=dble(k)/3000.d0*mco(1)*3.d0
c         write(55,*) k,ptmp,dsrdwx(ptmp)
c      enddo
c      close(55)
c      dpminr=dpminr/10.d0
c      endif
      
      

      end
