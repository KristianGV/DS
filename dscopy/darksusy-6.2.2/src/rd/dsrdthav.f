      real*8 function dsrdthav(x,wrate)
c_______________________________________________________________________
c  the thermal average of the effective annihilation cross section.
c  input:
c    x - mass/temperature (real)
c        NB: this temperature does not have to be the photon temperature!
c    wrate - invariant annihilation rate (real)
c  output:
c    dsrdthav - thermal averged cross section
c  common:
c    'dsrdcom.h' - included common blocks
c  uses qrkck or dgadap.
c  called by dsrdrhs
c  author: joakim edsjo (edsjo@fysik.su.se) 98-05-01
c  modified 2018-04-26 (torsten bringmann) : added external wrate as argument
c                                            to dsrdfuncs (and gpindp2 & dgadap2),
c                                            and updated integration to dqagp
c=======================================================================
      implicit none
      include 'dsrdcom.h'
      real*8 wavs,dsrdfuncs
      real*8 x,wrate,dsrdfunc,dsai,dsbi,
     &  epsin,epsout,gpindp2
      integer i,iop,k
      external dsrdfunc,wrate,dsrdfuncs,gpindp2
      real*8 wav,xmin,wtmp,utop,tmp,ptmp
      logical dsisnan

c required for dqagp2
      integer npts2, neval, ier, leniw,lenw,last,iwork(10000)
      real*8 points(10000), epsabs,epsrel,result,abserr,work(10000)


c-----------------------------------------------------------------------
c..added by je 97-01-17 for gadap integration
      real*8 rdx
      common /gadint2/ rdx

c-------------- thermally averaged cross section times relative velocity

c     integrate as far as umax in dsrdtab allows  ! je 97-01-07
      utop=umax*0.9999d0    ! must be <= umax, or rather the maximal u
                            ! corresponding to pmax set in dsrdtab

c...define xmin. xmin >= xstart in dsrdens required.
      xmin=max(xinit,1.0001d0*mco(1)/tgev(1))  ! je corr 97-03-14
      if (xmin.ne.xinit) rdwar=ibset(rdwar,3)

      nlo=1
      nhi=2
      

      if (thavint.eq.1) then ! gadap integration (thavint=1)
        rdx=x
        wav=0.0d0
        do i=1,nlim
          dsai=sqrt(2.0d0*x*(sqrt(1.0d0+(plow(i))**2/
     &          mco(1)**2)-1.0d0))
c...      The following if statement adds a linear interpolation over threshold
c...      to avoid missing the threshold when integration split at it.          
          if (i.gt.1.and.dsai.gt.dsbi) then
             wavs=0.5d0*(dsrdfuncs(wrate,dsai)+dsrdfuncs(wrate,dsbi))*(dsai-dsbi)
             wavs=wavs/mco(1)**2
             wav=wav+wavs/1.0d15
          endif
          dsbi=min(sqrt(2.0d0*x*(sqrt(1.0d0+(phigh(i))**2
     &      /mco(1)**2)-1.0d0)),utop)

          call dgadap2_thav(dsai,dsbi,dsrdfuncs,wrate,0.0001d0,wavs)
          
c TB debug: changing integration routine as below results in identical results
c note that accuracy had to be decreased to 1.0d-3
          
c          npts2=2
c          wavs=dsrdfuncs(wrate,dsai) ! why is this necessary ???
c          wavs=dsrdfuncs(wrate,dsbi) ! why is this necessary ???
          
c          leniw=100
c          lenw=200          
c         call dqagp2(dsrdfuncs,wrate,dsai,dsbi,
c     &     npts2,points,1.0d0, 1.0d-3, ! absolute accuracy not taken into account so far
c     &     result,abserr, neval,ier,leniw,lenw,last,iwork,work) ! output

cc           write(*,*) 'COMP: ', wavs, result

c          if (ier.ne.0) then
c            write(*,*) result,abserr, neval,ier,leniw,lenw,last !,iwork,work
c            stop
c          endif
     
c          wavs = result
          
          
c...divide out mco(1)^2 introduced in dsrdfunc
          wavs=wavs/mco(1)**2
          wav=wav+wavs/1.0d15
          if (dsbi.ge.0.9999d0*utop) goto 130
        enddo

      elseif (thavint.eq.2) then ! runge-kutta integration (thavint=2)

        wav=0.0d0
        do i=1,nlim
          dsai=sqrt(2.0d0*x*(sqrt(1.0d0+(plow(i))**2/
     &        mco(1)**2)-1.0d0))
          if (i.gt.1.and.dsai.gt.dsbi) then
             wtmp=0.5d0*(dsrdfuncs(wrate,dsai)+dsrdfuncs(wrate,dsbi))*(dsai-dsbi)
             wtmp=wtmp/mco(1)**2
             wav=wav+wtmp/1.0d15
          endif
          dsbi=min(sqrt(2.0d0*x*(sqrt(1.0d0+(phigh(i))**2
     &      /mco(1)**2)-1.0d0)),utop)
          call dsrdqrkck(dsrdfunc,x,wrate,dsai,dsbi,wtmp)
c...divide out mco(1)^2 introduced in dsrdfunc
          wtmp=wtmp/mco(1)**2
          wav=wav+wtmp
          if (dsbi.ge.0.9999d0*utop) goto 130
        enddo

      elseif (thavint.eq.3) then  ! gpindp (thavint=3)

        rdx=x
        wav=0.0d0
        epsin=0.001d0
        epsout=0.001d0
        iop=1
        do i=1,nlim
          dsai=sqrt(2.0d0*x*(sqrt(1.0d0+(plow(i))**2/
     &        mco(1)**2)-1.0d0))
          if (i.gt.1.and.dsai.gt.dsbi) then
             wavs=0.5d0*(dsrdfuncs(wrate,dsai)+dsrdfuncs(wrate,dsbi))*(dsai-dsbi)
             wavs=wavs/mco(1)**2
             wav=wav+wavs/1.0d15
          endif
          dsbi=min(sqrt(2.0d0*x*(sqrt(1.0d0+(phigh(i))**2
     &      /mco(1)**2)-1.0d0)),utop)
          wavs = gpindp2(dsai,dsbi,epsin,epsout,dsrdfuncs,wrate,iop)
c...divide out mco(1)^2 introduced in dsrdfunc
          wavs=wavs/mco(1)**2

          wav=wav+wavs/1.0d15
          if (iop.eq.0) then
            rderr=ibset(rderr,8)
            write(*,*) 'error in dsrdrhs: gpindp integration failed'
            write(*,*) '  for model ',rdtag
            goto 130
          endif
          if (dsbi.ge.0.9999d0*utop) goto 130
        enddo

      elseif (thavint.eq.4) then ! dqagp integration
         rdx=x
         wav=0.0d0
         do k=1,2               ! do split into two regions just to be safe
            if (k.eq.1) then
               dsai=sqrt(2.0d0*x*(sqrt(1.0d0+(plow(1))**2/
     &              mco(1)**2)-1.0d0))
               dsbi=min(sqrt(2.0d0*x*(sqrt(1.0d0+(mco(1))**2
     &              /mco(1)**2)-1.0d0)),utop)
            else
               dsai=dsbi
               dsbi=min(sqrt(2.0d0*x*(sqrt(1.0d0+(phigh(nlim))**2
     &              /mco(1)**2)-1.0d0)),utop)
            endif

         wavs=dsrdfuncs(wrate,dsai) ! to be sure we have endpoints
         wavs=dsrdfuncs(wrate,dsbi) ! to be sure we have endpoints

         npts2=2

c TB: to add N intermediate points, increase npts2 by N and define location
c as points(1),...,points(N)

c...JE add break points, only add at resonances we know and care about
c      do i=1,nres
c         if (resinc(i)) then
c            ptmp=sqrt(rgev(i)**2/4-mco(1)**2)
c            tmp=sqrt(2.0d0*x*(sqrt(1.0d0+(ptmp)**2/
c     &           mco(1)**2)-1.0d0))
c            if (tmp.gt.dsai.and.tmp.lt.dsbi) then
c               npts2=npts2+1
c               points(npts2-2)=tmp
c            endif
c         endif
c      enddo

c...Add integration limits as breakpoints here
         do i=1,nlim
            ptmp=phigh(i)
            tmp=sqrt(2.0d0*x*(sqrt(1.0d0+(ptmp)**2/
     &           mco(1)**2)-1.0d0))
            if (tmp.gt.dsai.and.tmp.lt.dsbi) then
               npts2=npts2+1
               points(npts2-2)=tmp
            endif
         enddo

         leniw=3000
         lenw=6000

         call dqagp2(dsrdfuncs,wrate,dsai,dsbi,
     &        npts2,points,1.0d0,1.0d-3, ! absolute accuracy not taken into account so far
     &        result,abserr, neval,ier,leniw,lenw,last,iwork,work) ! output

         if (ier.ne.0.and.ier.ne.2) then
            write(*,*) 'dqagp2 debug...'
            write(*,*) npts2
            write(*,*) dsai,dsbi
            do i=1,npts2-2
               write(*,*) i,points(i)
            enddo
            write(*,*) result,abserr, neval,ier,leniw,lenw,last !,iwork,work
            stop
         endif

         wavs=result

c...divide out mco(1)^2 introduced in dsrdfunc
         wavs=wavs/mco(1)**2
         wav=wav+wavs/1.0d15
         enddo
         
      else                      ! invalid thavint
         write(*,*)
     &     'DS ERROR in dsrdthav. Invalid integration method choice, ',
     &     'thavint =',thavint
         stop
      endif

 130  dsrdthav=wav

      return
      end








