      subroutine dsrdbw_setup(wrate,xmin)
**********************************************************************
*** This routine will go through all supplied resonances and see if
***   a) they are important enough to actually include in the calculation
***   b) if they can be fit with a Breit-Wigner form which can be used
***      instead of calling the invariant rate routine. If this is the
***      case, a fit with a momentum-dependent correction factor will be
***      applied and used and it will be determined how far away from the
***      resonance that this is useful
***
*** The routine takes input from and write outpout to common blocks
*** in dsrdcom.h. The routine dsrdwx will use this info when calculating
*** the effective invariant annihilation rate. Currently this is used
*** only for the fast=20 option.      
***
*** Author: Joakim Edsjo, edsjo@fysik.su.se
*** Date: 2018-10-26      
**********************************************************************      

      implicit none
      include 'dsrdcom.h'
      include 'dsio.h'
      real*8 xmin
      real*8 wrate
      external wrate
      real*8 pa,ptop,pb,wa,wtop,wb,wabw,wtopbw,wbbw
      real*8 kres
      real*8 ala,alb,al
      real*8 nga
      real*8 e
      real*8 tmpnorm
      integer i,j
      logical good,overlapping
      real*8 wresratio, wdev
      real*8 ngammamax
      real*8 dsrdbreit_wigner
      ! JE FIX: automatically set wresratio?
      parameter(wresratio=1.5d0) ! rate of on/off peak to treat as important
      parameter(wdev=0.005d0)   ! deviation from fit to calculation to include
      parameter(ngammamax=30.0d0) ! maximum number of widths to include

c...This code is copied from dsrdtab      
      pmax=mco(1)*umax*sqrt((1.0d0+0.25d0*umax**2/xmin)/xmin)
      pmax=min(10.0d0*mco(1),pmax)
      
c...Loop through the resonances
c      write(*,*) 'Checking resonances: ',nres
      do i=1,nres
         nga=0.5d0! 1.d0               ! number of gamma off peak to start with
         if (rgev(i).lt.2.0d0*mco(1)) then
            resinc(i)=.false.
            goto 99 ! end of loop, take next i
         endif
         overlapping=.false.
         do j=1,nres
            if (i.ne.j) then
               if (abs(rgev(j)-rgev(i))/rwid(i).lt.2.d0)
     &              overlapping=.true.
            endif
         enddo
         if (overlapping) then
            resinc(i)=.false.
            goto 99             ! end loop, take next i
         endif
            
         ptop=sqrt(rgev(i)**2/4-mco(1)**2)
         resinc(i)=.true.       ! default is to include as difficult point
         if (ptop.lt.pmax) then ! resonance in interesting p range
            e=rgev(i)-nga*rwid(i)
            if (e.lt.2.d0*mco(1)) then
               resfit(i)=.false.
               goto 99
            endif
            pa=sqrt((rgev(i)-nga*rwid(i))**2/4-mco(1)**2)
            pb=sqrt((rgev(i)+nga*rwid(i))**2/4-mco(1)**2)
            wa=wrate(-pa) ! -pa to force rate routine call, i.e. no interpol.
            wtop=wrate(-ptop)
            wb=wrate(-pb)
            if (wtop/min(wa,wb).gt.wresratio) then ! include resonance
               resinc(i)=.true.
c...Now fit it
               wabw=dsrdbreit_wigner(rgev(i),rwid(i),0.d0,mco(1),pa)
               wtopbw=dsrdbreit_wigner(rgev(i),rwid(i),0.d0,mco(1),ptop)
               wbbw=dsrdbreit_wigner(rgev(i),rwid(i),0.d0,mco(1),pb)
               tmpnorm=wtop/wtopbw ! only fit to top point
               if ((wa+wb).gt.(wabw*tmpnorm+wbbw*tmpnorm)) then ! need constant
                  tmpnorm=(wa+wb-2.d0*wtop)/(wabw+wbbw-2.d0*wtopbw)
                  kres=wtop-wtopbw*tmpnorm
               else ! no constant fit
                  ! tmpnorm=wtop/wtopbw ! already done above
                  kres=0.d0
               endif
c               write(*,*) 'Resonance number ',i
c               write(*,*) 'Before alpha fit: tmpnorm = ',tmpnorm
c               write(*,*) 'const: ',kres
c               write(*,*) 'a:   ',pa,wa,wabw*tmpnorm+kres
c               write(*,*) 'top: ',ptop,wtop,wtopbw*tmpnorm+kres
c               write(*,*) 'b:   ',pb,wb,wbbw*tmpnorm+kres
               ala=log(wa/(wabw*tmpnorm+kres))/log(pa/ptop)
               alb=log(wb/(wbbw*tmpnorm+kres))/log(pb/ptop)
               al=(ala+alb)/2.d0
c               write(*,*) 'ala alb al: ',ala,alb,al
               tmpnorm=tmpnorm/ptop**al
               wabw=dsrdbreit_wigner(rgev(i),rwid(i),al,mco(1),pa)
               wtopbw=dsrdbreit_wigner(rgev(i),rwid(i),al,mco(1),ptop)
               wbbw=dsrdbreit_wigner(rgev(i),rwid(i),al,mco(1),pb)
c               write(*,*) 'After alpha fit: tmpnorm = ',tmpnorm
c               write(*,*) 'const: ',kres
c               write(*,*) 'a:   ',pa,wa,wabw*tmpnorm+kres
c               write(*,*) 'top: ',ptop,wtop,wtopbw*tmpnorm+kres
c               write(*,*) 'b:   ',pb,wb,wbbw*tmpnorm+kres
c...Check if this is good enough
               good=.true.
               if (abs((wabw*tmpnorm+kres-wa)/wa).gt.wdev) good=.false.
               if (abs((wbbw*tmpnorm+kres-wb)/wb).gt.wdev) good=.false.
               ! JE FIX. Add check on reasonable alpha
c               write(*,*) abs((wabw*tmpnorm+kres-wa)/wa),wdev
c               write(*,*) abs((wbbw*tmpnorm+kres-wb)/wb),wdev
c               write(*,*) good
               if (good) then       ! Good fit, find out for how many gamma it is good
                  resfit(i)=.true.
                  resalpha(i)=al
                  resnorm(i)=tmpnorm
                  resgamma(i)=nga
                  resconst(i)=kres
 10               nga=nga*1.5d0 ! nga+1.d0  ! Might want to change this, nga*1.5d0?
                  e=rgev(i)-nga*rwid(i)
                  if (e.lt.2.d0*mco(1)) goto 88 ! exit here
                  pa=sqrt((rgev(i)-nga*rwid(i))**2/4-mco(1)**2)
                  pb=sqrt((rgev(i)+nga*rwid(i))**2/4-mco(1)**2)
c                  write(*,*) 'pa pb: ',pa,pb
                  wa=wrate(-pa)
                  wb=wrate(-pb)
                  wabw=dsrdbreit_wigner(rgev(i),rwid(i),al,mco(1),pa)
                  wbbw=dsrdbreit_wigner(rgev(i),rwid(i),al,mco(1),pb)
                  good=.true.
                  if (abs((wabw*tmpnorm+kres-wa)/wa).gt.wdev) good=.false.
                  if (abs((wbbw*tmpnorm+kres-wb)/wb).gt.wdev) good=.false.
c                  write(*,*) '  ',abs((wabw*tmpnorm-wa)/wa),wdev
c                  write(*,*) '  ',abs((wbbw*tmpnorm-wb)/wb),wdev
c                  write(*,*) '  ',good
                  if (good) then
                     resgamma(i)=nga
                     if (nga.lt.ngammamax) goto 10
                  endif
 88               continue
               else
                  resfit(i)=.false.
               endif
            else
               resinc(i)=.false.
            endif
         endif
         if (resinc(i).and.resfit(i)) then
            if (prtlevel.ge.3) then
               write(*,*) 'dsrdbw_setup: Found a resonance, i=',i
               write(*,*) 'sqrt(s)=',rgev(i),' width=',rwid(i)
               write(*,*) 'alpha=',resalpha(i),' ngamma=',resgamma(i)
               write(*,*) 'const=',resconst(i)
               pa=sqrt((rgev(i)-resgamma(i)*rwid(i))**2/4-mco(1)**2)
               pb=sqrt((rgev(i)+resgamma(i)*rwid(i))**2/4-mco(1)**2)
               write(*,*) 'BW with good fit in p-range: ',pa,pb
            endif
         endif
 99      continue
      enddo


      return
      end
      
         
         
            
         


         
         
