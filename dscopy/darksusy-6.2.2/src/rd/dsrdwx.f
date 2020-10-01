      real*8 function dsrdwx(p)
c_______________________________________________________________________
c  Wrapper routine for dsanwx. This routine will call dsanwx to get the
c  invariant annihilation rate, but instead of calling it directly it will
c     a) call dsanwx when needed and store results
c     b) use previous results (with interpolation) if close enough in
c        momentum
c     c) consider resonances and thresholds when it determines if
c        a call to dsanwx is needed instead of tabulation
c     d) use a Breit-Wigner for for resonances iff it gives a very
c        good fit to the resonance (this is setup with dsrdbw_setup)      
c  input:
c    p - momentum (GeV) in center of momentum frame
c      if p<0, invariant rate will be calculated for |p| and a call
c      to dsanwx is forced, i.e. it will not interpolate tabulated results
c  common:
c    'dsrdcom.h' - included common blocks
c  output:
c    dsrdwx - invariant annihilation rate (i.e. the same as dsanwx)
c  Author: Joakim Edsjo (edsjo@fysik.su.se), 2018
c=======================================================================
      implicit none
      include 'dsrdcom.h'
      real*8 p,ppl,ptmp
      real*8 dsanwx
      real*8 dpmin,rrr,mindp,mindp_crit
      real*8 mindp_res,mindp_th
      real*8 sqrts,tmp,tmp1,tmp2
      real*8 lx1,lx2,ly1,ly2
      real*8 wx1,wx2
      real*8 corr
      integer i,j,ifound
      logical found,force
      save dpmin
      real*8 tmpwx
      integer istat
      real*8 dsrdbw_get
      real*8 reslin,err,dsrdquad

c-----------------------------------------------------------------------

      force=.false. ! 
      
      if (p.lt.0.d0) then
         force=.true.
      endif
      ptmp=abs(p)
      
c...Determine correction factor for point insertion checks
c...The correction factor is 1 at mco(1)/3 and then increases further away
c...to a maximum of 5
      if (ptmp.lt.mco(1)/3.d0) then
         corr=exp(-3.d0*(ptmp/mco(1)-0.33d0)**2)
         corr=min(2.d0,1.d0/corr)
      elseif (ptmp.gt.2.d0*mco(1)/3.d0) then
         corr=exp(-4.d0*(ptmp/mco(1)-0.67d0)**2)
         corr=min(5.d0,1.d0/corr)
      else
         corr=1.d0
      endif
      
c      write(62,*) p/mco(1),corr
      
      sqrts=sqrt(mco(1)**2 + ptmp**2)*2.d0
      if (nrd.eq.0) then        ! first call for this model         
         dpmin=dpminr*mco(1)*corr
         ppp(0)=0.d0
c...Set up two endpoints as a starting point
         ppp(1)=0.d0
         rdwx(1)=dsanwx(ppp(1))
         ppp(2)=mco(1)*5.d0
         rdwx(2)=dsanwx(ppp(2))
         nrd=2
      endif

      tmpwx=dsrdbw_get(ptmp,istat)
      if (istat.eq.0) then      ! found resonance fit, use this
         dsrdwx=tmpwx
         return
      endif
      
c...Find smallest distance to already existing point      
c      mindp=mco(1)*1.d10
c      do i=1,nrd
c         tmp=abs(p-ppp(i))
c         if (tmp.lt.mindp) then
c            mindp=tmp
c            if ((p-ppp(i)).ge.0.d0) then
c               ifound=i
c            else
c               ifound=i-1
c            endif
c         endif
c      enddo

      call dshunt(ppp(1),nrd,ptmp,ifound)
      if (ptmp.eq.ppp(ifound+1)) ifound=ifound+1 ! if p=upper
      tmp1=ptmp-ppp(ifound)
      tmp2=ppp(ifound+1)-ptmp
      mindp=min(tmp1,tmp2)

      if (force) then
         if (mindp.eq.0.d0) then
            found=.true.
            goto 90
         endif
         goto 100
      endif
      
c...Now determine criterion on mindp, mindp_crit. This is the one mindp will
c...be compared with later      
      
c...General statement
      if (ptmp.le.mco(1)) then
         mindp_crit=dpmin
      else
         mindp_crit=dpmin
      endif

c...Now check resonances
      do i=1,nres
         tmp=abs(sqrts-rgev(i))
c...mindp_res is now set to be small at res. and then increase as we
c...move away from the peak
c...For narrow resonances we might need to play with how fast it grows
c...We need to go to maybe 50 widths or even more.         
c         mindp_res=0.05d0*rwid(i)+0.25d0*tmp ! *(tmp/rwid(i))**0.2d0
         mindp_res=0.05d0*rwid(i)+0.10d0*tmp ! add 0.1gamma per gamma
c         if (p.gt.mco(1)) mindp_res=mindp_res*2.d0
c         if (tmp/rwid(i).lt.100.d0) then
            mindp_crit=min(mindp_res,mindp_crit) ! JE , no corr here
c         endif
      enddo   

      
c...Now check thresholds
      do i=1,nth
         tmp=abs(ptmp-pth(i))
c...Determine minimum p for thresholds, mindp_th
         mindp_th=max(dpmin/100.d0,0.5d0*tmp)
         if (tmp/mco(1).lt.0.05d0) then ! only do this close to threshold
            mindp_crit=min(mindp_th*corr,mindp_crit) ! JE, corr added
         endif
      enddo
      
c...Now check if there is reason to insert a new point
c...Defult is to not insert. Now check if there is reason to
c...do so anyway      
      found=.true.

c...Now check against criterion
      if (mindp.gt.mindp_crit) found=.false.

c...If found=.false., check if interpolation between nearest and next-nearest
c...neighboors is more or less the same in which case we probably won't
c...need to add a new point
      if ((.not.found).and.ifound.gt.1.and.ifound.lt.nrd-1) then      ! check if we really need to add this point
         ppl=(ptmp-ppp(ifound))/(ppp(ifound+1)-ppp(ifound))
         wx1=rdwx(ifound)*(1.0d0-ppl)+rdwx(ifound+1)*ppl
         ppl=(ptmp-ppp(ifound-1))/(ppp(ifound+2)-ppp(ifound-1))
         wx2=rdwx(ifound-1)*(1.0d0-ppl)+rdwx(ifound+2)*ppl
         if (abs(wx2-wx1)/wx1.lt.0.01d0) found=.true.
      endif

c...add point if above last tabulated point
      if (ptmp.gt.ppp(nrd)) found=.false.


 90   continue  
      
      if (found) then           ! interpolate
         if (ifound+2.le.nrd) then  ! can go up two steps for quad interpolation
            dsrdwx=dsrdquad(ppp(ifound),rdwx(ifound),
     &                      ppp(ifound+1),rdwx(ifound+1),
     &           ppp(ifound+2),rdwx(ifound+2),ptmp,reslin,err)
         else ! need to go down for interpolation
            dsrdwx=dsrdquad(ppp(ifound-1),rdwx(ifound-1),
     &                      ppp(ifound),rdwx(ifound),
     &           ppp(ifound+1),rdwx(ifound+1),ptmp,reslin,err)
         endif
c...Option to use linear interpolation if error is large, not used
c         if (abs(err).gt.0.05d0) then ! use linear interpolation instead
c            dsrdwx=reslin
c         endif
            
c         write(60,*) p,dsrdwx,reslin,err
         
            
c...Lin-lin interpolation        
c         ppl=(ptmp-ppp(ifound))/(ppp(ifound+1)-ppp(ifound))
c         dsrdwx=rdwx(ifound)*(1.0d0-ppl)+rdwx(ifound+1)*ppl
c...Log-log
c         lx1=log(max(ppp(ifound),1.d-10))
c         lx2=log(ppp(ifound+1))
c         ly1=log(rdwx(ifound))
c         ly2=log(rdwx(ifound+1))
c         ppl=(log(ptmp)-lx1)/(lx2-lx1)
c         dsrdwx=exp(ly1*(1.d0-ppl)+ly2*ppl)
c...Lin-log
c         ly1=log(rdwx(ifound))
c         ly2=log(rdwx(ifound+1))         
c         ppl=(ptmp-ppp(ifound))/(ppp(ifound+1)-ppp(ifound))
c         dsrdwx=exp(ly1*(1.d0-ppl)+ly2*ppl)
         return
      endif

c...Not found, or forced to add, make new call and add to table
 100  nrd=nrd+1
      if (nrd.gt.nrmax) then
         write(*,*) 'DS ERROR in dsrdwx:'
         write(*,*) '  Tried to add more tabulation points than the',
     &        ' allowed maximum: ',nrmax
         write(*,*) '  You need to increase nrmax in dsrdcom.h'
         write(*,*) '  DarkSUSY stopping.'
         do i=1,nrd
            write(47,*) i,ppp(i),rdwx(i)
         enddo
         
         stop
      endif
      ppp(nrd)=ptmp
      rrr=dsanwx(ptmp)
      rdwx(nrd)=rrr

c...Add this entry to table and mare sure it is sorted
      if (nrd.gt.2) then 
         do i=0,nrd-2
            if (ppp(i).lt.ppp(nrd).and.ppp(i+1).gt.ppp(nrd)) then ! here
               do j=nrd,i+2,-1
                  ppp(j)=ppp(j-1)
                  rdwx(j)=rdwx(j-1)
               enddo
               ppp(i+1)=ptmp
               rdwx(i+1)=rrr
               goto 20
            endif
         enddo
 20      continue
      else if (nrd.eq.2) then
         if (ppp(2).lt.ppp(1)) then ! change order of first points
            ppp(2)=ppp(1)
            rdwx(2)=ppp(1)
            ppp(1)=ptmp
            rdwx(1)=rrr
         endif
      endif
      
      dsrdwx=rrr

      return
      end
