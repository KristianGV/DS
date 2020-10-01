*******************************************************************************
*** Subroutine dshealpixave sums over the value of EXTERNAL function f(l,b) ***
*** at the center of all HEALPIX pixels of level n, divided by their size.  ***
*** longitude and latitude in rad, respectively. The integration is         ***
***                                                                         ***
***  Input:                                                                 ***
***    f   - EXTERNAL function of longitude l and latitude b (both in rad)  ***
***    n   - HEALPIX level                                                  ***
***    how - sum over all pixels at level n (how=0), or only those that     ***
***          are neighbouring or inside those that gave a non-zero value    ***
***          at the revious call (How=1; this assmues that the previous     ***
***          call was done with HEALPIX level n-1                           ***
***                                                                         ***
***  Output:                                                                ***
***    ave   - average value of f on unit sphere                            ***
***                                                                         ***
*** author: Torsten.Bringmann.fys.uio.no                                    ***
*** date 2018-11-28                                                         ***
*******************************************************************************
      subroutine dshealpixave(f,n,how,ave)
      use pix_tools
      implicit none
      include 'dsmpconst.h'
      integer n, how
      real*8 f,ave
      integer i,ppmax, nside
      real*8 sum,tmp,l,b,l0,b0
      integer plist(2000000),nlist(2000000),pmax,nmax
      save nlist, nmax

      nside=2**n      
      ppmax=12*nside*nside
      sum=0.d0

c      goto 10 ! uncomment to avoid any bookkeeping, which is faster if you
               ! always want to use how=0. This is typically better for
               ! very large regions with Omega >~ 4*pi / 10.

c... if this is the first time we call this routine, we need to look
c... at all points   
      if (how.eq.0) then
         nmax=12*nside*nside
         do i=1,nmax
           nlist(i) = i-1
         enddo
      endif
c... for subsequent calls, only the potentially relevant points
c... contained in plist are considered (saved in nlist from previous call)    
      pmax=nmax
      do i=1,nmax
        plist(i) = nlist(i)
      enddo

      nmax=0
      call pix2ang_nest(nside, 0, l0, b0)  ! to rotate the HP system to (0,0)
      do i=1,pmax
        call pix2ang_nest(nside, plist(i), l, b)
        tmp=f(-l+l0,b-b0)
        if (tmp.ne.0.0d0) then
          call dsaddhpxpt(nside,plist(i),nmax,nlist)      
          sum=sum+tmp
        endif
      enddo
      
      goto 20 ! skip the "always how=0" implementation
10    sum=0.d0
      call pix2ang_nest(nside, 0, l0, b0)  ! to rotate the HP system to (0,0)
      do i=1,ppmax
         call pix2ang_nest(nside, i-1, l, b)
         sum=sum+f(-l+l0,b-b0)
      enddo
20    ave=sum/real(ppmax)
      end



*******************************************************************************
*** Auxiliary routine needed by dshealpixave, adding points to nlist        ***
*******************************************************************************
      subroutine dsaddhpxpt(nside,ipix,nmax,nlist)
      use pix_tools
      implicit none
        integer nside,ipix,nmax,nlist(10000)
        integer i,j,nneigh,lneigh(9),cand(36)
        integer klo,khi,k
        
c... create a list with 4*(nneigh+1) candidate pixels to be added 
c... to the list of relevant pixels at next level
        call neighbours_nest(nside, ipix, lneigh, nneigh)        
        lneigh(nneigh+1) = ipix ! add central point to list
        do i=0,nneigh
          do j=0,3
            cand(4*i+j+1)=4*lneigh(i+1)+j
          enddo
        enddo

c... now sort these points into the existing nlist
        do i=1,4*(nneigh+1)
c... first determine where the new point would fit        
          if (nmax.eq.0) then
            k=1
          elseif(cand(i).ge.nlist(nmax)) then
            if (cand(i).eq.nlist(nmax)) goto 50 ! point already exists
            k=nmax+1
          elseif(cand(i).le.nlist(1)) then
            if (cand(i).eq.nlist(1)) goto 50 ! point already exists
            k=1
          else
            klo=1
            khi=nmax
 10         if (khi-klo.gt.1) then
              k=(khi+klo)/2
              if (nlist(k).gt.cand(i)) then
                 khi=k
              elseif (nlist(k).lt.cand(i)) then
                 klo=k
              else
                 goto 50 ! point already exists
              endif
              goto 10
            endif
            k=khi
          endif
c... then add new point at location k
          do j=nmax,k,-1
            nlist(j+1)=nlist(j)
          enddo
          nlist(k)=cand(i)
          nmax=nmax+1          
 50       continue        
        enddo
      return  
      end
