*******************************************************************************
*** Subroutine dshealpixint integrates a function f(l,b) on the unit sphere ***
*** by means of pixelization in the HEALpix scheme, where l and b are       ***
*** longitude and latitude in rad, respectively. The integration is         ***
*** optimized for functions f centered on (l,b)=0.0d0 (but not necessarily  ***
*** spheically symmetric).                                                  ***       
***                                                                         ***
***  Input:                                                                 ***
***    f      - function to be integrated over; must be declared EXTERNAL   ***
***             in calling main program / subroutine                        ***
***    nmin   - HEALPIX levels to start integration (at least 5 recomm.)    ***
***    nmax   - HEALPIX level where integration is stopped even if accuracy ***
***             goal is not met                                             ***
***    eps    - relative accuracy goal                                      ***
***    epsabs - absolute accuracy goal                                      ***
***                                                                         ***
***  Output:                                                                ***
***    res   - integration result                                           ***
***    niter - number of iterations (HEALPIX levels) actually used          ***
***    ierr  - 1 if nmax is reachd before accuracy goal, else 0             ***
***                                                                         ***
*** author: Torsten.Bringmann.fys.uio.no (based on previous DS version)     ***
*** date 2018-11-28                                                         ***
*******************************************************************************
      subroutine dshealpixint(f,nmin,nmax,eps,epsabs,res,niter,ierr)
      implicit none
      include 'dsio.h'
      real*8 f,res
      external f
      integer ierr,niter
      real*8 aveo,ave,eps,epsabs,FourPI
      integer n,nmin,nmax
c      parameter (nmin=4,nmax=400)
c      parameter (eps=1.d-4,epsabs=1.d-4)
      ierr=0
      FourPI = 16.0 * datan(1.d0)      
      call dshealpixave(f,nmin,0,aveo)
      do n=nmin+1,nmax
         call dshealpixave(f,n,1,ave)
         if (prtlevel.gt.1) write(*,*) 'dshealpixint: n, res = ', n,ave*FourPI
         if (abs(ave-aveo).lt.eps*abs(ave)+epsabs) then
            res=ave*FourPI
            niter=n
            return
         endif
         aveo=ave
      enddo
      niter=n
      res=ave*FourPI
      if (prtlevel.gt.0) 
     &    write (*,*) 'dshealpixint: max number of sides reached ',nmax
      ierr=1
      return
      end
