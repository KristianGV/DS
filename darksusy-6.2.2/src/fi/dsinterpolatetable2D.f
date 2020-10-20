*******************************************************************************
*** Function dsInterpolateTable2D takes tabulated values of a 2D real       ***
*** function, f(x_i,y_j), and returns an interpolated value  for f(x,y).    ***
***                                                                         ***
***  Input:                                                                 ***
***    x, y     - independent variables                                     ***
***    functab  - Array containing f(x_i,y_j)                               ***
***    Nx, Ny   - size of functab (Nx x Ny)                                 ***
***    Nymax    - 1D array containing the number of relevant bins in y,     ***
***               for a given value of x_i  (length: Nx)                    ***
***    xtab     - 1D array containing x_i  (length: Nx)                     ***
***    ytab     - 2D array containing y_ij (length: Nx x Ny)                ***
***    tabID    - integer from 1-15 (see below)                             ***
***    how      - linear interpolation (how=1)                              ***
***               log-interpolation  (how=2)                                ***
***                                                                         ***
***  Output: f(x,y)                                                         ***
***                                                                         ***
*** ASSUMPTIONS about input values (NOT tested, to increase performance!)   ***
***     - xtab is ordered as x(1) < x(2) < ... < x(Nx)                      ***
***     - for each value of x(i), there are at least a few tabulated values ***
***       of y, with y(i,1) < y(i,2) < ... < y(i,Nymax(i))                  ***
***     - Nymax(i) <= Ny                                                    ***
***     - any values for f(x_i,y_j), as well as y(i,j),                     ***
***       with y_ij > y(i,Nymax(i)) can be disregarded                      ***
***       (as indicated, Nymax(i) can be chosen arbitrary -- but            ***
***        performance is best for a uniform Nymax(i) = Ny)                 ***
***                                                                         ***
*** (Typical usage would make use of a wrapper function, creating a         ***
***  table satisfying the above assumptions. For an example, see            ***
***  dsanyield_contrib_A.)                                                  ***
***                                                                         ***
***  Note: 'tabID' does not affect the result, but increases performance    ***
***        in case of subsequent calls to dsInterpolateTable with similar   ***
***        values of x. If dsInterpolateTable is simultaneously used on     ***
***        different tables, each of those tables should be assigned a      ***
***        unique (and different) tabID.                                    ***
***        A call with tabID<0 will enforce a "re-set" of this label, to    ***
***        to ensure that the behaviour is identical to the first function  ***
***        call with abs(tabID).                                            ***
***                                                                         ***
*** author: Torsten.Bringmann.fys.uio.no                                    ***
*** date 2020-05-26                                                         ***
*******************************************************************************
      real*8 function dsInterpolateTable2D(x,y,functab,Nx,Ny,Nymax,xtab,ytab,
     &                                     tabID,how)
      implicit none
  
      integer Nx, Ny, Nymax(Nx), tabID, how
      real*8 x, y, functab(Nx,Ny), xtab(Nx), ytab(Nx,Ny)
  
      integer tabIDmax, tabIDp
      parameter (tabIDmax=15)
      logical initialized(tabIDmax)
      data initialized /tabIDmax*.false./

      integer i, j, k, dk, kx(2)
      integer xlo(tabIDmax), xhi(tabIDmax),ylo(2,tabIDmax), yhi(2,tabIDmax)
      real*8 res, resy(2), funcbest(2,2), xbest(2), ybest(2,2), xfit,yfit
      real*8 dslogcut
      save initialized, xlo, xhi, ylo, yhi

      tabIDp = abs(tabID)
      if (tabIDp.lt.1.or.tabIDp.gt.tabIDmax) then
        write(*,*) 'Fatal ERROR in dsInterpolateTable2D: tabID = ',tabID
        write(*,*) 'Please provide a value between 1 and ',tabIDmax
        write(*,*) '(or multiplied by -1, to re-set that tabel ID)'
        stop
      endif

c TB debug
c      write(*,*) 'input: x,y,Nx,Ny,tabID,how =',x,y,Nx,Ny,tabID,how

      res = 0.0d0
      if (tabID.lt.0.or.(.not.initialized(tabIDp))) then
        xlo(tabIDp) = 1
        xhi(tabIDp) = Nx
        ylo(1,tabIDp) = 1
        yhi(1,tabIDp) = Ny
        ylo(2,tabIDp) = 1
        yhi(2,tabIDp) = Ny
        initialized(tabIDp) = .true.
      endif

c... start by finding smallest intervall around x
      if (x.le.xtab(1)) then
         xlo(tabIDp) = 1
         xhi(tabIDp) = 1
      elseif (x.ge.xtab(Nx)) then
         xlo(tabIDp) = Nx
         xhi(tabIDp) = Nx
      else
        dk = 1
 100    if (x.gt.xtab(xhi(tabIDp))) then
           xhi(tabIDp) = xhi(tabIDp) + dk
           dk = 3*dk
           if (xhi(tabIDp).lt.Nx) goto 100
           xhi(tabIDp) = Nx
        endif
 110    if (x.lt.xtab(xlo(tabIDp))) then
           xlo(tabIDp) = xlo(tabIDp) - dk
           dk = 2*dk
           if (xlo(tabIDp).gt.1) goto 110
           xlo(tabIDp) = 1
        endif
 120    if (xhi(tabIDp)-xlo(tabIDp).gt.1) then
           k = (xlo(tabIDp)+xhi(tabIDp))/2
           if (x.lt.xtab(k)) then
              xhi(tabIDp) = k
           else
              xlo(tabIDp) = k
           endif
           goto 120
        endif
      endif

c TB debug
c      write(*,*) 'x, [xlo,xhi] = ',x,xtab(xlo(tabIDp)),xtab(xhi(tabIDp))
      

c... now find smallest intervall around y, both at xtab(xlo) and xtab(xhi)
      kx(1) = xlo(tabIDp)
      kx(2) = xhi(tabIDp)
      dk=1
      do i = 1,2
        if (y.le.ytab(kx(i),1)) then
           ylo(i,tabIDp) = 1
           yhi(i,tabIDp) = 1
        elseif (y.ge.ytab(kx(i),Nymax(kx(i)))) then
           ylo(i,tabIDp) = Nymax(kx(i))
           yhi(i,tabIDp) = Nymax(kx(i))
        else
  200     if (y.gt.ytab(kx(i),yhi(i,tabIDp))) then
             yhi(i,tabIDp) = yhi(i,tabIDp) + dk
             dk = 3*dk
             if (yhi(i,tabIDp).lt.Nymax(kx(i))) goto 200
             yhi(i,tabIDp) = Nymax(kx(i))
          endif
  210     if (y.lt.ytab(kx(i),ylo(i,tabIDp))) then
             ylo(i,tabIDp) = ylo(i,tabIDp) - dk
             dk = 2*dk
             if (ylo(i,tabIDp).gt.1) goto 210
             ylo(i,tabIDp) = 1
          endif
  220     if (yhi(i,tabIDp)-ylo(i,tabIDp).gt.1) then
             k = (ylo(i,tabIDp)+yhi(i,tabIDp))/2
             if (y.lt.ytab(kx(i),k)) then
                yhi(i,tabIDp) = k
             else
                ylo(i,tabIDp) = k
             endif
             goto 220
          endif
        endif
c TB debug
c        write(*,*) 'y, [ylo,yhi] = ',y,ytab(kx(i),ylo(i,tabIDp)),ytab(kx(i),yhi(i,tabIDp))
      enddo


c... These are the coordinates and function values at the four nearest tabulated points:
      do i = 1,2
        xbest(i)      = xtab(kx(i))
        ybest(i,1)    = ytab(kx(i),ylo(i,tabIDp))
        ybest(i,2)    = ytab(kx(i),yhi(i,tabIDp))
        funcbest(i,1) = functab(kx(i),ylo(i,tabIDp))
        funcbest(i,2) = functab(kx(i),yhi(i,tabIDp))
      enddo

c TB debug
c      do i=1,2
c        do j=1,2
c          write(*,*) i,j,xbest(i),ybest(i,j),funcbest(i,j)
c        enddo
c      enddo

c... convert everything to log (else a linear interpolation is performed)
      xfit = x
      yfit = y
      if (how.eq.2) then
        xfit = dslogcut(xfit)
        yfit = dslogcut(yfit)
        do j = 1,2
          xbest(j) = dslogcut(xbest(j))
          do i = 1,2
            ybest(i,j)    = dslogcut(ybest(i,j))
            funcbest(i,j) = dslogcut(funcbest(i,j))
          enddo
        enddo
      endif

c... first interpolate in y (for both values of x), then in x
      do i = 1,2
        resy(i) = funcbest(i,1)
        if (ybest(i,2).ne.ybest(i,1)) resy(i) = resy(i) +
     &           (funcbest(i,2)-funcbest(i,1))*(yfit-ybest(i,1))/
     &           (ybest(i,2)-ybest(i,1))
      enddo
      res = resy(1)
      if (xbest(2).ne.xbest(1)) res = res +
     &     (resy(2)-resy(1))*(xfit-xbest(1))/(xbest(2)-xbest(1))

      if (how.eq.2) res = exp(res)

      dsInterpolateTable2D = res

c TB debug
c      write(*,*) resy(1),resy(2),res


      return
      end

*********************************************************
*** auxiliary function to prevent too small or large logs
*********************************************************
      real*8 function dslogcut(x)
      implicit none
      real*8 x
      
      if (x.lt.1.d-50) then
        dslogcut = -115.129
      elseif (x.gt.1.d50) then
        dslogcut = 115.129
      else
        dslogcut = log(x)
      endif

      return
      end
