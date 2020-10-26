*******************************************************************************
*** Function dsInterpolateTable takes tabulated values of a real function,  ***
*** (x_i, f(x_i)), and returns an interpolated value for f(x).              ***
*** NB: This routine assumes that the x-values are ordered as x_1<x_2<... ! ***
***                                                                         ***
***  Input:                                                                 ***
***    x        - independent variable                                      ***
***    functab  - Array containing (x_i, f(x_i))                            ***
***    Npoints  - size of functab                                           ***
***    tabID    - integer from 1-10                                         ***
***    how      - linear interpolation (how=1)                              ***
***               log-interpolation  (how=2)                                ***
***                                                                         ***
***  Output: f(x)                                                           ***
***                                                                         ***
***  Note: 'tabID' does not affect the result, but increases performance in ***
***        case of subsequent calls to dsInterpolateTable with similar      ***
***        values of x. If dsInterpolateTable is simultaneously used on     ***
***        different tables, each of those tables should be assigned a      ***
***        unique (and different) tabID.                                    ***
***                                                                         ***
*** author: Torsten.Bringmann.fys.uio.no                                    ***
*** date 2019-10-09                                                         ***
*******************************************************************************
      real*8 function dsInterpolateTable(x,functab,Npoints,tabID,how)
      implicit none
  
      integer Npoints, tabID, how
      real*8 x, functab(2,Npoints) 
      real*8 T,lnt
      ! integer k, dk, klo(tabIDmax), khi(tabIDmax)
      integer maxdat, k, dk, klo(2), khi(2)
      parameter (maxdat=100)
      
      integer tabIDmax
      parameter (tabIDmax=10)

      logical initialized
      data initialized /.false./

      real*8 res, LogTinTzdat(maxdat,2,2)
      common /soilcom/ LogTinTzdat
      save /soilcom/
      ! logical initialized
      ! data initialized /.false./
      save initialized, klo, khi




  
      ! real*8


      if (tabID.lt.1.or.tabID.gt.tabIDmax) then
        write(*,*) 'Fatal ERROR in dsInterpolateTable: tabID = ',tabID
        write(*,*) 'Please provide a value between 1 and ',tabIDmax
        stop
      endif

      ! integer how  ! 1: T=T0->Tz; 2: T=Tz->T0
      

      lnt = log(T)  
      res = 0.0d0
      if (.not.initialized) then
        klo(1) = 1
        klo(2) = 1
        khi(1) = maxdat
        khi(2) = maxdat
        initialized = .true.
      endif

      if (lnt.le.LogTinTzdat(1,how,1)) then
         klo(how) = 1
         res = LogTinTzdat(1,how,2)
      elseif (lnt.ge.LogTinTzdat(maxdat,how,1)) then
         khi(how) = maxdat
         res = LogTinTzdat(maxdat,how,2)
      else
        dk=maxdat/50
 100    if (LogTinTzdat(khi(how),how,1).lt.lnt) then
           khi(how) = khi(how) + dk
           dk = 2*dk
           if (khi(how).lt.maxdat) goto 100
           khi(how)=maxdat
        endif
        dk=maxdat/50
 110    if (LogTinTzdat(klo(how),how,1).gt.lnt) then
           klo(how) = klo(how) - dk
           dk = 2*dk
           if (klo(how).gt.1) goto 110
           klo(how)=1
        endif
 120    if (khi(how)-klo(how).gt.1) then
           k=(khi(how)+klo(how))/2
           if (LogTinTzdat(k,how,1).gt.lnt) then
              khi(how)=k
           else
              klo(how)=k
           endif
           goto 120
        endif
         res = LogTinTzdat(klo(how),how,2)+
     &               (LogTinTzdat(khi(how),how,2)-LogTinTzdat(klo(how),how,2))
     &               *(lnt-LogTinTzdat(klo(how),how,1))/
     &               (LogTinTzdat(khi(how),how,1)-LogTinTzdat(klo(how),how,1))
      endif
      dsInterpolateTable = exp(res)
      return
      end ! dsgetTzratio
      
