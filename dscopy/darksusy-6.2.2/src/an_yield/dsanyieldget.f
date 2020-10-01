*****************************************************************************
*** function dsanyieldget gives the information in the differential and
*** integrated arrays phiit and phidiff for given array indices.
*** compared to getting the results directly from the array, this
*** routine performs smoothing if requested.
*** the smoothing is controlled by the parameter ansmooth in the
*** following manner.
*** ansmooth = 0 - no smoothing
***            1 - smoothing of zi-1,zi and zi+1 bins if z>.3
***            2 - smoothing of zi-2,zi-1,zi,zi+1 and zi+2 if z>0.3
*** Author: Joakim Edsjo, edsjo@fysik.su.se
*****************************************************************************

      real*8 function dsanyieldget(zi,mxi,ch,fi,ftype,istat)
      implicit none
      include 'dsanyieldcom.h'

c------------------------ variables ------------------------------------

      integer zi,mxi,ch,fi,ftype,istat,zsmstart,zimin,kind
      real*8 dsanyieldgetint,dsanyieldgetdiff

c-----------------------------------------------------------------------

      if (ftype.le.2) then
         kind=1                 ! yields
      else
         kind=2                 ! error on yields
      endif

      dsanyieldget=0.d0
      zsmstart = int(0.3*zn)
      if (ftype.eq.1) then
        zimin=0
      else
        zimin=-1
      endif

      if (zi.lt.zimin) then
        write(*,*) 'error in dsanyieldget: zi = ',zi
        write(*,*) '  fi = ',fi
        write(*,*) '  ftype = ',ftype
        stop
      endif

      if (fi.lt.9.or.fi.gt.13) then ! not dbar
        if (zi.gt.zn) then
          write(*,*) 'error in dsanyieldget: zi = ',zi
          write(*,*) '  fi = ',fi
          write(*,*) '  ftype = ',ftype
          stop
        endif
      else
        if (zi.gt.zndb) then
          write(*,*) 'error in dsanyieldget: zi = ',zi
          write(*,*) '  fi = ',fi
          write(*,*) '  ftype = ',ftype
          stop
        endif
      endif

      if (ansmooth.eq.0.or.zi.le.zsmstart) then
  
        if (ftype.eq.1.or.ftype.eq.3) then
          dsanyieldget=dble(dsanyieldgetint(zi,mxi,ch,fi,kind))
        else
          dsanyieldget=dble(dsanyieldgetdiff(zi,mxi,ch,fi,kind))
        endif

      elseif (ansmooth.eq.1) then

        if (ftype.eq.1) then
          if (zi.ge.1.and.zi.le.(zn-1)) then
             dsanyieldget=
     &        .25d0*dble(dsanyieldgetint(zi-1,mxi,ch,fi,kind))
     &        +0.5d0*dble(dsanyieldgetint(zi,mxi,ch,fi,kind))
     &        +0.25d0*dble(dsanyieldgetint(zi+1,mxi,ch,fi,kind))
          elseif (zi.eq.0) then
             dsanyieldget=
     &        .75d0*dble(dsanyieldgetint(zi,mxi,ch,fi,kind))+
     &        0.25d0*dble(dsanyieldgetint(zi+1,mxi,ch,fi,kind))
          elseif (zi.eq.zn) then
             dsanyieldget=
     &        .75d0*dble(dsanyieldgetint(zi,mxi,ch,fi,kind))+
     &        0.25d0*dble(dsanyieldgetint(zi-1,mxi,ch,fi,kind))
          endif
        else
          if (zi.ge.1.and.zi.le.(zn-1)) then
             dsanyieldget=
     &        .25d0*dble(dsanyieldgetdiff(zi-1,mxi,ch,fi,kind))
     &        +0.5d0*dble(dsanyieldgetdiff(zi,mxi,ch,fi,kind))
     &        +0.25d0*dble(dsanyieldgetdiff(zi+1,mxi,ch,fi,kind))
          elseif (zi.eq.0) then
             dsanyieldget=
     &        .75d0*dble(dsanyieldgetdiff(zi,mxi,ch,fi,kind))+
     &        0.25d0*dble(dsanyieldgetdiff(zi+1,mxi,ch,fi,kind))
          elseif (zi.eq.zn) then
             dsanyieldget=
     &        .75d0*dble(dsanyieldgetdiff(zi,mxi,ch,fi,kind))+
     &        0.25d0*dble(dsanyieldgetdiff(zi-1,mxi,ch,fi,kind))
          endif
        endif

      else  ! ansmooth.eq.2

        if (ftype.eq.1) then
          if (zi.le.(zn-2)) then
             dsanyieldget=
     &        .10d0*dble(dsanyieldgetint(zi-2,mxi,ch,fi,kind))
     &        +0.225d0*dble(dsanyieldgetint(zi-1,mxi,ch,fi,kind))
     &        +0.350d0*dble(dsanyieldgetint(zi,mxi,ch,fi,kind))
     &        +0.225d0*dble(dsanyieldgetint(zi+1,mxi,ch,fi,kind))
     &        +0.100d0*dble(dsanyieldgetint(zi+2,mxi,ch,fi,kind))
          elseif (zi.eq.zn-1) then
             dsanyieldget=
     &        .10d0*dble(dsanyieldgetint(zi-2,mxi,ch,fi,kind))
     &        +0.225d0*dble(dsanyieldgetint(zi-1,mxi,ch,fi,kind))
     &        +0.450d0*dble(dsanyieldgetint(zi,mxi,ch,fi,kind))
     &        +0.225d0*dble(dsanyieldgetint(zi+1,mxi,ch,fi,kind))
          else
             dsanyieldget=
     &        .10d0*dble(dsanyieldgetint(zi-2,mxi,ch,fi,kind))
     &        +0.225d0*dble(dsanyieldgetint(zi-1,mxi,ch,fi,kind))
     &        +0.675d0*dble(dsanyieldgetint(zi,mxi,ch,fi,kind))
          endif
        else
          if (zi.le.(zn-2)) then
             dsanyieldget=
     &        .10d0*dble(dsanyieldgetdiff(zi-2,mxi,ch,fi,kind))
     &        +0.225d0*dble(dsanyieldgetdiff(zi-1,mxi,ch,fi,kind))
     &        +0.350d0*dble(dsanyieldgetdiff(zi,mxi,ch,fi,kind))
     &        +0.225d0*dble(dsanyieldgetdiff(zi+1,mxi,ch,fi,kind))
     &        +0.100d0*dble(dsanyieldgetdiff(zi+2,mxi,ch,fi,kind))
          elseif (zi.eq.zn-1) then
             dsanyieldget=
     &        .10d0*dble(dsanyieldgetdiff(zi-2,mxi,ch,fi,kind))
     &        +0.225d0*dble(dsanyieldgetdiff(zi-1,mxi,ch,fi,kind))
     &        +0.450d0*dble(dsanyieldgetdiff(zi,mxi,ch,fi,kind))
     &        +0.225d0*dble(dsanyieldgetdiff(zi+1,mxi,ch,fi,kind))
          else
             dsanyieldget=
     &        .10d0*dble(dsanyieldgetdiff(zi-2,mxi,ch,fi,kind))
     &        +0.225d0*dble(dsanyieldgetdiff(zi-1,mxi,ch,fi,kind))
     &        +0.675d0*dble(dsanyieldgetdiff(zi,mxi,ch,fi,kind))
          endif

        endif

      endif

      end
