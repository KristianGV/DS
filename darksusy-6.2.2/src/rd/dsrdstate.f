**********************************************************************
*** dsrdstate saves or resets the common blocks in dsrdcom.h
*** It is used to make dsrdomega flexible by allowing for a
*** fast flag to conveniently change some parameters, without
*** permanently overwriting common block variables.
***
*** Input: c='save'  saves the values of the common block varaiables
***                  in dsrdcom
***        c='reset' resets the common block variables in dsrdcom
***                  from saved values
*** Author: Joakim Edsjo, edsjo@fysik.su.se
*** Date: June 19, 2019
**********************************************************************

      subroutine  dsrdstate(c)
      implicit none
      include 'dsrdcom.h'
      real*8 rpar(16)
      integer ipar(10)
      save rpar,ipar

      character*(*) c

      if (c.eq.'save') then

         rpar(1)=cosmin
         rpar(2)=waccd
         rpar(3)=dpminr
         rpar(4)=dpthr
         rpar(5)=wdiffr
         rpar(6)=wdifft
         rpar(7)=hstep
         rpar(8)=hmin
         rpar(9)=compeps
         rpar(10)=xinit
         rpar(11)=xfinal
         rpar(12)=umax
         rpar(13)=cfr
         rpar(14)=pdivr
         rpar(15)=dpres
         rpar(16)=rdt_max
         ipar(1)=thavint
         ipar(2)=rdluerr
         ipar(3)=rdlulog
         ipar(4)=rdprt
         ipar(5)=nlow
         ipar(6)=nhigh
         ipar(7)=npres
         ipar(8)=nthup
         ipar(9)=cthtest
         ipar(10)=spltest

      elseif (c.eq.'reset') then

         cosmin=rpar(1)
         waccd=rpar(2)
         dpminr=rpar(3)
         dpthr=rpar(4)
         wdiffr=rpar(5)
         wdifft=rpar(6)
         hstep=rpar(7)
         hmin=rpar(8)
         compeps=rpar(9)
         xinit=rpar(10)
         xfinal=rpar(11)
         umax=rpar(12)
         cfr=rpar(13)
         pdivr=rpar(14)
         dpres=rpar(15)
         rdt_max=rpar(16)
         thavint=ipar(1)
         rdluerr=ipar(2)
         rdlulog=ipar(3)
         rdprt=ipar(4)
         nlow=ipar(5)
         nhigh=ipar(6)
         npres=ipar(7)
         nthup=ipar(8)
         cthtest=ipar(9)
         spltest=ipar(10)

      else

         write(*,*) 'DS ERROR in dsrdstore: Incorrect option c=',c
         stop

      endif
      return
      end

