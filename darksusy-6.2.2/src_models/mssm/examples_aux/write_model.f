      subroutine write_model(lunit)
c
c     To write MSSM7 model parameters to unit lunit
c
      implicit none
      include 'dsidtag.h'
      include 'dsmssm.h'
      real*8 at,ab,mtop,mqtild
      integer lunit
      integer first
      save first
      data first/0/
      if (first.eq.0) then
         write (lunit,1001)
         first = 1
      endif
      mtop = mass(kt)
      mqtild=dsqrt(mass2q(3))
      at = asoftu(3)/mqtild
      ab = asoftd(3)/mqtild
      write (lunit,2001) 
     &     idtag,mu,m2,ma,tanbe,mqtild,at,ab
      return
 1001 format ('#','.....id.....',1x,'......mu......',1x,
     &     '......m2......',1x,'......ma......',1x,
     &     '.....tanbe....',1x,'......m.......',1x,
     &     '.....at/m.....',1x,'.....ab/m.....')
 2001 format (1x,a12,7(1x,e14.8))
      end

