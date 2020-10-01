*******************************************************************************
*** Read in quenching factors and store in common block tables.             ***
***                                                                         ***
*** author: Torsten.Bringmann.fys.uio.no                                    ***
*** date 2018-07-04                                                         ***
*******************************************************************************
      subroutine dsddreadQF(filename1,filename2)
      implicit none
      character*(*) filename1,filename2
      character*200 scr
      include 'dsddcom.h'
      real*8 tmp
      integer i
c-----------------------------------------------------------------------
      open (unit=98,file=filename1,err=2000)
      open (unit=99,file=filename2,err=2001)
      
      do i=1,4
        read (98,'(A)',err=2000) scr
        read (99,'(A)',err=2001) scr
      enddo
      i=1
 100  continue
      read (98,*,end=1000,err=2000) lntqdat(i),lnTrecdat(i)
      read (99,*,end=1000,err=2001) tmp,lndTdTdat(i)
      if (tmp.ne.lntqdat(i)) then
        write (*,*) 'ERROR: Table entries for Tquench do not agree in these files:'
        write (*,*) filename1
        write (*,*) filename2
        stop
      endif
      i=i+1
      goto 100
 1000 continue
      close(98)
      close(99)
      nqdat=i-1
      klo=1
      khi=nqdat
      return
 2000 continue
      write (*,*) 'dsrdreaddof: error while reading ',filename1
      stop
 2001 continue
      write (*,*) 'dsrdreaddof: error while reading ',filename2
      stop
      end

