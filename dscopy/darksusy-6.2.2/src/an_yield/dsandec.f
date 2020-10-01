*****************************************************************
*** suboutine dsandec decomposes yieldcode yieldk to flyxtype
*** fltype and fi
*****************************************************************

      subroutine dsandec(yieldk,fltyp,fi)
      implicit none

      integer yieldk,fltyp,fi

      if (yieldk.lt.100) then
        fltyp=1      ! integrated yields
        fi=yieldk-50
      elseif (yieldk.ge.100.and.yieldk.lt.200) then
        fltyp=2      ! differential yields
        fi=yieldk-150
      elseif (yieldk.ge.1000.and.yieldk.lt.1100) then
         fltyp=3                ! error on integrated yields
         fi=yieldk-1050
      elseif (yieldk.ge.1100.and.yieldk.lt.1200) then
         fltyp=4                ! error on differential yields
         fi=yieldk-1150
      else
         write(*,*) 'DS ERROR on dsandec: ',
     &     'invalid yieldk = ',yieldk
      endif

      if (fi.eq.9) fi=4 ! use pbar for dbar fluxes
      
      if (fi.lt.1.or.fi.eq.10.or.(fi.gt.13.and.fi.lt.21)
     &  .or.fi.gt.23) then
         write(*,*) 'DS error in dsandec: ',
     &        'inconsistent yield code: ',yieldk
        stop
      endif


      end
