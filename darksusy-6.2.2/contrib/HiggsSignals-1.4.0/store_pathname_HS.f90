!******************************************************************
module store_pathname_HS
!******************************************************************
 implicit none

 integer,parameter:: pathname_length= 57
 character(len=pathname_length),parameter :: pathname_HS= &
     &     "/home/kriss/DS/darksusy-6.2.2/contrib/HiggsSignals-1.4.0" // &
     &     "/"

end module store_pathname_HS
!******************************************************************
