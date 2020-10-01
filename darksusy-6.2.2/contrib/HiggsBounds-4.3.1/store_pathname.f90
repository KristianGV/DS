!******************************************************************
module store_pathname
!******************************************************************
 implicit none

 integer,parameter:: pathname_length= 56
 character(len=pathname_length),parameter :: pathname= &
     &     "/home/kriss/DS/darksusy-6.2.2/contrib/HiggsBounds-4.3.1" // &
     &     "/"

end module store_pathname
!******************************************************************
