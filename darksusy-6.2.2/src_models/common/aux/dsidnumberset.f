	logical function dsidnumberset()
************************************************************************
*** The function dsidnumber returns an integer number that is 
*** unique for each new particle physics model. This function, dsidnumberset,
*** works as a check if the function dsidnumber is properly initialized in the
*** respective particle module. If it returns false, the output of dsidnumber 
*** should not be used
***
***    Type: interface function
***	
*** Author: Torsten Bringmann, torsten.bringmann@fys.uio.no
*** Date: 10/10/20174	
************************************************************************
      implicit none
      include 'dsidnumber.h'

      if (checksum.eq.31415) then
        dsidnumberset=.true.
      else
        dsidnumberset=.false.      
      endif

      return
      end
