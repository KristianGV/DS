      real*8 function dsdmdrho0mw(labhalo)
************************************************************************
***  Local halo density for the halo model with access label 'labhalo',
***  within the set of models defined in the halo model repository.
***  In case you are linking to a model which has not been defined as
***  suitable for the Milky Way, the program stops.      
************************************************************************
      implicit none
      include 'dsdmdcom.h'  ! to interface dmd halo parameters
      character(*) labhalo
ccc
      call dsdmdselect_halomodel(labhalo)
ccc
ccc check whether you have called this function for a halo model which
ccc has been defined for the Milky Way or not
ccc      
      if(.not.dmdmw) then
        write(*,*) 'DS: call to dsdmdrho0mw with the halo label: ',
     &       labhalo 
        write(*,*) 'DS: which has not been initialized as suitable'
        write(*,*) 'DS: for the Milky Way. Program stopped'
        stop
      endif
      dsdmdrho0mw=dmdrho0
      return
      end
