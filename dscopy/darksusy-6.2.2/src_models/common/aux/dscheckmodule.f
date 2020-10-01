*****************************************************************************
*** subroutine dscheckmodule checks whether a routine or function is called 
*** for the correct particle module. This internal consistency check should 
*** (at least) be included for all required functions.
***
*** author: Torsten.Bringmann.fys.uio.no
*** date: 2015-06-11
*****************************************************************************

      subroutine dscheckmodule(mymodule,functionname)
      implicit none
      include 'dsidtag.h'
      character*(*) functionname,mymodule
      character*20 testmodule

      testmodule=mymodule

      if (testmodule.ne.moduletag) then
        write(*,*)
        write(*,*) '===================' 
        write(*,*) '!!! FATAL ERROR !!!' 
        write(*,*) '===================' 
        write(*,*)
        write(*,*) 'You have linked to particle module ''', moduletag, 
     &             ''' but ', functionname,' has been called from module ',
     &              testmodule, '.'
        write(*,*)
    
        stop
      endif

      return
      end
