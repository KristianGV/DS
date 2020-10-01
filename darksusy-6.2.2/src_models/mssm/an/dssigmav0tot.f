**********************************************************************
*** function dssigmav0tot returns the *total* annihilation cross section
*** sigma v at p=0 for neutralino-neutralino annihilation.
*** This is obtained by summing over all implemented 2-body channels
*** (as returned by dssigmav) plus contributions from final states with
*** more particles.
***                                                                         
***  type : interface                                                       
***                                                                         
*** Units of returned cross section: cm^3 s^-1
***
*** author: Torsten.Bringmann.fys.uio.no
*** date: 2014-11-14
**********************************************************************

      real*8 function dssigmav0tot()
      implicit none
      include 'dsanyieldmodelcom.h'
      include 'dssvcom.h'
      
      real*8 res, dssigmav0
      integer ch

      integer dsidnumber
      integer idold
      data idold/-123456789/
      save idold

      real*8 svsave
      data svsave/-333.d33/
      save svsave
      
c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('MSSM','dssigmav0tot')

      res=0d0
      if (idold.ne.dsidnumber()) then
        do ch=1,numanch2b
          res= res + dssigmav0(anch_2body(ch,1),anch_2body(ch,2))
        enddo


c... Add here corrections to the total annihilation rate from final states
c... with at least 3 particles

c... the above call to dssigmav0 has ensured that we have also calculated
c... the qqg contribution, hence we can simply add it here explicitly
        res = res + sigv(27)    

        svsave=res 
        idold=dsidnumber()
      endif
      
      dssigmav0tot=svsave

      return

      end
