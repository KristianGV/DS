
c_______________________________________________________________________
c  Function dsrddofDS returns the effective number of energy degrees of 
c  freedom in the dark sector. For simplicity, it is assumed that all 
c  dark sector particles are in thermal equilibirum with each other. In
c  this approximation, energy and entropy degrees of freedom are identical.
c
c  type : interface                                                   
c
c  desc : Relativistic BSM degrees of freedom 
c
c  input:
c    Td - dark sector temperature [GeV]
c
c   NB: This specific implementation only applies to the vdSIDM module, 
c       but it is of course straightforward to extend to other field
c       contents in the dark sector.
c
c  author: Torsten Bringmann (torsten.bringmann@fys.uio.no), 2018-05-31
c=======================================================================
      real*8 function dsrddofDS(Td)
      implicit none
      include 'dsvdSIDM.h'
      real*8 Td
      
      real*8 res, dsrdsingledof

c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('vdSIDM','dsrddofDS')

      res = 0d0 

c... DM fermion, vector mediator:
      if (DMtype.eq.1.and.Mediatortype.eq.2) then
        
        res = 3.5d0 ! = (7/8)*2*2 : we always have DR, assumed to be Dirac fermions
        res = res + 3.*dsrdsingledof(mass(kmed)/Td,1)                 
        res = res + 4.*dsrdsingledof(mass(kdm)/Td,2)                   

c... DM fermion, scalar mediator
      elseif (DMtype.eq.1.and.Mediatortype.eq.3) then

        res = 3.5d0 ! = (7/8)*2*2 : we always have DR, assumed to be Dirac fermions
        res = res + 1.*dsrdsingledof(mass(kmed)/Td,1)                 
        res = res + 4*dsrdsingledof(mass(kdm)/Td,2)                 
            
      else
        write(*,*) 'ERROR in dsanwx: called with unsupported '
        write(*,*) 'DM/DR/mediator spin combination: ',  
     &              DMtype, Mediatortype, DRtype   
        write(*,*) 'program stopping...'
        stop
      endif
   
      dsrddofDS = res
 
      return
      end
