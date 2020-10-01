
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
c   NB: See the vdSIDM module for examples on how to set up this routine.
c
c   NB: This is just a DUMMY routine in the empty module.    
c
c  author: Torsten Bringmann (torsten.bringmann@fys.uio.no), 2018-05-31
c=======================================================================
      real*8 function dsrddofDS(Td)

      real*8 Td
      
      write(*,*) 'WARNING: you are calling routine dsrddofDS in',
     &     ' module empty.'
      write(*,*) 'You need to modify this routine.'
      
      dsrddofDS = 1.d0
 
      return
      end
