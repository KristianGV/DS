*******************************************************************************
***  subroutine dsmodelsetup computes various particle physics quantities   ***
***                                                                         ***
***  output :                                                               ***
***    ierr  - error, should not continue if ierr<>0                        ***
***    iwarn - warning, can continue despite iwarn<>0                       ***
***                                                                         ***
***                                                                         ***
*** author: Paolo Gondolo                                                   ***
*** date 2016-05-19                                                         ***
*** Modified: Joakim Edsjo, 2019-12-13 (added SM width, fixed typo)         ***
*******************************************************************************
      subroutine dsmodelsetup(ierr,iwarn)
      implicit none
      include 'dssilveira_zee.h'
      integer ierr,iwarn

      real*8 higgs_inv_width,exp_higgs_inv_width,aux
      real*8 dsgammah

c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('silveira_zee','dsmodelsetup')

      call dsnewidnumber ! This should always be the first call in dsmodelsetup,
                         ! and makes sure the model has a new unique ID number.
                         ! Several routines in src_models/ and in src/ 
                         ! *require* this to be set anew for every new model.


      exp_higgs_inv_width=0.19d0/(1.d0-0.19d0)*4.07d-3 ! Belanger et al 1306.2941

      ierr=0
      iwarn=0
      
c-----------------------------------------------------------------------

      aux=1.d0-4.d0*(mass(khs)/mass(khsm))**2
      if (aux.le.0.d0) then
        higgs_inv_width = 0.d0
      else
        higgs_inv_width = (lambda*v0)**2/(32.d0*pi*mass(khsm))*sqrt(aux)
        if (higgs_inv_width.gt.exp_higgs_inv_width) then
          iwarn=1
        endif
      endif

c      write (*,*) 'PG dsmodelsetup > ',higgs_inv_width,exp_higgs_inv_width
c      write (*,*) 'PG dsmodelsetup > ',ierr,iwarn


c...We need to set width of SM Higgs including decay H -> S + S
      width(khsm)=dsgammah(mass(khsm))
      
      return
      end
