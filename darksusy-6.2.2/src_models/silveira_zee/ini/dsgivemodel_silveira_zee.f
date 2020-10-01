*******************************************************************************
***  subroutine dsgivemodel_silveira_zee sets the model parameters          ***
***                                                                         ***
***  input:                                                                 ***
***                                                                         ***
***    lambda - quartic coupling                                            ***
***    ms     - singlet Higgs mass in GeV                                   ***
***                                                                         ***
*** author: Paolo Gondolo                                                   ***
*** date 2016-05-19                                                         ***
*******************************************************************************
      subroutine dsgivemodel_silveira_zee(mylambda,myms)
      implicit none
      include 'dssilveira_zee.h'

      real*8 mylambda,myms

c-----------------------------------------------------------------------

      mhs=myms
      lambda=mylambda
      mass(khs)=mhs

      return
      end


