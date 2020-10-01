*******************************************************************************
*** Function dsddTrmax returns the maximal (kinetic) recoil energy Trmax of ***
*** a target particle at rest, when hit by an incoming particle of kinetic  ***
*** energy Tkin. The maximal recoil energy is given in the lab frame (where ***
*** the target particle is initially at rest), with no assumption made      ***
*** concerning whether the incident particle is relativistic or not.        ***
*** For isotropic scattering (in the CMS), the actual recoil energy is      ***
*** equally distributed between 0 and the value returned by this function.  ***
***                                                                         ***
***  Input:                                                                 ***
***    Tkin    - kinetic energy of incident particle [GeV]                  ***
***    min     - mass of incident particle [GeV]                            ***
***    mtarget - mass of target particle [GeV]                              ***
***                                                                         ***
***  Units of output: [GeV]                                                 ***
***                                                                         ***
*** author: Torsten.Bringmann.fys.uio.no                                    ***
*** date 2018-06-22                                                         ***
*******************************************************************************
      real*8 function dsddTrmax(Tkin, min, mtarget)
      implicit none

      real*8 Tkin, min, mtarget
      
      dsddTrmax=0.0d0

      if (mtarget.le.0) return    
      if (Tkin.le.0) return    
      
      dsddTrmax = 2.*mtarget*(Tkin**2 + 2.*min*Tkin) /
     &            ((min+mtarget)**2 + 2.*mtarget*Tkin) 

      return
      end


