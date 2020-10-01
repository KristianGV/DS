      real*8 function dslosifunsph(radius,ilosi)
************************************************************************
*** function to be integrated over the line of sight and, eventually, 
*** over solid angle in case of spherically symmetric emissivities,
*** including (eventually) an absorption factor
*** inputs:
***   radius = radial coordinate in kpc
***   ilosi = index for link to dslosidriver
************************************************************************
      implicit none
      include 'dsdmdcom.h' ! for ksoupow
      include 'dsdvcom.h'
      real*8 radius
      integer ilosi
      real*8 out
ccc
      integer powerloc
      common/dslosifunsphcom/powerloc
ccc
      ksoupow=powerloc
      dvsph(idvrsph)=radius
      call dslosidriver(ilosi,ntysph,dvsph,out)
      dslosifunsph=out
      return
      end
