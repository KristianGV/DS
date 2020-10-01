      subroutine dsdmdinit
c_______________________________________________________________________
c subroutine adding to the database of currently available halo models
c a set of default profiles which have been defined as 'defaults' in
c dsdmddriver
c_______________________________________________________________________
      implicit none
      include 'dsdmdcom.h'
      call dsdmddriver(idmdinit,ndummy,redummy,ndummy,chdummy,labdummy,
     &  redummy2)
      return
      end
