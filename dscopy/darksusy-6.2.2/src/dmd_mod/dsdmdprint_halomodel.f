      subroutine dsdmdprint_halomodel(labin)
c_______________________________________________________________________
c subroutine printing name and parameters for the halo with label equal
c to 'labin' within the database of currently set halo models
c if labin='printall' all models in the database are printed      
c if labin refers to a temporary halo model the last temporary setting
c is returned.      
c_______________________________________________________________________
      implicit none
      include 'dsdmdcom.h'
      character(*) labin
      call dsdmddriver(idmdparprint,ndummy,redummy,ndummy,chdummy,labin
     &  ,redummy2)
      return
      end
