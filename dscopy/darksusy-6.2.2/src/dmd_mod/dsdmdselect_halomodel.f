      subroutine dsdmdselect_halomodel(labin)
c_______________________________________________________________________
c subroutine selecting from the database of currently set halo models
c the one corresponding to the label 'labin'.
c if labin refers to a temporary halo model the last temporary setting
c is returned.      
c this subroutine also make sure that the following quantities, needed
c for detection rates are set by dsdmddriver consistently with 'labin': 
c   - dmdmw: logical variable saying whether the model should be used 
c     for Milky Way rates or not
c   - dmdlabel: internal dmd label for tabulation files      
c   - dmdobjdist: distance from the center of dark matter object (the 
c     Sun galactocentric distance in case of the Milky Way)
c   - dmdradintr: inner density truncation radius
c   - dmdradouttr: outer radius at which the density is set to zero
c   - dmdrho0: if dmdmw=.true. this is the local halo density      
***
***   type : commonly used
***   desc : select between halo models in the halo repository
***
c_______________________________________________________________________
      implicit none
      include 'dsdmdcom.h'
      character(*) labin
      call dsdmddriver(idmdload,ndummy,redummy,ndummy,chdummy,labin
     &  ,redummy2)
      return
      end
