      subroutine dsdmdset_halomodel(labin,nin,chin,rein)
c_______________________________________________________________________
c subroutine adding to the repository of currently available halo models
c the one corresponding to the label 'labin'.
c at this stage you need to know the rules that have been defined in
c dsdmddriver on how to link 'labin' to the corresponding object and
c profile as well as what is the list of input parameters and relative
c tags dsdmddriver expects. if such rules are not respected an error
c message is printed and the program stops.
***
***   type : commonly used
***   desc : add an existing halo model to halo repository
***
c_______________________________________________________________________
      implicit none
      include 'dsdmdcom.h'
      character(*) labin
      integer nin
      real*8 rein(nin)
      character*10 chin(nin)
      call dsdmddriver(idmdcreate,nin,rein,nin,chin,labin,redummy2)
      return
      end
