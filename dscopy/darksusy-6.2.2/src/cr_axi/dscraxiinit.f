      subroutine dscraxiinit
************************************************************************
*** initialization routine for cosmic ray propagation routine
*** the subroutine contains settings of options or control parameters 
*** that not necessarily needed to be changed when changing propagation
*** parameters, but however may require different optimizations.      
*** a few propagation parameters are also set here, however are
*** parameters that the user would seldomly change 
*** finally a switch on whether you want to treat propagation of
*** antiprotons and antideuterons in the approximation of delta
*** function for the gas disc or finite thickness is set here
***
*** the meaning of options/paramters specified before their setting in
*** the file
***      
************************************************************************
      implicit none
      include 'dscraxicom.h'
ccc
ccc when calling cosmic ray axisymmetric propagation routines,
ccc confinement times or green functions are tabulated at given position
ccc in the galaxy this is associated to a label, also used to generate
ccc the file names where tables are stored; if the corresponding label
ccc does not exist already in label file 'axirzlabel.dat', the two
ccc specifications below make the generation of the label automatic out
ccc of the suffix rzaxisuf
ccc at the moment the code does not allow for a non-automatic
ccc generation of this label, so setting rzaxiaddlab=.false. works only
ccc if the user provides also an appropriate setting routine. rzaxisuf
ccc can instead be initialized to something else overwriting the setting
ccc below      
ccc      
      rzaxiaddlab=.true.
      rzaxisuf='ds'     ! maximum length = 4 characters 
ccc
ccc tabulations depend also on cosmic ray propagation parameters. set 
ccc of parameters are also tracked via a label:
ccc epcraxitag and pbcraxitag are used for the tabulation, respectively,
ccc of the positron green function and of the antiproton and
ccc antideuteron confinement time. such label can be specified by the
ccc user or automatically generated, see the structure of the file
ccc dscraxiset_model. the switchs below allows for automatic generator      
ccc in case this names of labels needs to be tracked set 'addlab'
ccc variables to .false. and make sure you are linking to existing
ccc labels and/or you are properly providing new ones. the suffixes
ccc epcraxisuf & pbcraxisuf can be changed by overwritingthe setting
ccc below after the call to dsinit      
ccc      
      epcraxiaddlab=.true.
      pbcraxiaddlab=.true.
      epcraxisuf='epau' ! maximum length = 4 characters 
      pbcraxisuf='pbau' ! maximum length = 4 characters 
ccc
ccc the computation of the pb and db diffusion times is performed summing
ccc a series of zeros for the Bessel function J0; for singular dark matter
ccc profiles or if the position at which the diffusion time is computed is
ccc close to R=0 the series may converge slowly. Parameters checking whether
ccc convergence has been reached are set below; note that too strict values
ccc may significantly increase the running times, hence if lower accuracy
ccc is acceptable it could be preferable to relax them
ccc
      pbprec=1.d-3 ! value setting the truncation of the sum 
c      pbnprec=80 ! number of zeros within which the convergence is
                 ! checked in the first run, otherwise a central
                 ! cylinder is cut out of the diffusion region
                 ! the choice pbnprec=50 is slightly less conservative
      pbnprec=50
ccc
ccc here are propagation parameters (pb & db) which the user can change but
ccc that is not recommended. they can be changed by overwriting on the 
ccc setting below after the call to dsinit    
ccc      
      diffhg=0.1d0    ! half thickness of the disk in kpc
      diffng=1.d0     ! density of hydrogen in the disc in cm^-3
      diffnh=0.d0     ! density of hydrogen in the halo in cm^-3
ccc
ccc radial cores for the source function implemented in view of the 
ccc fact that it is unlikely that the approximation of spatially 
ccc constant diffusion coefficient or positron energy loss term holds 
ccc in the innest part of the galaxy. In practice these matter only for
ccc very singular source functions. they can be changed by overwriting on the 
ccc setting below after the call to dsinit
ccc
      diffrcpb=0.01d0    ! inner radial cut in pb and db routines in kpc
      diffrcep=0.01d0    ! inner radial cut in ep routines in kpc    
ccc
ccc option to link to antiproton and antideuteron propagation libraries
ccc treating the gas disc as a delta function or finite thickness. results
ccc are marginally different in the two cases, however higher accuracy
ccc requires the finite thickness option, increasing the running time
ccc the option can be changed by overwriting on the setting below after
ccc the call to dsinit
ccc
      isopt=1 ! delta function
c      isopt=2 ! finite thickness
ccc      
      return
      end
