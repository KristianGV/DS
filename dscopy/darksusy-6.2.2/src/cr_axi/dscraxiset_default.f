      subroutine dscraxiset_default
************************************************************************
*** a sample default setting of propagation parameters. this subroutine 
*** is called in the program inizialization. to change the propagation
*** parameters to another set call the routine dscraxiset_model.
*** users CAN NOT change hardcoded values of paramaters in this routine
*** changing parameters here could possibly give wrong results 
***      
************************************************************************
      implicit none
      include 'dscraxicom.h'
ccc
ccc model "A", setting the dskdiff function (pb, db & ep):
ccc
      nkdiff=2            ! 2nd functional form hardcoded in dskdiff
      k0halo=25.d0        ! normalization of the diffusion coefficient 
                          ! in the halo in 10^27 cm^2 s^-1
      k0gasdisc=k0halo    ! normalization of the diffusion coefficient 
                          ! in the disc in 10^27 cm^2 s^-1
      kdiffrig0=4.d0      ! rigidity break in GV
      kdiffdeltalow=0.d0  ! spectral index below the break
      kdiffdelta=0.6d0    ! spectral index above the break
      kdiffeta=1.d0    ! multiply the diffusion coefficient by beta
ccc
ccc other propagation parameters (pb, db & ep):
ccc
      diffhh=4.d0     ! 1/2 of vertical size of the diff. region in kpc
ccc
ccc other propagation parameters (pb & db only):
ccc
      diffRh=30.d0    ! radial size of the diff. region in kpc
      diffcvel=10.d0   ! galactic wind velocity in km/s
ccc
ccc for this hardcoded set of parameters labels for proper linking of
ccc tables for pb & db confinement times and the ep green function
ccc already exist. you cannot generate a label here since the
ccc setting/checking routines are replaceable functions, while here
ccc only the the default setup is linked. in case this default setting
ccc is mantained, a negative value of the two integer below forces the
ccc initialization of the tables on the first call to table themselves
ccc see dspbtdaxitab, dsdbtdaxitab & dsepgreenaxitab       
ccc       
ccc      
      dspbcraxino=-1
      dsepcraxino=-1

ccc
ccc model "H" for energy loss parameters in the mean field model
ccc
      ivopt=1  ! the matching between  pp and the variable v is done 
               ! assuming fully general momentum loss rate and spatial
               ! diffusion coefficient, with a numerical integral + a
               ! tabulation involved
ccc
ccc some mean from moskalenko talk or porter paper
ccc
      starlightmean=1.d0-0.25d0 ! optical + IR, eV cm^-3
ccc
ccc some 'local' mean in 8*dexp(-(r-r0)/50kpc) *dexp(-|z|/3kpc)
ccc random + regular, strong talk 
ccc
      bfieldmean=6.d0      ! \muG
ccc
      return
      end
