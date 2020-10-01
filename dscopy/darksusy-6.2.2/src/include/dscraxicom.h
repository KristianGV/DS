*         -*- mode: fortran -*-

ccc
ccc propagation model label for the tabulation of the positron green 
ccc function and of the antiproton "confinement time" and number of
ccc characters defining crproptag
ccc
      character*10 epcraxitag,pbcraxitag
ccc
ccc logical variable to (dis)allow to modify the repository label
ccc files adding new automatically generated labels
ccc 
      logical epcraxiaddlab,pbcraxiaddlab
ccc
ccc propagation model suffix label in case you allow for automatic 
ccc genereration of propagation model label for new labels
ccc
      character*4 epcraxisuf,pbcraxisuf
      common/propaxitagcom/epcraxiaddlab,pbcraxiaddlab,
     &  epcraxitag,pbcraxitag,epcraxisuf,pbcraxisuf
      
ccc
ccc position label, logical variable and suffix label for the (R,z) 
ccc cylindrical coordinates, analogously as above
ccc
      character*10 rzaxitag
      logical rzaxiaddlab
      character*4 rzaxisuf
      common/posaxitagcom/rzaxiaddlab,rzaxitag,rzaxisuf

ccc
ccc nparkdiff is the number of parameters needed to define the diffusion 
ccc coefficient and which needs to be stored in label file; nparkdiff 
ccc needs to be lower or equal to nparkdiffdef, see below.
ccc
      integer nparkdiff 
      common/kdiffdefcom/nparkdiff

ccc
ccc spatial diffusion coefficient parameters, for the DS hardcoded diffusion
ccc coefficient, see the header of the function dskdiff; this function
ccc is defined by 7 parameters and nparkdiff above is set to nparkdiffdef  
ccc
      real*8 k0halo,k0gasdisc,kdiffrig0,kdiffdelta,kdiffdeltalow,
     &  kdiffeta
      integer nkdiff
      common/kdiffparcom/k0halo,k0gasdisc,kdiffrig0,kdiffdelta,
     &  kdiffdeltalow,kdiffeta,nkdiff
      integer nparkdiffdef 
      parameter (nparkdiffdef=7)

ccc
ccc parameters defining the size of the axisymmetric diffusion region,
ccc plus other parameters relevant for propagation (except kdiff),
ccc and eventual inner cutoffs for the dark matter density profile
ccc below which it is not reasonable to either extrapolate the source
ccc or add a source while still treating diffusion under the
ccc approximation of spatially constant diffusion coefficient and/or
ccc energy loss rate (mean values, mostly extrapolated from values
ccc close to the observer) 
ccc

      real*8 diffRh, ! radial size of the diff. region in kpc
     &  diffhh,      ! 1/2 of vertical size of the diff. region in kpc
     &  diffhg,      ! 1/2 of the gas disc height in kpc
     &  diffnh,      ! hydrogen density in the halo in cm^-3
     &  diffng,      ! hydrogen density in the disk in cm^-3
     &  diffcvel,    ! spatially constant convective velocity in km s^-1
     &  diffrcpb,    ! inner radial cut in pb and db routines in kpc    
     &  diffrcep     ! inner radial cut in ep routines in kpc

      common/diffpropparcom/diffRh,diffhh,diffhg,diffnh,diffng,diffcvel,
     &  diffrcpb,diffrcep

ccc
ccc electron momentum loss parameters, see the header of the function
ccc pdotminus and of functions called therein
ccc
      real*8 starlightmean,bfieldmean,denHImean,denHIImean,denH2mean,
     & blossmean
      integer ivopt
      common/pdotmeancom/starlightmean,bfieldmean,denHImean,denHIImean,
     &  denH2mean,blossmean,ivopt

ccc
ccc antiproton and antideutron parameters involved in the computation 
ccc of the diffusion time, when using the tabulated functions rather
ccc than linking function itself - they are assumed to be the same for
ccc both species
ccc 
      real*8 pbprec
      integer pbnprec
ccc
ccc is option (delta function for the disc / finite thickness)
ccc
      integer isopt
      common/pbtaboptcom/pbprec,pbnprec,isopt

ccc
ccc integer checking whether propagation parameters entering in pb/db
ccc confinement times are changed       
ccc      
      integer dspbcraxino
      common/dspbcraxinocom/dspbcraxino

ccc
ccc integer checking whether propagation parameters entering in ep
ccc green functions are changed       
ccc
      integer dsepcraxino
      common/dsepcraxinocom/dsepcraxino

ccc
ccc parameters setting the default orbit, see the header of the 
ccc subroutine dsorbitps
ccc
      integer horbit ! orbit choice
      real*8 y0in,dminin,vxin,vzin  ! straight line orbit
      real*8 vcircin,radorb,rosb  ! circular orbit
      integer norbit              ! circular orbit
      common/orbitcom/y0in,dminin,vxin,vzin,vcircin,radorb,rosb,horbit,
     &  norbit

      save /propaxitagcom/,/posaxitagcom/,/kdiffdefcom/,/kdiffparcom/,
     &  /diffpropparcom/,/pdotmeancom/,/pbtaboptcom/,/orbitcom/,
     &  /dspbcraxinocom/,/dsepcraxinocom/

