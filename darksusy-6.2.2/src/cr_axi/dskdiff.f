************************************************************************
*** Function dskdiff returns the spatial diffusion coefficient
***
***  type : REPLACEABLE                                                   
***                                                                       
*** inputs:
***     rig - particle rigidity in GV
***     n = 1 - diffusion coefficient in the halo
***     n = 2 - diffusion coefficient in the gas disc
***     beta - particle velocity in c units (this value affects
***            dskdiff only in case betalabel=.true.)
***
*** additional inputs are through common blocks, in particular the 
*** spectral is set by nkdiff:
***     nkdiff = 1 - K0 * (rig/kdiffrig0)**kdiffdelta
***     nkdiff = 2 - K0 * (rig/kdiffrig0)**kdiffdelta 
***                    iff rig>kdiffrig0
***                  K0 * (rig/kdiffrig0)**kdiffdeltalow 
***                    iff rig<kdiffrig0
***     nkdiff = 3 - K0 * (1+rig/kdiffrig0)**kdiffdelta 
***     nkdiff = 4 - K0 * [1+(rig/kdiffrig0)**kdiffdelta] 
*** K0 is eqal to k0halo or k0gasdisc depending on n. For 
*** dabs(kdiffeta).gt.1.d-7 all the above are multiplied by 
***   beta**kdiffeta 
***
*** output: 10^27 cm^2 s^-1
***
***  author: Piero Ullio
***  modified: Torsten Bringmann, 11/06/2015 (adapted to replaceable   
***                                           function concept)
************************************************************************
      real*8 function dskdiff(rig,n,beta)
      implicit none
      include 'dscraxicom.h'
      integer n
      real*8 rig,beta
ccc
      if(nkdiff.eq.1) then
ccc
ccc single power law
ccc
        dskdiff=(rig/kdiffrig0)**kdiffdelta
      elseif(nkdiff.eq.2) then
ccc
ccc double power law with a break at the rigidity kdiffrig0 
ccc
        if(rig.lt.kdiffrig0) then
          dskdiff=(rig/kdiffrig0)**kdiffdeltalow
        else
          dskdiff=(rig/kdiffrig0)**kdiffdelta
        endif
      elseif(nkdiff.eq.3) then
ccc
ccc power law with low energy cutoff 
ccc
        dskdiff=(1.d0+rig/kdiffrig0)**kdiffdelta
      elseif(nkdiff.eq.4) then
ccc
ccc power law with low energy cutoff, second version 
ccc
        dskdiff=1.d0+(rig/kdiffrig0)**kdiffdelta
      else
        write(*,*) 'DS: dskdiff called with invalid nkdiff = ',
     &              nkdiff
        write(*,*) 'DS: program stopped'
        stop
      endif
      if(n.eq.1) then
ccc
ccc dskdiff equal to the spatial diffusion coefficient in the halo
ccc
        dskdiff=dskdiff*k0halo
      else if(n.eq.2) then
ccc
ccc dskdiff equal to the spatial diffusion coefficient in the gas 
ccc disc
ccc
        dskdiff=dskdiff*k0gasdisc
      else
        write(*,*) 'DS: dskdiff called with invalid n = ',n
        write(*,*) 'DS: program stopped'
        stop
      endif
ccc
ccc multiply by beta if dabs(kdiffeta).gt.1.d-7
ccc
      if(dabs(kdiffeta).gt.1.d-7) then
        dskdiff=dskdiff*beta**kdiffeta
      endif 
      return
      end
