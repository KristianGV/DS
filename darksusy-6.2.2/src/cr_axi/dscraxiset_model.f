      subroutine dscraxiset_model(gparin,kparin,eparin,labpbin,labepin)
************************************************************************
*** subroutine setting propagation parameters.
***
***   type : REPLACEABLE
***      
*** it assumes that the diffusion coefficient dskdiff and the mean
*** positron energy loss dseppdotmmean are in the hardcoded versions
*** given in the cr_axi directory. if that are replaced you should
*** change accordingly this subroutine:
***
*** input:      
*** - gparin -> vector of propagation parameters for antiprotons,
**    antideuterons and positrons (for the latter only diffhh
***      gparin(1) = diffhh -> 1/2 of vertical size of the diff. region
***              in kpc
***      gparin(2) = diffRh -> radial size of the diff. region in kpc
***      gparin(3) = diffcvel -> galactic wind velocity in km/s
*** - kparin -> vector of parameters for the diffusion coefficient as
***   given in the function dskdiff -- dskdiff can be replaced by some
***   other function and then you need to change this routine      
***      kparin(1) = dble(nkdiff) -> index for functional form hardcoded
***               in dskdiff
***      kparin(2) = k0halo, k0gasdisc -> normalization of the diffusion
***      coefficient in the halo and the disc in 10^27 cm^2 s^-1 
***      kparin(3) = kdiffrig0 -> rigidity break in GV
***      kparin(4) = kdiffdelta ->  spectral index (above kdiffrig0 if
***              nkdiff=2)   
***      kparin(5) = kdiffeta -> exponent of beta**kdiffeta prefactor, 
***              in case kdiffeta is larger than 1.d-7
***      kparin(6) = kdiffdeltalow ->  spectral index below kdiffrig0 
***              if nkdiff=2    
*** - eparin -> vector of parameters for the positron energy loss
***   function dseppdotmmean (ivopt=1) or dseppdotmana (ivopt=1).
***   dseppdotmmean can be replaced by some other function and then you
***   need to change this routine, dseppdotmana cannot be replaced
***   since it is the only form for which the analytic matching pp<->v
***   works. ivopt=2 also works only in case nkdiff = 1 or 4       
***      eparin(1) = dble(ivopt) -> option for energy losses and tabulation
***           of the positron green function in the variable v 
***      ivopt=1 -> the matching between  pp and the variable v is done 
***           assuming fully general momentum loss rate and spatial
***           diffusion coefficient, with a numerical integral + a
***           tabulation involved; you link to dskdiff and dseppdotmmean
***           which can be replaced by arbitry functions of, respectively
***           rigidity and energy 
***      ivopt=2 -> the matching between  pp and the variable v is done 
***           assuming the momentum loss rate scales with pp^2, and the
***           diffusion coefficient has a form which allows for an
***           analitical inversion between v and pp, see Eq.~(6) & (7) in
***           Baltz & Edsjo (1998)
***      NOTE: ivopt is a common block variable on a fixed range, user
***      cannot other ivopt values!     
***      if ivopt = 1 
***        eparin(2) = starlightmean -> mean optical + IR light density,
***           eV cm^-3
***        eparin(3) = bfieldmean -> mean magnetic field, \muG
***      if ivopt = 2 
***        eparin(2) = blossmean -> some mean value pdot0 assuming the
***           momentum loss rate scales with pdot0 * pp^2
*** - labpbin = label for pb & db confinement time; if set to 'none'
***   a label is searched in the label file and automatically generated
***   if it does not exist -- default option as set in dscraxiinit, if  
***   that option has been changed then output modified accordingly
*** - labepin = label for ep green functions; if set to 'none'
***   a label is searched in the label file and automatically generated
***   if it does not exist -- default option as set in dscraxiinit, if  
***   that option has been changed then output modified accordingly
***      
************************************************************************
      implicit none
      include 'dscraxicom.h'
      real*8 gparin(3),kparin(6),eparin(3)
      character*10 labpbin,labepin,labelout
      logical ivck
      real*8 dspbcraxisetlabel,dsepcraxisetlabel
      external dspbcraxisetlabel,dsepcraxisetlabel
ccc
ccc propagation parameters (pb, db & ep):
ccc
      diffhh=gparin(1)     ! 1/2 of vertical size of the diff. region in kpc
ccc
ccc other propagation parameters (pb & db only):
ccc
      diffRh=gparin(2)    ! radial size of the diff. region in kpc
      diffcvel=gparin(3)   ! galactic wind velocity in km/s
ccc
ccc diffusion coefficient parameters:
ccc
      nkdiff=int(kparin(1)+0.001d0) ! functional form hardcoded in dskdiff
      k0halo=kparin(2) ! normalization of the diffusion coefficient 
                      ! in the halo in 10^27 cm^2 s^-1k0gasdisc
      k0gasdisc=kparin(2)  ! normalization of the diffusion coefficient 
                          ! in the disc in 10^27 cm^2 s^-1
      kdiffrig0=kparin(3)  ! rigidity break in GV
      kdiffdelta=kparin(4)  ! spectral index (above the break if nkdiff=2)
      kdiffeta=kparin(5)    ! multiply the diffusion coefficient by beta
      if(nkdiff.eq.2) then   ! diffusion coefficient with a break
        kdiffdeltalow=kparin(6) ! spectral index below the break
      endif
ccc
ccc read/store label for pb & db confinement time
ccc      
      if(labpbin(1:4).eq.'none') then
ccc
ccc check whether a label exists, if not generate it automatically
ccc        
        call dsgetlabel2(dspbcraxisetlabel,labelout)
      else
ccc
ccc add the label to labelfile in case it does not exist already; if it 
ccc exists but corresponds to another set of parameters the program stops
ccc        
        call dsaddlabel2(dspbcraxisetlabel,labpbin)
      endif  
ccc
ccc call the subroutine which ckecks whether any of the parameters
ccc entering in the pb & db confinement time has actually been changed
ccc in case it has it increments the integer dspbcraxino by one unit
ccc
ccc NB: USER MUST DO THIS CALL to get proper linking of confinement time
ccc tables, implement this check in any other version of this replaceble
ccc      
      call dspbcraxick
ccc
ccc read/store label for ep green functions:
ccc      
      if(labepin(1:4).eq.'none') then
ccc
ccc check whether a label exists, if not generate it automatically
ccc        
        call dsgetlabel2(dsepcraxisetlabel,labelout)
      else
ccc
ccc add the label to labelfile in case it does not exist already; if it 
ccc exists but corresponds to another set of parameters the program stops
ccc        
        call dsaddlabel2(dsepcraxisetlabel,labepin)
      endif
ccc
ccc call the subroutine which ckecks whether any of the parameters
ccc entering in the ep green function has actually been changed
ccc in case it has it increments the integer dsepcraxino by one unit
ccc
ccc NB: USER MUST DO THIS CALL to get proper linking of green function
ccc tables, implement this check in any other version of this replaceble
ccc      
      call dsepcraxick
ccc
      ivopt=int(eparin(1)+0.001d0)  ! 1 or 2 pp to v matching
ccc
      if(ivopt.eq.1) then
ccc
ccc some mean from moskalenko talk or porter paper
ccc
        starlightmean=eparin(2) ! optical + IR, eV cm^-3
ccc
ccc some 'local' mean in 8*dexp(-(r-r0)/50kpc) *dexp(-|z|/3kpc)
ccc random + regular, strong talk 
ccc
        bfieldmean=eparin(3)      ! \muG
      elseif(ivopt.eq.2) then
ccc
ccc some coefficient in from of pp^2
ccc
        blossmean=eparin(2) 
ccc
ccc check that, in case ivopt=2, you are using a compatible diffusion 
ccc coefficient. note that this check refers to the hard-coded set of
ccc diffusion coefficients. 
ccc
        ivck=.true.
        if((nkdiff.eq.1).and.(dabs(kdiffeta).le.1.d-7)) ivck=.false.
        if((nkdiff.eq.4).and.(dabs(kdiffeta).le.1.d-7)) ivck=.false.
        if(ivck) then
          write(*,*) 'DS: in  dscraxiset_model ivopt = ',ivopt
          write(*,*) 'DS: but with wrong nkdiff, kdiffeta = ',nkdiff, 
     &               kdiffeta
          write(*,*) 'DS: program stops'
          stop
        endif
      else
        write(*,*) 'DS: in  dscraxiset_model ivopt = ',ivopt
        write(*,*) 'DS: out of the allowed range 1-2 '
        write(*,*) 'DS: program stops'
        stop
      endif    
      return
      end
