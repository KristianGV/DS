      real*8 function dseppdotm(pp,iv)
************************************************************************
*** electron/positron momentum loss rate -pdot(pp) -- average value in
*** the milky way, or value applying close to the place where the 
*** equilibrium number density for electron/positrons is computed.
*** input: pp - momentum (GeV)
***        iv = 1 - links to a function including a full set of 
***                 energy loss effects
***        iv = 2 - links to a function with the simple scaling
***                 blossmean *pp^2
*** output: 10^-16 GeV s^-1
************************************************************************
      real*8 pp
      integer iv
      real*8 dseppdotmmean,dseppdotmana
ccc
      if(iv.eq.1) then
        dseppdotm=dseppdotmmean(pp) ! 10^-16 GeV s^-1
      elseif(iv.eq.2) then
        dseppdotm=dseppdotmana(pp) ! 10^-16 GeV s^-1
      else
        write(*,*) 'DS: in dseppdotm invalid iv = ',iv
        write(*,*) 'DS: program stopped'
        stop
      endif
      return
      end


      real*8 function dseppdotmana(pp)
************************************************************************
*** electron/positron momentum loss rate in the approximate form:
***    -pdot(pp) = bloss *pp^2
*** input: pp - momentum (GeV)
*** output: 10^-16 GeV s^-1
************************************************************************
      implicit none
      include 'dscraxicom.h'
      real*8 pp ! momentum in GeV
ccc
      dseppdotmana=blossmean*pp**2  ! 10^-16 GeV s^-1
      return
      end


      real*8 function dseppdotmmean(pp)
************************************************************************
*** electron/positron momentum loss rate -pdot(pp) -- average value in
*** the milky way, or value applying close to the place where the 
*** equilibrium number density for electron/positrons is computed.
*** input: pp - momentum (GeV)
*** additional inputs through common blocks:
***    starlightmean - mean starlight energy density in ev/cm
***    bfieldmean - mean interstellar magnetic field in \mu G
*** FIX THE OTHERS !!!!!!!!!!!!!!!!!!!!!
*** output:  10^-16 GeV s^-1
************************************************************************
      implicit none
      include 'dscraxicom.h'
      include 'dsmpconst.h' 
      real*8 pp
      real*8 b0_ic,b0_syn,b0_coul,b0_brem
      real*8 b_ic,b_syn,b_brem,b_bremn,b_bremi,b_coul,b_ion,energy,gam
      real*8 dssphcoulombloss,dssphionizloss,dssphbremloss
      real*8 hmass,hionizen,radlength,h2mass,h2ionizen,radlength2

      b0_ic=7.64d-1 ! assuming 1 ev/cm as energy density
      b0_syn=2.54d-2 ! assuming 1 \muG as magnetic field
ccc
ccc inverse compton, CMB + starlight:
ccc
      b_ic=b0_ic*(0.25d0+starlightmean)*pp**2 

ccc
ccc synchrotron:
ccc
      b_syn=b0_syn*bfieldmean**2*pp**2
ccc
ccc temporary fix: skip the others at this stage
ccc
      dseppdotmmean=b_ic+b_syn
      return
ccc
ccc end of temporary fix
ccc
      energy=dsqrt(pp**2+m_e**2)
      gam=energy/m_e
      b0_coul=7.6d-2
      b0_brem=7.0d-1
ccc
ccc hydrogen atom properties
ccc
      hmass=m_p ! proton mass
      hionizen=13.6d-9  ! ionization potential
      radlength=62.8d0*5.60947d23  !GeV/cm^2
ccc
ccc hydrogen molecule properties
ccc
      h2mass=2.d0*hmass ! mass
      h2ionizen=4.52d-9  !ionization potential
      radlength2=61.3d0*5.60947d23  !GeV/cm^2
ccc
ccc interstellar medium densities: fix this later, maybe move to main 
ccc file 
ccc
      denHImean=0.d0
      denHIImean=0.d0
      denH2mean=0.d0
ccc
ccc ionization (on neutral gas)
ccc
      b_ion=b0_coul*(denHImean*dssphionizloss(energy,hmass,hionizen)
     &        +denH2mean*dssphionizloss(energy,h2mass,h2ionizen))
ccc
ccc coulomb (on ionized cold plasma)        
ccc
      b_coul=b0_coul*denHIImean*dssphcoulombloss(energy,denHIImean)
ccc   
ccc bremstrahlung
ccc bremstrahlung in a neutral gas
ccc
      b_bremn=b0_brem*(denHImean*dssphbremloss(energy,hmass,radlength)
     &        +denH2mean*dssphbremloss(energy,h2mass,radlength2))
ccc
ccc good approximation for loss due to ee + ep bremstrahlung in an 
ccc ionized gas
ccc 
      b_bremi=b0_brem*denHIImean*energy*(dlog(gam)+0.36d0)
      b_brem=b_bremn+b_bremi
      dseppdotmmean=dseppdotmmean+b_brem+b_coul+b_ion
      return
      end



      real*8 function dssphionizloss(energy,amass,ionizen)
ccc  dssphionizloss gives the energy dependent part of the electron 
ccc  energy loss rate by ionisation in a neutral gas of atoms with 
ccc  mass amass and ionization potential ionizen (ref: Longair, pag.60)
ccc  input: energy,amass,ionizen in GeV
ccc  output: adimensional 
      implicit none
      include 'dsmpconst.h' 
      real*8 energy,amass,ionizen
      real*8 gamma,beta,emax,gammamin

      gamma=energy/m_e
      beta=dsqrt(gamma**2-1.d0)/gamma
c  the formula below is valid only if e- momentum >> ionizen
      gammamin=1.0001d0
      if(gamma.lt.gammamin) then
        gamma=gammamin
        beta=dsqrt(gamma**2-1.d0)/gamma
      endif
      emax=2.d0*gamma**2*amass**2*m_e*beta**2/(m_e**2+amass**2+2.d0
     &     *gamma*m_e*amass)
      dssphionizloss=(dlog(gamma**2*m_e*beta**2*emax/(2.d0*ionizen**2))
     &                -(2.d0/gamma-1.d0/gamma**2)*dlog(2.d0)+
     &                1.d0/gamma**2+1.d0/8.d0*(1.d0-1.d0/gamma)**2)/beta
      return
      end


      real*8 function dssphcoulombloss(energy,eden)
ccc  dscoulombloss gives the energy dependent part of the Coulomb energy
ccc  loss rate in a fully ionized plasma (in the cold plasma limit)
ccc  with electron number density eden  (ref: Ginzburg, pag. 358-361)
ccc  input: energy in GeV , eden in cm^-3
ccc  output: adimensional 
      implicit none
      include 'dsmpconst.h' 
      real*8 energy,eden
      real*8 gamma,beta,echarge,cmGeV,gammamin

      if(eden.lt.1.d-10) then
        dssphcoulombloss=0.d0
        return
      endif
      gamma=energy/m_e
      beta=dsqrt(gamma**2-1.d0)/gamma
      gammamin=1.0001d0
      if(gamma.lt.gammamin) then
        gamma=gammamin
        beta=dsqrt(gamma**2-1.d0)/gamma
      endif
      cmGeV=5.06765d13
      echarge=0.085425d0
      dssphcoulombloss=(dlog(gamma*m_e**3*beta**2/(4.d0*pi
     &         *eden*echarge**2)*cmGeV**3)-2.d0*beta**2+1.25d0)/beta
      return
      end


      real*8 function dssphbremloss(energy,mass,radl)
ccc  dssphbremloss gives the energy dependent part of the electron 
ccc  energy loss rate by bremstrahlung in a neutral gas
ccc  (ref: Ginzburg, pag.408 and Strong-Moskalenko astro-ph/9807150)
ccc  input: energy,mass in GeV, radlength in GeV/cm^2
ccc  output: adimensional 
      implicit none
      include 'dsmpconst.h' 
      real*8 energy,mass,radl
      real*8 gam,br1,br2,acoeff,bcoeff

      gam=energy/m_e
      if(gam.lt.100.d0) then
         dssphbremloss=energy*(dlog(gam)+0.36d0)
      elseif(gam.gt.800.d0) then
        dssphbremloss=c_light*1.d5*energy*mass/radl/7.0d-17
      else
        br1=100.d0*m_e*(dlog(100.d0)+0.36d0)
        br2=c_light*1.d5*800.d0*m_e*mass/radl/7.0d-17
        acoeff=(br2-br1)/(800.d0-100.d0)
        bcoeff=br1-100.d0*acoeff
        dssphbremloss=acoeff*gam+bcoeff
      endif
      return
      end

