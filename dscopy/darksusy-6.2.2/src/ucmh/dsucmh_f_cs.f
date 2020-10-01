!Function for calculating the theoretical differential fraction of DM in UCMHs from cosmic string loops  
!
!Authors: Madeleine Anthonisen
!          maddyanthonisen@hotmail.com
!         Pat Scott
!          pscott@imperial.ac.uk
!Date: 2014/15
!
!Inputs:  r_cs               cosmic string loop radius (kpc)
!         Gmu                cosmic string tension (dimensionless) 
!         z_c                latest allowed redshift of UCMH collapse (dimensionless) 
!         zstop              redshift at which UCMHs stop growing (dimensionless)
!         alpha              ratio of loop length to horizon scale at formation (dimensionless)
!         beta               ratio of loop length to loop radius (dimensionless)
!         gamma              cosmic string loop decay constant (dimensionless)
!         N                  number of loops formed per Hubble 4-volume (dimensionless)
!         K                  number of loop radii loops can travel before UCMH formation is impossible (dimensionless)
!         sigma              stdev of loop velocity distribution / speed of light (dimensionless)
!
!Output:  dsucmh_f_cs        differential fraction of DM in UCMHs (M_Sun^-1)

      double precision function dsucmh_f_cs(r_cs, Gmu, z_c, zstop, alpha, beta, gamma, N, K, sigma)

      implicit none

      include 'dsucmh.h'
      include 'dsmpconst.h'

      double precision, intent(IN) :: r_cs, Gmu, z_c, zstop, alpha, beta, gamma, N, K, sigma
      double precision :: C, const, constRD, constMD, xstop, x_c, xiRD, xiMD, xdRD, xdMD, logPart, S, v_i
      double precision, parameter :: ct_eq = 1.d3 * c_light * t_eq  !c * teq in m

      xstop = (zeq+1.d0)/(zstop+1.d0)   !x corresponding to zstop
      x_c   = (zeq+1.d0)/(z_c+1.d0)     !x corresponding to zcollapse      

      !Loop decay and formation times in RD era and MD era
      xiRD = (beta*r_cs*mperkpc / (alpha*ct_eq))**(1.d0/2.d0) 
      xiMD = (beta*r_cs*mperkpc / (alpha*ct_eq))**(2.d0/3.d0) 
      xdRD = (beta*r_cs*mperkpc / (Gmu*gamma*ct_eq))**(1.d0/2.d0) 
      xdMD = (beta*r_cs*mperkpc / (Gmu*gamma*ct_eq))**(2.d0/3.d0) 

      !Cosmic string velocity effects
      logPart = log(alpha/(gamma*Gmu))
      v_i = K*(alpha*gamma*Gmu)**(1.d0/2.d0) / (logPart*beta)  !maximum initial velocity
      S = (v_i/sigma)**(3.d0)*(2.0/pi)**(0.5)/3.0              !velocity correction factor (approximated)

      constRD = (alpha*ct_eq/(beta*r_cs*mperkpc))**(1.d0/2.d0)
      constMD = 1.d0

      !Determine which regime CS loop seed belongs to    
      if (xiRD .le. 1.d0 .and. xdRD .le. 1.d0) then
        C=(2+3*xdRD) / (3+3*xdRD)
        const =  constRD
      endif 

      if (xiRD .le. 1.d0 .and. xdMD .gt. 1.d0) then
        C = 1.d0
        const = constRD
      endif

      if (xiMD .gt. 1.d0 .and. xdMD .le. x_c)then
         C = (6*xstop*(xdMD-xiMD)/(xdMD*xiMD)-9)/(2*xstop*(xdMD-xiMD)
     &     /(xdMD*xiMD)-9)
        const = constMD
      endif

      if (xiMD .gt. 1.d0 .and. xdMD .gt. x_c)then
        C = (6*xstop-15*xiMD)/(2*xstop-15*xiMD)
        const = constMD
      endif

      !Differential mass fraction
      dsucmh_f_cs = 1.d-12*S*C*const*16.d0*pi*N*alpha**(2.d0)
     &  *GNewton*SolarMass/
     & (3.d0*dmfrac*kappa*(beta*c_light)**(2.d0)*r_cs*mperkpc)

      end function dsucmh_f_cs
