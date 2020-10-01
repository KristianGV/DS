!Subroutine for calculating UCMH mass and regime of collapse (0-6) for a given cosmic string radius and value of Gmu
!
!Authors: Madeleine Anthonisen
!          maddyanthonisen@hotmail.com
!         Pat Scott
!          pscott@imperial.ac.uk
!Date: 2014/15
!
!All equation numbers refer to v1 of the arXiv posting:
!  Anthonisen, Brandenberger & Scott
!  'Constraints on cosmic strings from ultracompact minihalos', arXiv:1504.something
!
!Inputs:  r_cs    cosmic string radius (kpc)
!         Gmu     cosmic string tension (dimensionless)
!         z_c     latest allowed redshift of UCMH collapse (dimensionless) 
!         zstop   redshift at which UCMHs stop growing (dimensionless)
!         alpha   ratio of loop length to horizon scale at formation (dimensionless)
!         beta    ratio of loop length to loop radius (dimensionless)
!         gamma   cosmic string loop decay constant (dimensionless)
!
!Outputs: mucmh   mass of UCMH (M_solar)
!         ronteq  r_cs / (c*teq) (dimensionless)
!         regime  0-VI: regimes of UCMH collapse (0=No UCMHs)
!         gmu_c1  minimal G\mu for regime I, Scenario A, Eq28
!         gmu_c2  minimal G\mu for regime I, Scenario B, Eq29
!         gmu_c3  G\mu separating regimes III & "No UCMHs", Scenario A, Eq49
!         gmu_c4  G\mu separating regimes V & IV, Scenario B, Eq59(x_f=x_decay)
!         gmu_c5  G\mu separating regimes IV & "No UCMHs", Scenario C, Eq59(x_f=x_c)
!         gmu_c6  G\mu separating regimes VI & VII, Scenario D, Eq63 (x_f=x_decay)
!         gmu_c7  G\mu separating regimes VI & "No UCMHs", Scenario E, Eq63 (x_f=x_c)
!         gmu_c8  G\mu separating regimes VII & "No UCMHs", Scenario D, Eq72
!         gmu_c9  G\mu separating regimes II & III, Scenario A, Eq32

      subroutine dsucmh_mass_cs(r_cs, Gmu, z_c, zstop, alpha, beta, gamma,
     &                          mucmh, ronteq, regime, gmu_c1,gmu_c2,gmu_c3,gmu_c4,gmu_c5,gmu_c6,gmu_c7,gmu_c8,gmu_c9)

      implicit none

      include 'dsmpconst.h'

      double precision, intent(IN) :: r_cs, Gmu, zstop, z_c, alpha, beta, gamma
      double precision, intent(OUT) :: mucmh, ronteq, gmu_c1, gmu_c2, gmu_c3, gmu_c4, gmu_c5, gmu_c6, gmu_c7, gmu_c8, gmu_c9
      integer, intent(OUT) :: regime
      double precision :: Mloop, xstop, x_c, xiRD, xiMD, xdRD, xdMD    
      double precision, parameter :: ct_eq = 1.d3 * c_light * t_eq ! c * teq in m

      xstop = (zeq+1.d0)/(zstop+1.d0)   !x corresponding to zstop
      x_c   = (zeq+1.d0)/(z_c+1.d0)     !x corresponding to zcollapse      

      regime = -1                       !Initial value to be overwritten

      !Loop decay and formation times in RD era and MD era
      xiRD = (beta*r_cs*mperkpc / (alpha*ct_eq))**(1.d0/2.d0) 
      xiMD = (beta*r_cs*mperkpc / (alpha*ct_eq))**(2.d0/3.d0) 
      xdRD = (beta*r_cs*mperkpc / (Gmu*gamma*ct_eq))**(1.d0/2.d0) 
      xdMD = (beta*r_cs*mperkpc / (Gmu*gamma*ct_eq))**(2.d0/3.d0) 

      !Loop mass (Msun)
      Mloop = 1.d12*Gmu*c_light**2*beta*r_cs*mperkpc/(GNewton*SolarMass) 
      
      !Critical values of Gmu delineating different collapse regimes
      gmu_c1 = gamma*(kappa*dmfrac)**(2)*beta**(-6)*alpha**(3)/64 
      gmu_c2 = kappa*dmfrac*beta**(-5.d0/2.d0)*alpha**(3.d0/2.d0)*(r_cs*mperkpc/ct_eq)**(1.d0/2.d0)/8 
      ronteq = r_cs*mperkpc/ct_eq 
      gmu_c9 = 4*Gmu*gamma/(9*beta)-8*beta**2*alpha**(-3.d0/2.d0)*gamma**(1.d0/2.d0)*Gmu**(3.d0/2.d0)/
     &         (kappa*dmfrac)+1024*beta**(5.d0)*alpha**(-3)*(Gmu / (kappa*dmfrac))**(2.d0)
      gmu_c3 = gamma*((2+3*xdRD)*kappa*dmfrac)**(2)*beta**(-6)*alpha**(3)*x_c**(-2)/1600 
      gmu_c4 = kappa*dmfrac*beta**(-3)*alpha**(2)*xiRD/(2*xdMD) 
      gmu_c5 = kappa*dmfrac*beta**(-3)*alpha**(2)*xiRD/(2*x_c) 
      gmu_c6 = 5*kappa*dmfrac*beta**(-3)*alpha**(2)*xiMD**(2)/(8*xdMD-20*xiMD) 
      gmu_c7 = 5*kappa*dmfrac*beta**(-3)*alpha**(2)*xiMD**(2)/(8*x_c-20*xiMD) !Eq 62 xf=x_c 
      gmu_c8 = 5*kappa*dmfrac*beta**(-3)*alpha**(2)*xiMD**(2)*xdMD/(8*x_c*xdMD-8*x_c*xiMD-12*xiMD*xdMD)

      !Determination of UCMH collapse regime.  Each if block corresponds to a different scenario 
      if (xiRD .lt. 1.d0 .and. xdRD .lt. 1.d0)then
        mucmh = 10*xstop*xdRD*Mloop/(2+3*xdRD)                            !Eq48
        if (Gmu .le. gmu_c3) regime = 0
        if (Gmu .gt. gmu_c3 .and. gmu_c9 .le. ronteq) regime = 3
        if (Gmu .le. gmu_c1 .and. gmu_c9 .gt. ronteq) regime = 2
        if (Gmu .gt. gmu_c1) regime = 1
        if (regime .lt. 0) then
          write(*,*) 'Gmu, ronteq, gmu_c1, gmu_c3', Gmu, ronteq, gmu_c1, gmu_c3
          write(*,*) 'No matching regime in if block 1 of dsucmh_mass_cs!'
          stop
        endif
        return
      endif

      if (xiRD .lt. 1.d0 .and. xdMD .gt. 1.d0 .and. xdMd .lt. x_c)then
        mucmh = 2*xstop*Mloop                                             !Eq58
        if (Gmu .le. gmu_c5) regime = 0
        if (Gmu .gt. gmu_c5 .and. Gmu .le. gmu_c4) regime = 5
        if (Gmu .gt. gmu_c4 .and. Gmu .le. gmu_c2) regime = 4
        if (Gmu .gt. gmu_c2) regime = 1
        if (regime .lt. 0) then
          write(*,*) 'Gmu, gmu_c2, gmu_c4, gmu_c5', Gmu, gmu_c2, gmu_c4, gmu_c5
          write(*,*) 'No matching regime in if block 2 of dsucmh_mass_cs!'
          stop
        endif
        return
      endif

      if (xiRD .lt. 1.d0 .and. xdMd .gt. x_c)then
        mucmh = 2*xstop*Mloop                                             !Eq58
        if (Gmu .le. gmu_c5) regime = 0
        if (Gmu .gt. gmu_c5 .and. Gmu .le. gmu_c2) regime = 4
        if (Gmu .gt. gmu_c2) regime = 1
        if (regime .lt. 0) then
          write(*,*) 'Gmu, gmu_c2, gmu_c5', Gmu, gmu_c2, gmu_c5
          write(*,*) 'No matching regime in if block 3 of dsucmh_mass_cs!'
          stop
        endif
        return
      endif

      if (xiMD .gt. 1.d0 .and. xdMd .lt. x_c)then
        ! TB corr 9/15 -> 9./15., 2019-06-01
        mucmh = (2*xstop*(xdMD-xiMD)/(5*xdMD*xiMD)-9./15.)*Mloop            !Eq73
        if (Gmu .le. gmu_c8) regime = 0
        !if (Gmu .gt. gmu_c8 .and. Gmu .lt. gmu_c6) regime = 7
        if (Gmu .gt. gmu_c8) regime = 6
        if (regime .lt. 0) then
          write(*,*) 'Gmu, gmu_c8', Gmu, gmu_c8
          write(*,*) 'No matching regime in if block 4 of dsucmh_mass_cs!'
          stop
        endif
        return
      endif

      if (xiMD .gt. 1.d0 .and. xdMd .gt. x_c)then
        mucmh = (2*xstop/(5*xiMD)-1)*Mloop                                !Eq64
        if (Gmu .le. gmu_c7) regime = 0
        if (Gmu .gt. gmu_c7) regime = 6
        if (regime .lt. 0) then
          write(*,*) 'Gmu, gmu_c7', Gmu, gmu_c7
          write(*,*) 'No matching regime in if block 5 of dsucmh_mass_cs!'
          stop
        endif
        return
      endif

      write(*,*) ' Error in dsucmh_mass_cs.f: did not find valid regime.'
      stop
      
      end subroutine
