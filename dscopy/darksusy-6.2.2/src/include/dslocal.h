      !Common variables and constants for dsplocalmass
      !Pat Scott Jan 2011 patscott@physics.mcgill.ca

      real*8 ly,pc,Msun,rhobg,M_200,c,r_virial
      real*8 r_s3,r_sun_gal,del,theeta,dmin_internal
            
      parameter(ly=9.461d15)            ! 1 light year in meters
      parameter(pc=30.857d15)           ! 1 pc in meters
      parameter(Msun=1.9891d30)         ! solar mass in kg
      parameter(c=18.d0)                ! Conc. parameter of halo
      parameter(M_200=9.4d11*Msun)      ! Halo's mass in kg
      parameter(rhobg=2.67d-27)         ! bkgr density in kg/m3
      parameter(r_sun_gal=26d3*ly)      ! distance btw Sun and galactic centre in meters
      parameter(r_s3=(M_200/(16.0d0*datan(1.d0)*rhobg*200d0/3d0*c**3))) !scale radius^3 i meter
      parameter(del=200d0/3d0*c**3/(log(c+1d0)-c/(c+1d0)))
      parameter(r_virial=c*r_s3**(1.d0/3.d0))   ! virial radius of the Milky Way in meters


      common /localm/dmin_internal
      save /localm/
      common /fourkpc/theeta
      save /fourkpc/

