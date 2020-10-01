*****************************************************************************
***   subroutine dsaninit initializes and loads (from disk) the common
***   block variables needed by the other halo yield routines.
***   author: joakim edsjo
***   edsjo@physto.se date: 96-10-23 (based on dsmuinit.f version
***   3.21) 
***   modified: 98-01-26
***   modified: 09-10-20 pat scott pat@fysik.su.se
***   modified: May, 2014, Joakim Edsjo
***   modified: April, 2016, Joakim Edsjo, to read new simulations tables
***     (including dbars)      
*****************************************************************************

      subroutine dsanyield_init
      implicit none
      include 'dsio.h'
      include 'dsanyieldcom.h'

c------------------------ variables ------------------------------------

      integer i
 
      logical first
      data first/.true./
      save first

c...startup values
      data dsanyieldinitcalled/.false./

c...dbar setup
c...dbp0bar is deviations from best-fit p0 in units of sigma, i.e.
c...p0 = p0_best_fit + p0_sigma*dbp0bar      
      data dbp0bar/0.0d0/   !Default is to use best fit p0 value for dbar calcs
c...Default yield code for dbars, otpions are
c...     yieldk = 59 - old spherical coalescence model
c...     yieldk = 61 Pythia 6 MC
c...     yieldk = 62 Pythia 8 MC (not implemented yet)
c...     yieldk = 63 Herwig++ MC (not implemented yet)
      data dbflxk/61/ ! Default dbar MC, Pythia 6
      

c----------------------------------------- set-up common block variables

      dsanyieldinitcalled=.true.

      if (first) then

c...Setups
         anftype='a' ! ascii input files
         ansmooth=0

c...masses of annihilation products used in the simulation
c...masses of annihilation products d, u, s, c, b, t, glue, W+W-, Z0Z0, 
c...mu+mu- tau+ tau-
        msim(1)=1.0       ! d
        msim(2)=1.0       ! u
        msim(3)=1.0       ! s
        msim(4)=1.35      ! c
        msim(5)=4.8       ! b
        msim(6)=173.5     ! top
        msim(7)=0.0       ! glue
        msim(8)=80.25     ! w
        msim(9)=91.2      ! z
        msim(10)=0.10566  ! mu
        msim(11)=1.7841   ! tau

c...lower and upper bounds on WIMP masses for which different
c...annihilation channels can be used
        lb(1)=3.0
        lb(2)=3.0
        lb(3)=3.0
        lb(4)=3.0
        lb(5)=6.0
        lb(6)=174.0
        lb(7)=3.0
        lb(8)=80.2
        lb(9)=91.2
        lb(10)=3.0
        lb(11)=3.0

        do i=1,11
          ub(i)=20000.0
        enddo

c...masses for simulation corresponding to mass index i
        mi(1)=3.0
        mi(2)=4.0
        mi(3)=5.0
        mi(4)=6.0
        mi(5)=8.0
        mi(6)=10.0
        mi(7)=15.0
        mi(8)=20.0
        mi(9)=25.0
        mi(10)=35.0
        mi(11)=50.0
        mi(12)=80.2
        mi(13)=91.2
        mi(14)=100.0
        mi(15)=150.0
        mi(16)=174.0
        mi(17)=200.0
        mi(18)=250.0
        mi(19)=350.0
        mi(20)=500.0
        mi(21)=750.0
        mi(22)=1000.0
        mi(23)=1500.0
        mi(24)=2000.0
        mi(25)=3000.0
        mi(26)=5000.0
        mi(27)=7500.0
        mi(28)=10000.0
        mi(29)=15000.0
        mi(30)=20000.0

c...lowest mass index for channel j
        milow(1)=1   ! d d-bar
        milow(2)=1   ! u u-bar
        milow(3)=1   ! s s-bar
        milow(4)=1   ! c c-bar
        milow(5)=4   ! b b-bar
        milow(6)=16  ! t t-bar
        milow(7)=1   ! gluons
        milow(8)=12  ! w+ w-
        milow(9)=13  ! z z
        milow(10)=1  ! mu+ mu-
        milow(11)=1  ! tau+ tau-

c...initialize eindex array where the energies for the bins are stored
c...integrated yields (lower end of each bin)
        do i=-1,zn
          zindex(i,1)=dble(i)/dble(zn)
        enddo

c...differential yields (center of each bin)
        do i=-1,zn
          zindex(i,2)=dble(i)/dble(zn)+0.5d0/dble(zn)
          dz(i)=1.0d0/dble(zn)
        enddo

c...initialize dbzindex array where the energies for the bins are stored
c...integrated yields (lower end of each bin), for dbar which has diff. binning
        do i=-1,zndb
          dbzindex(i,1)=dble(i)/dble(zndb)
        enddo

c...differential yields (center of each bin)
        do i=-1,zndb
          dbzindex(i,2)=dble(i)/dble(zndb)+0.5d0/dble(zndb)
          dbdz(i)=1.0d0/dble(zndb)
        enddo


c...initialize dbdpindex array where the upper end of each bin is stored
c ...(the tables are accumulative)
        do i=-1,dpn
          dbdpindex(i)=dble(i+1)/dble(dpn)*300 ! MeV
        enddo

        first=.false.
      endif

      return

      end

