      subroutine dsdddrde(t,e,n,a,z,stoich,rsi,rsd,modulation)
c_______________________________________________________________________
c
c  type : commonly used
c  desc : Differential WIMP-nucleus recoil rates
c
c  differential recoil rate
cc
c  input:
c    e : real*8              : nuclear recoil energy in keV
c    n : integer             : number of nuclear species
c    a : n-dim integer array : mass numbers
c    z : n-dim integer array : atomic numbers
c    stoich : n-dim integer array : stoichiometric coefficients
c    t : real*8 : time in days from 12:00UT Dec 31, 1999
c    modulation : integer : 0=no modulation 1=annual modulation
c      with no modulation (ie without earth velocity), t is irrelevant
c  output:
c    rsi : real*8 : spin-independent differential rate
c    rsd : real*8 : spin-dependent differential rate
c  units: counts/kg-day-keV
c
c  NB: This assumes that neutralinos make up 100% of the local DM density
c      If not, rsi and rsd need to be multiplied (by hand) by the fraction
c      of local DM that consists of neutralinos!
c
c  author: paolo gondolo (paolo@physics.utah.edu) 2004
c
c  2018-01-19: removed rescaling factor for local neutralino density (rhox)
c              from common block [TB]
c=======================================================================
      implicit none
      include 'dshmcom.h'
      include 'dsnuclides.h'
      include 'dsmpconst.h'

      integer n,a(n),z(n),stoich(n)
      real*8 e,t,rsi,rsd,siff,sdff,cw(n),cwtmp
      real*8 sigij(27,27)
      real*8 mx,dsmwimp
      integer i,ii,dsnucldindx,ierr
      real*8 mni,muxi,norm,vmin,eta,v
      ! norm = (kg day keV)/(cm GeV^2 km/s) 
      !      = 86400/10 * [kg/(GeV/c^2)] * [c/(km/s)]^2
      parameter (norm=4.355983308d41)
      integer modulation,modulatio
      common /ddmodul/ modulatio

      cwtmp = 0.d0
      do i=1,n
        cw(i) = a(i)*stoich(i)
        cwtmp = cwtmp+cw(i)
      enddo
      do i=1,n
        cw(i) = cw(i)/cwtmp
      enddo

      mx=dsmwimp()
      
      rsi = 0.d0
      rsd = 0.d0
      do i=1,n
        ii=dsnucldindx(a(i),z(i))
        if (ii.ne.0) then
          mni=nucldm(ii)
          muxi = mni*mx/(mni+mx)
          vmin = sqrt(mni*e*1.d-6/2.d0/muxi**2)*c_light
          call dsddeta(vmin,t,eta)
          ! We now need to calculate the cross section including form factor.
          ! The form factor depends on the recoil energy, not the velocity
          ! but we need to set the velocity such that Er can be reached
          v=vmin*1.1d0  ! to make sure we can reach Er, siff and sdff below
                        ! would otherwise be zero.
          call dsddsigma(v,e,a(i),z(i),sigij,ierr)
          siff=sigij(1,1)
          sdff=sigij(4,4)
          if (ierr.ne.0) then
            rsi=0.d0
            rsd=0.d0
            write (*,*) 'WARNING in dsdddrde: dsddsigma returned ierr=',ierr
            return
          endif
          rsi = rsi + cw(i)*siff/2.d0/mx/muxi**2*eta
          rsd = rsd + cw(i)*sdff/2.d0/mx/muxi**2*eta
        endif
      enddo
c... see comment in header: this assumes that neutralino density = local DM density !       
      rsi = norm*rsi*rho0  
      rsd = norm*rsd*rho0

      return
      end
