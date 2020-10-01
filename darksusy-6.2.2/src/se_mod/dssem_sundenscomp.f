***********************************************************************
*** dssem_sundenscomp gives the number density of nucleons of atomic
*** number Z and isotope number i (as given in the definitions/files
*** in dssem_sunread) per cm^3
*** input: radius - in meters
***        z: atomic number of element
***        iso: isotope number (usually 1 for the most common and higher
***             for more rare isotopes)
*** Author: joakim edsjo
*** Date: 2003-11-26
*** Modified: 2006-03-21 (atomic mass unit fix (was off by 6%)) JE
***           2009-11-10 new more general routine with 'all' elements
***             and different isotopes
***********************************************************************

      real*8 function dssem_sundenscomp(r,z,iso)
      implicit none
      include 'dsmpconst.h'
      real*8 r,dssem_sundens,fraction,dssem_sunmfrac,m_a
      integer z,iso

      include 'dssem_sun.h'

      if (z.lt.1.or.z.gt.zmax.or.iso.lt.0.or.iso.gt.isomax) then
        write(*,*) 'WARNING in dssem_sundenscomp: illegal element type: ',
     &    'z=',z,'  isotope=',iso
      endif

c...m_a is the mass in atomic mass units, amu (where C12 has mass 12 amu)
c...sdma(4) is the mass of C12
      m_a = sdma(z,iso)/sdma(6,1)*12.0d0  ! JE corr 060321
      fraction=dssem_sunmfrac(r,z,iso)

c...m_a is now in atomic units u which is what we want

      dssem_sundenscomp=fraction*dssem_sundens(r)/m_a*n_avogadro ! number/cm^3

      return
      end
