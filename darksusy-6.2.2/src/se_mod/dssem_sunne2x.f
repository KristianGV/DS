***********************************************************************
*** dssem_sunne2x takes an input number density of electrons and
*** converts this to a fractional solar radius, x.
*** Input: n_e [cm^-3]
*** Output: x = r/r_sun [0,1]
*** See dssem_sunread for information about which solar model is used.
*** Author: Joakim Edsjo, edsjo@fysik.su.se
*** Date: 2006-03-27
***********************************************************************

      real*8 function dssem_sunne2x(ne)
      implicit none

      include 'dssem_sun.h'
      include 'dsmpconst.h'

      real*8 ne,nepl,lne
      integer j

      j=0                       ! initialize
      
c...Check if data file is loaded
      call dssem_sunread

      lne=log10(ne/n_avogadro) ! this is what is tabulated
      if (lne.ge.sdne(1)) then
        dssem_sunne2x=0.0d0
        return
      endif

      if (lne.le.sdne(sdnne)) then
        dssem_sunne2x=sdrne(sdnne)
        return
      endif

      call dshunt(sdne,sdnne,lne,j)
      if (j.lt.sdnne) goto 20
      
      dssem_sunne2x=0.0d0
      return

 20   nepl=(lne-sdne(j))/(sdne(j+1)-sdne(j))

      dssem_sunne2x=sdrne(j)*(1.0d0-nepl)+sdrne(j+1)*nepl

      return

      end
