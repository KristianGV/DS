      subroutine dsddffsi(q,a,z,ff,ierr)
c_______________________________________________________________________
c  Spin-independent form factor for direct detection.
c  input:
c    q : real*8  : momentum transfer in GeV ( q=sqrt(2*M*E) )
c    a : integer : mass number
c    z : integer : atomic number
c  output:
c    ff : |F(q)|^2, the square of the form factor
c  author: paolo gondolo (paolo@physics.utah.edu) 2004-2008
c=======================================================================
      implicit none
      include 'dsddcom.h'
      include 'dsio.h'
      include 'dsmpconst.h'

      real*8 q
      integer a,z,ierr
      real*8 ff
      
      ierr=0
      
      if (a.eq.1) then
         ff = 1.d0

c best available form factor

      else if (ddsf(1)(2:).eq.'best') then
         call dsddfffb(q,a,z,ff,ierr)
         if (ierr.ne.0) then
            if (prtlevel.gt.1) write (*,*) 
     &           'dsddffsi: switching to SOG....'
            call dsddffsog(q,a,z,ff,ierr)
            if (ierr.ne.0) then
               if (prtlevel.gt.1) write (*,*) 
     &              'dsddffsi: switching to Fermi....'
               call dsddfffermi(q,a,z,ff,ierr)
               if (ierr.ne.0) then
                  if (prtlevel.gt.1) write (*,*) 
     &                 'dsddffsi: switching to Lewin-Smith....'
                  call dsddffls(q,a,z,ff,ierr)
                  ierr=0
               endif
            endif
         endif
         
c form factor at q=0

      else if (ddsf(1)(1:1).eq.'0') then
         ff=1.d0

c Helm form factor with Lewin-Smith parameters

      else if (ddsf(1)(2:).eq.'L-S') then
         call dsddffls(q,a,z,ff,ierr)

c Helm form factor as in DarkSUSY 4.1 (with corrected expansion of f, PG 20080216)

      else if (ddsf(1)(2:).eq.'ds4.1') then
         call dsddffh41(q,a,z,ff,ierr)

c old obsolete (incorrect) form factor from Gould
      else if (ddsf(1)(2:).eq.'gould') then
        call dsddffgould(q,a,z,ff,ierr)

c gaussian form factor

      else if (ddsf(1)(2:).eq.'gauss') then
        call dsddffgauss(q,a,z,ff,ierr)

cc Sum of Gaussian form factor                                cc

      else if (ddsf(1)(2:).eq.'SOG') then
         call dsddffsog(q,a,z,ff,ierr)

cc Fourier-Bessel form factor                                 cc

      else if (ddsf(1)(2:).eq.'FB') then
         call dsddfffb(q,a,z,ff,ierr)

cc Fermi Integration form factor                              cc

      else if (ddsf(1)(2:).eq.'Fermi') then
         call dsddfffermi(q,a,z,ff,ierr)

c unrecognized form factor label

      else
         write (*,*) 
     &        'dsddffsi : Unrecognized spin-independent form factor',ddsf(1)
         ierr=1
         return
      endif

      if (abs(ff).gt.1.001d0) then
        write (*,*) 
     &       'dsddffsi : SI form factor too large. Suspect numerical extrapolation problem.'
        ierr=2
        return
      endif
      
      return
      end
