      real*8 function dsrdbw_get(p,istat)
**********************************************************************
*** This routine returns the Breit-wigner fit to the invariant
*** annihilation rate. It only returns a value if
***   a) the resonance has been found and fit within the required
***      accuracy
***   b) the p value is within the range where the fit is good      
***
*** The setup of the resonances is done with the routine dsrdbw_setup.
*** that has to be called before this routine is called.
***
*** Upon return, istat=0: within p range of well-fit resonance
***              istat=1: outside or p range, returns 0      
***
*** Author: Joakim Edsjo, edsjo@fysik.su.se
*** Date: 2018-10-26      
**********************************************************************      

      implicit none
      include 'dsrdcom.h'

      integer istat,i
      real*8 p
      real*8 dsrdbreit_wigner
      real*8 pa,pb
      real*8 ea,eb

      istat=1
      dsrdbw_get=0.d0
      
      do i=1,nres
         if (resinc(i).and.resfit(i)) then
            ea=rgev(i)-resgamma(i)*rwid(i)
            eb=rgev(i)+resgamma(i)*rwid(i)
            pa=sqrt(ea**2/4-mco(1)**2)
            pb=sqrt(eb**2/4-mco(1)**2)
            if (p.ge.pa.and.p.le.pb) then ! within resonance
               dsrdbw_get=resnorm(i)*dsrdbreit_wigner(rgev(i),rwid(i),
     &              resalpha(i),mco(1),p)+resconst(i)
               istat=0
            endif
         endif
      enddo

      return
      end
      
         
         
            
         


         
         
