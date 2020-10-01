*****************************************************************
*** suboutine dsanyield_dbset sets options for antideuteron
*** yields
*** Input:
***    yieldk = 59: old spherical coalescence model (should only
***                 be used for comparison, obsolete
***             61: Pythia 6 MC results (default)      
***             62: Pythia 8 MC results (not yet implemented)
***             63: Herwig++ MC results (not yet implemented)
***              0: default, i.e. the same as 61
***           -100: do not change yieldk      
***    Note: integrated or differential yields are chosen
***    explicitly when calling dsanyield_sim.
***      
***    p0bar: coalescence momentum chosen. This is the *deviation*
***           from the best fit value from Aleph data (for MC yields)
***           Hence, the used coalescence momentum is given by
***               p0 = p0_best-fit + p0bar*sigma_p0      
***           p0bar is hence the dimenionless number of sigma away (+ or -)
***           from the best fit value. The best fit value is determined
***           from comparisons with Aleph measurements and is determined
***           for each MC separately. sigma_p0 is defined as
***           p0_high - p0_best-fit for p0>p0_best-fit and
***           p0_best-fit-p0_low for p0<p0_best-fit
***           Default: 0.0 (i.e. best fit p0 will be used)
***           If -100.d0, p0bar is not changed in this routine.      
***
*** Author: Joakim Edsjo, edsjo@fysik.su.se
*** Date: April, 2016      
*****************************************************************

      subroutine dsanyield_dbset(yieldk,p0bar)
      implicit none

      include 'dsanyieldcom.h'
      integer yieldk
      real*8 p0bar

      if (yieldk.ge.59.and.yieldk.le.63) then
         dbflxk=yieldk
      elseif (yieldk.eq.0) then
         dbflxk=61
      endif

      if (p0bar.gt.-99.99d0) then
         dbp0bar=p0bar
      endif
      
      return

      end
