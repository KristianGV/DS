***********************************************************************
*** dskdcfull returns the integrand for the integration over the SM 
*** (or other scattering partner) momenta that appears in the collision 
*** term. Called by dskdgammarate. itegrand is multiplied by 1d20 to work 
*** better with the integration routine
***
*** author: torsten bringmann (troms@physto.se), 2010-01-23
*** updates: 2014-05-09 (removed model-dependence)
***          2017-04-28 (added two options for quark scattering)
***          2018-05-14 (added option for non-SM scattering)
***********************************************************************

      real*8 function dskdcint(kint)
      implicit none

      include 'dsmpconst.h'
      include 'dskdcom.h'

      real*8 kint,omega,tmpres, sign
      real*8 scatterdist,dskdm2, amp
      integer Stype

      Stype = 1
      tmpres = 0.0d0


 10   omega=kint
      amp = dskdm2(omega,Stype) ! NB: kint->omega (energy) on return!

      if (amp.gt.0) then 

        sign = 1.d0  ! +1.d0 (-1.d0) for fermionic (bosonic) scattering partners
        if (Stype.gt.12) then
          sign = -1.d0
          if (Stype.gt.100) then
            if (BSMscattfermion(Stype-100)) sign = 1.d0
          endif
        endif        
        scatterdist = 1.d0/(exp(omega/Tint)+sign)  ! assume thermal distribution of scattering partners
        
c... this would be the expression for the simplified amplitude
c      call dskdm2simp(m0,Stype,cn,n)
c      tmpres=tmpres+cn*(omega/m0)**2*scatterdist*(1.d0-sign*scatterdist)/omega
        tmpres = tmpres + amp*scatterdist*(1.d0-sign*scatterdist)/omega
      endif
  
      Stype = Stype+1

      if (Stype.le.6) goto 10

c... add all quarks above tqcdmin 
      if (quark_how.eq.1.and.Stype.le.12.and.Tint.gt.tqcdmin) goto 10

c... OR add only light quarks above tqcdmax 
      if (quark_how.eq.2.and.Stype.le.9.and.Tint.gt.tqcdmax) goto 10

c... now add SM bosons as scattering partners
      if (Stype.le.13) Stype = 13
      if (Stype.le.13) goto 10 ! presently only photons supported

c... now add contribution from non-SM scattering partners
      if (nBSMscatt.gt.0) then
        if (Stype.lt.100) Stype=101
        if (Stype.le.(100+nBSMscatt)) goto 10
      endif

c... this is the *remaining* part of Eq. A7 in 1603.04884, taking into account 
c... the erratum to hep-ph/0612238
c... NB: The 'missing' factor of 1/kdof(dm) is already taken into account in 
c... the definition of |M|^2!
      dskdcint = tmpres*kint**5/(48.*pi**3*m0**3*Tint)
      dskdcint=dskdcint*1.d25 ! rescale for integration to work fine
                      

      return
      end

