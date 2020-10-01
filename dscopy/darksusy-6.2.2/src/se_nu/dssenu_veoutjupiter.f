**********************************************************************
*** Function dssenu_veoutjupiter determines an approximate escape
*** velocity cut to be used instead of 0 to include effects of Jupiter
*** during capture in the Sun. The idea is that while the WIMPs are 
*** being captured by the Sun, Jupiter can deflect the orbits and throw
*** the WIMPs out of the solar system, before they get captured again.
*** How efficient this is depends on the mass of the WIMP and its
*** scattering cross sections. We here use approximate results in
*** A. Peter, PR D79 (2009) 103532, arXiv:0902.1347.
*** More specifically, we interpolate in Fig. 2 to determine how far
*** out the WIMPs can safely reach after the first scatter not to be 
*** disturbed by Jupiter before scattering again. We then translate this
*** to an escape velocity veout of this distance from the Sun, that is
*** the returned and used in the calculations. This is not superprecise,
*** as we do not e.g. distinguish the different regions a21 and a22 in
*** Fig. 2 (i.e. we treat long-lived orbits as captured quickly).
*** Input: mx = WIMP mass in GeV
***        sigsi = spin-independent scattering cross section (cm^2)
***        sigsd = spin-dependent scattering cross section (cm^2)
*** Output: escape velocity of orbits where to cut capture rate
***         calculation (km/s)       
*** Author: Joakim Edsjo, edsjo@fysik.su.se
*** Date: May 18, 2010
**********************************************************************

      real*8 function dssenu_veoutjupiter(mx,sigsi,sigsd)
      implicit none

      real*8 mx,sigsi,sigsd,veout,vcirc,vearth
      real*8 dssenu_jupa1_si,dssenu_jupa2_si,
     &  dssenu_jupa1_sd,dssenu_jupa2_sd
      real*8 fudge
      real*8 asi,asd,amax,ainf
      real*8 a1sig,a2sig,a12sig

c...  Setup
      fudge=100.d0 ! range over which upper a goes from 5.2 (Jupiter) to infinity
      ainf=1000.d0 ! practially infinity
      vearth=29.78d0 ! circ. vel. earth, km/s

c...Determine maximal semi-major axis for SI scattering
      a1sig=dssenu_jupa1_si(mx)
      a2sig=dssenu_jupa2_si(mx)
      a12sig=10**((log10(a1sig)+log10(a2sig))/2.d0) ! mid a2 region
c      if (sigsi.gt.fudge*a1sig) then
c         asi=ainf ! practially infinity
c      elseif (sigsi.gt.a1sig) then
c         asi=(ainf-5.2d0)/(log10(fudge*a1sig)-log10(a1sig))
c     &     *(log10(sigsi)-log10(a1sig)) + 5.2d0
      if (sigsi.gt.a1sig) then
         asi=6.2d0*(sigsi/a1sig)**2 ! dep. on sigsi is just a guess JE
      elseif (sigsi.gt.a12sig) then
         asi=(6.2d0-5.2d0)/(log10(a1sig)-log10(a12sig))
     &     *(log10(sigsi)-log10(a12sig)) + 5.2d0
      elseif (sigsi.gt.a2sig) then
         asi=(5.2d0-4.6d0)/(log10(a12sig)-log10(a2sig))
     &     *(log10(sigsi)-log10(a2sig)) + 4.6d0
      elseif (sigsi.gt.a2sig/10.d0) then
         asi=(4.6d0-3.d0)/(log10(a2sig)-log10(a2sig/10.d0))
     &     *(log10(sigsi)-log10(a2sig/10.d0)) + 3.0d0
      else
         asi=3.0d0
      endif

c...Determine maximal semi-major axis for SD scattering
      a1sig=dssenu_jupa1_sd(mx)
      a2sig=dssenu_jupa2_sd(mx)
      a12sig=10**((log10(a1sig)+log10(a2sig))/2.d0) ! mid a2 region
c      if (sigsd.gt.fudge*a1sig) then
c         asd=ainf ! practially infinity
c      elseif (sigsd.gt.a1sig) then
c         asd=(ainf-5.2d0)/(log10(fudge*a1sig)-log10(a1sig))
c     &     *(log10(sigsd)-log10(a1sig)) + 5.2d0

      if (sigsd.gt.a1sig) then
         asd=6.2d0*(sigsd/a1sig)**2 ! dep. on sigsd is just a guess JE
      elseif (sigsd.gt.a12sig) then
         asd=(6.2d0-5.2d0)/(log10(a1sig)-log10(a12sig))
     &     *(log10(sigsd)-log10(a12sig)) + 5.2d0
      elseif (sigsd.gt.a2sig) then
         asd=(5.2d0-4.8d0)/(log10(a12sig)-log10(a2sig))
     &     *(log10(sigsd)-log10(a2sig)) + 4.8d0
      elseif (sigsd.gt.a2sig/10.d0) then
         asd=(4.8d0-3.0d0)/(log10(a2sig)-log10(a2sig/10.d0))
     &     *(log10(sigsd)-log10(a2sig/10.d0))+3.0d0
      else
         asd=3.0d0
      endif

c...Now pick the largest one, because this is the most effective
c...capture process for this model
      amax=max(asi,asd)

c...Convert maximal semi-major axis amax to an escape velocity  at this
c...radius
      vcirc=vearth/sqrt(amax) ! circular velocity at amax
      veout=sqrt(2.0d0)*vcirc ! escape velocity at amax

c      veout=sqrt(2.0d0)*13.07d0 ! escape velocity at Jupiter

      write(*,*) 'YYY: ',mx,amax,veout
      dssenu_veoutjupiter=veout
      return
      end



**********************************************************************
*** function dssenu_jupa1_si returns the lower boundary of region a1
*** in Fig 2a in A. Peter, PR D79 (2009) 103532, arXiv:0902.1347.
*** SI scattering.
*** Input: WIMP mass in GeV
*** Output: lower bound on cross section for region (cm^2)
*** Author: Joakim Edsjo (edsjo@fysik.su.se)
**********************************************************************

      real*8 function dssenu_jupa1_si(mx)
      implicit none
      real*8 mx
      dssenu_jupa1_si=9.31d-42
      return
      end

**********************************************************************
*** function dssenu_jupa2_si returns the lower boundary of region a2 (a22 and a21)
*** in Fig 2a in A. Peter, PR D79 (2009) 103532, arXiv:0902.1347.
*** SI scattering.
*** Input: WIMP mass in GeV
*** Output: lower bound on cross section for region (cm^2)
*** Author: Joakim Edsjo (edsjo@fysik.su.se)
**********************************************************************

      real*8 function dssenu_jupa2_si(mx)
      implicit none
      real*8 mx,ly
      if (mx.lt.6200.d0) then 
        dssenu_jupa2_si=6.37d-46
      else
        ly=(log10(8.47d-46)-log10(6.37d-46))/(4.d0-log10(6200.d0))
     &   *(log10(mx)-log10(6200.d0))+log10(6.37d-46)
        dssenu_jupa2_si=10**ly
      endif

      return
      end

**********************************************************************
*** function dssenu_jupa1_sd returns the lower boundary of region a1
*** in Fig 2b in A. Peter, PR D79 (2009) 103532, arXiv:0902.1347.
*** SI scattering.
*** Input: WIMP mass in GeV
*** Output: lower bound on cross section for region (cm^2)
*** Author: Joakim Edsjo (edsjo@fysik.su.se)
**********************************************************************

      real*8 function dssenu_jupa1_sd(mx)
      implicit none
      real*8 mx,ly
      if (mx.lt.1550.d0) then 
        dssenu_jupa1_sd=9.30d-40
      else
        ly=(log10(4.53d-39)-log10(9.30d-40))/(4.d0-log10(1550.d0))
     &   *(log10(mx)-log10(1550.d0)) +log10(9.30d-40)
        dssenu_jupa1_sd=10**ly
      endif
      return
      end

**********************************************************************
*** function dssenu_jupa2_sd returns the lower boundary of region a2 (a21 & a22)
*** in Fig 2b in A. Peter, PR D79 (2009) 103532, arXiv:0902.1347.
*** SI scattering.
*** Input: WIMP mass in GeV
*** Output: lower bound on cross section for region (cm^2)
*** Author: Joakim Edsjo (edsjo@fysik.su.se)
**********************************************************************

      real*8 function dssenu_jupa2_sd(mx)
      implicit none
      real*8 mx,ly
      if (mx.lt.803.d0) then 
        dssenu_jupa2_sd=6.47d-44
      else
        ly=(log10(5.93d-43)-log10(6.47d-44))/(4.d0-log10(803.d0))
     &   *(log10(mx)-log10(803.d0)) +log10(6.47d-44)
        dssenu_jupa2_sd=10**ly
      endif
      return
      end




