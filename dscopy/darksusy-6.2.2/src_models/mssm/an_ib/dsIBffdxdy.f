*****************************************************************************
*** The function dsIBffdxdy gives the full analytical expressions for
*** the differential IB photon yield (dxdy) from fermion final states, 
*** normalized to the annihilation rate into fermion pairs f fbar
***
*** The kinematic variables x,y are
***
*** x = E_gamma/mx
*** y = (p+k)^2/(4 mx^2),
***
*** where p denotes the fermion momentum and k the photon momentum.
*** (note that the expressions above and below only apply to the v->0 limit)
***
*** author: Torsten Bringmann (bringman@sissa.it)
*** date: 2007-07-05
*** update: 2016-04-04 (switched to analytic tree-level expressions to avoid
***                     recursive call to dssigmav0)
*****************************************************************************

      real*8 function dsIBffdxdy(IBch,x,y)

      implicit none
      include 'dsmssm.h'
      include 'dsmpconst.h'
      include 'dshmcom.h'
      include 'dsidtag.h'
      include 'dsio.h'

      real*8 x,y,tmpresult,tmpdecay
      real*8 msq2bodyds, msq2body
      integer IBch, kf, ksf1, ksf2                ! IB channel and particle codes
      real*8 dsIBffdxdy_1,dsIBffdxdy_2,dsIBffdxdy_3
      real*8 dsIBffdxdy_4,dsIBffdxdy_5,dsIBffdxdy_6
      real*8 dsIBffdxdy_7,dsIBffdxdy_8,dsIBfsrdxdy

      character*12 memory    ! to suppress multiple error messages
      save memory
      data memory /'____________'/


c------------- couplings and masses relevant for IB -------------------

      real*8 C11, C12, C151, C152, C21, C22, C251, C252, CZ5, CH
      real*8 m0,msf2, msf1, mz, mf, mh03, Gh03, GZ, zf, hf

        
      dsIBffdxdy=0d0
      tmpresult=0d0
      tmpdecay=0d0

c...determine fermion and sfermion particle codes
     
      if (IBch.eq.4) then
        kf   = ke
        ksf1 = ksl_flav(1,1)
        ksf2 = ksl_flav(1,2)
      elseif (IBch.eq.5) then
        kf   = kmu
        ksf1 = ksl_flav(2,1)
        ksf2 = ksl_flav(2,2)
      elseif (IBch.eq.6) then
        kf   = ktau
        ksf1 = ksl_flav(3,1)
        ksf2 = ksl_flav(3,2)
      elseif (IBch.eq.7) then
        kf   = ku
        ksf1 = ksqu_flav(1,1)
        ksf2 = ksqu_flav(1,2)
      elseif (IBch.eq.8) then
        kf   = kd
        ksf1 = ksqd_flav(1,1)
        ksf2 = ksqd_flav(1,2)
      elseif (IBch.eq.9) then
        kf   = kc
        ksf1 = ksqu_flav(2,1)
        ksf2 = ksqu_flav(2,2)
      elseif (IBch.eq.10) then
        kf   = ks
        ksf1 = ksqd_flav(2,1)
        ksf2 = ksqd_flav(2,2)
      elseif (IBch.eq.11) then
        kf   = kt
        ksf1 = ksqu_flav(3,1)
        ksf2 = ksqu_flav(3,2)
      elseif (IBch.eq.12) then
        kf   = kb
        ksf1 = ksqd_flav(3,1)
        ksf2 = ksqd_flav(3,2)
      else
    
        write(*,*) 'ERROR in dsIBffdxdy: called with IBch = ', IBch 
        return  

      endif

c...set up masses and widths

      m0   = mass(kn(1))
      mz   = mass(kz)
      mh03 = mass(kh3)
      Gh03 = width(kh3)
      GZ   = width(kz)

      mf   = mass(kf)
      msf1 = mass(ksf1)
      msf2 = mass(ksf2)
      Gh03 = width(kh3)
  
c...set up couplings
 
      C11   = abs(gl(ksf1,kf,kn(1)))**2 - 
     -        abs(gr(ksf1,kf,kn(1)))**2
      C12   = abs(gl(ksf2,kf,kn(1)))**2 - 
     -        abs(gr(ksf2,kf,kn(1)))**2
      C151  = abs(gl(ksf1,kf,kn(1)))**2 + 
     -        abs(gr(ksf1,kf,kn(1)))**2
      C152  = abs(gl(ksf2,kf,kn(1)))**2 + 
     -        abs(gr(ksf2,kf,kn(1)))**2
      C21   = imag(gr(ksf1,kf,kn(1))*conjg(gl(ksf1,kf,kn(1))))
      C22   = imag(gr(ksf2,kf,kn(1))*conjg(gl(ksf2,kf,kn(1))))
      C251  = dble(gr(ksf1,kf,kn(1))*conjg(gl(ksf1,kf,kn(1))))
      C252  = dble(gr(ksf2,kf,kn(1))*conjg(gl(ksf2,kf,kn(1))))
      CZ5   = dble(gl(kz,kn(1),kn(1)))*
     -        dble(gr(kz,kf,kf)-gl(kz,kf,kf))
      CH    = imag(gr(kh3,kn(1),kn(1)))*imag(gr(kh3,kf,kf))

c...import IB expressions for |M|**2 from form/mathematica; 
c...for light leptons and quarks, take the simplified expression with mf=0:

      if ((IBch.eq.4).or.(IBch.eq.5).or.(IBch.eq.7).or.
     -    (IBch.eq.8).or.(IBch.eq.10)) then
            tmpresult=
     -        -2*m0**6*(-1 + x)*(x**2 - 2*x*y + 2*y**2)*
     -       ((C11/((msf1**2 + m0**2*(1 - 2*y))*
     -             (msf1**2 + m0**2*(1 - 2*x + 2*y))) + 
     -         C12/
     -          ((msf2**2 + m0**2*(1 - 2*y))*
     -          (msf2**2 + m0**2*(1 - 2*x + 2*y))))**2 + 
     -        (C151/
     -          ((msf1**2 + m0**2*(1 - 2*y))*
     -            (msf1**2 + m0**2*(1 - 2*x + 2*y))) + 
     -          C152/
     -          ((msf2**2 + m0**2*(1 - 2*y))*
     -            (msf2**2 + m0**2*(1 - 2*x + 2*y))))**2)
           endif

c...for heavy leptons and quarks, add
c...contributions from all linear independent combinations of coupling
c...constants separately, utilizing various symmetries

      if ((IBch.eq.6).or.(IBch.eq.9).or.(IBch.eq.11)
     -    .or.(IBch.eq.12)) then

          zf=(mz**2*(16*m0**4+GZ**2*mz**2-8*m0**2*mz**2
     -       +mz**4))/(-4*m0**2 + mz**2)**2
          hf=((-16*m0**4-Gh03**2*mh03**2+8*m0**2*mh03**2
     -       -mh03**4)*mf)/(2.*m0*(4*m0**2 - mh03**2))

      tmpresult=
     -     C11*C12*dsIBffdxdy_1(x,y,m0,mf,msf1,msf2)
     -  +  C11**2*dsIBffdxdy_1(x,y,m0,mf,msf1,msf1)/2.
     -  +  C12**2*dsIBffdxdy_1(x,y,m0,mf,msf2,msf2)/2.
     -
     -  +  C151*C152*dsIBffdxdy_2(x,y,m0,mf,msf1,msf2)
     -  +  C151**2*dsIBffdxdy_2(x,y,m0,mf,msf1,msf1)/2.
     -  +  C152**2*dsIBffdxdy_2(x,y,m0,mf,msf2,msf2)/2
     -
     -  +  C21*C22*dsIBffdxdy_3(x,y,m0,mf,msf1,msf2)
     -  +  C21**2*dsIBffdxdy_3(x,y,m0,mf,msf1,msf1)/2.
     -  +  C22**2*dsIBffdxdy_3(x,y,m0,mf,msf2,msf2)/2.
     -
     -  +  C251*C252*dsIBffdxdy_4(x,y,m0,mf,msf1,msf2)
     -  +  C251**2*dsIBffdxdy_4(x,y,m0,mf,msf1,msf1)/2.
     -  +  C252**2*dsIBffdxdy_4(x,y,m0,mf,msf2,msf2)/2.
     -
     -  +  C151*C252*dsIBffdxdy_5(x,y,m0,mf,msf1,msf2)
     -  +  C151*C251*dsIBffdxdy_5(x,y,m0,mf,msf1,msf1)
     -  +  C152*C252*dsIBffdxdy_5(x,y,m0,mf,msf2,msf2)
     -  +  C152*C251*dsIBffdxdy_5(x,y,m0,mf,msf2,msf1)
     -
     -  +  C151*CZ5*dsIBffdxdy_6(x,y,m0,mf,msf1)/zf
     -  +  C152*CZ5*dsIBffdxdy_6(x,y,m0,mf,msf2)/zf
     -  +  C151*CH*dsIBffdxdy_6(x,y,m0,mf,msf1)/hf
     -  +  C152*CH*dsIBffdxdy_6(x,y,m0,mf,msf2)/hf
     -
     -  +  C251*CZ5*dsIBffdxdy_7(x,y,m0,mf,msf1)/zf
     -  +  C252*CZ5*dsIBffdxdy_7(x,y,m0,mf,msf2)/zf
     -  +  C251*CH*dsIBffdxdy_7(x,y,m0,mf,msf1)/hf
     -  +  C252*CH*dsIBffdxdy_7(x,y,m0,mf,msf2)/hf
     -
     -  +  dsIBffdxdy_8(x,y,CZ5,CH,m0,mf,mz,mh03,GZ,Gh03)

      endif   ! contribution from heavy fermions


c... this is the analytically obtained 2-body 
c... annihilation amplitude squared (in the same v->0 limit, 
c... and up to color factors, as considered for the 3-body case) :
c 
      msq2body=  
     -  4*(m0**4 - m0**2*mf**2)*
     -   (C21/(m0**2 - mf**2 + msf1**2) + 
     -      C22/(m0**2 - mf**2 + msf2**2))**2 + 
     -  m0**2*((2*C251*m0 + C151*mf)/(m0**2 - mf**2 + msf1**2) + 
     -      (2*C252*m0 + C152*mf)/(m0**2 - mf**2 + msf2**2))**2 + 
     -  (4*CZ5*m0**2*mf*((2*C251*m0 + C151*mf)/
     -        (m0**2 - mf**2 + msf1**2) + 
     -       (2*C252*m0 + C152*mf)/(m0**2 - mf**2 + msf2**2) + 
     -       (CZ5*mf)/mz**2)*(-4*m0**2 + mz**2)**2)/
     -   (mz**2*(GZ**2*mz**2 + (-4*m0**2 + mz**2)**2)) + 
     -  (8*CH*m0**3*(2*CH*m0 + 
     -       (-4*m0**2 + mh03**2)*
     -        ((2*C251*m0 + C151*mf)/(m0**2 - mf**2 + msf1**2) + 
     -          (2*C252*m0 + C152*mf)/(m0**2 - mf**2 + msf2**2)) + 
     -       (2*CZ5*mf*(-(Gh03*GZ*mh03*mz*(4*m0**2 - mz**2)) + 
     -            (-4*m0**2 + mh03**2)*(-4*m0**2 + mz**2)**2))/
     -        (mz**2*(GZ**2*mz**2 + (-4*m0**2 + mz**2)**2))))/
     -   (Gh03**2*mh03**2 + (-4*m0**2 + mh03**2)**2)

        msq2bodyds = msq2body


c... In a previous code version, the tree-level rates were obtained from dssigmav0
c... This is no longer possible because now dssigmav0 calls the IB routines 
c... (via dsandwdcosnn). The consistency check below still works if the conflicting
c... lines in dsandwdcosnn are commented out as well.

c      if (IBch.eq.4) msq2bodyds=dssigmav0(11,-11)  ! e- e+
c      if (IBch.eq.5) msq2bodyds=dssigmav0(13,-13)  ! mu- mu+
c      if (IBch.eq.6) msq2bodyds=dssigmav0(15,-15)  ! tau- tau+
c... remove color factor for quarks from tree--level result
c... (not yet taken into account in 3body results -> ratio OK) 
c      if (IBch.eq.7) msq2bodyds=dssigmav0(2,-2)/3. ! dd; color factor for quarks
c      if (IBch.eq.8) msq2bodyds=dssigmav0(1,-1)/3. ! u ubar
c      if (IBch.eq.9) msq2bodyds=dssigmav0(4,-4)/3. ! c cbar
c      if (IBch.eq.10) msq2bodyds=dssigmav0(3,-3)/3.! s sbar 
c      if (IBch.eq.11) msq2bodyds=dssigmav0(6,-6)/3.! t tbar
c      if (IBch.eq.12) msq2bodyds=dssigmav0(5,-5)/3.! b bbar

c      msq2bodyds=msq2bodyds*32*pi*m0**2/sqrt(1-mf**2/m0**2)/gev2cm3s

c... compare msq2body to msq2bodyds:
c      if (((1d0-msq2body/msq2bodyds)**2.ge.0.01)
c     -     .and.(memory.ne.idtag)) then
c        memory=idtag
c        write(*,*) '*****'
c        write(*,*) 'model ',idtag,' - ', 'warning from ',
c     -             'dsIBffdxdy:'
c        write(*,*) 'Analytic v=0 result for tree-level cross section in '.
c     -             'channel ',IBch, differs from the result from dssigmav0',
c     -             ' by a factor of',msq2bodyds/msq2body
c        write(*,*) '*****'
c       endif

c...check that result is positive
      if (0.gt.tmpresult) then
        if (m0**2*tmpresult.lt.(-1D-12)
     -     .and.(idtag.ne.memory).and.prtlevel.ge.1) then
          write(*,*) '*****'
          write (*,*) 'Error in dsIBffdxdy (channel:',IBch,
     -                '): negative |M|^2 for model ',idtag,'.'
          write (*,*) 'Setting corresponding contributions to zero...'
          write(*,*) '*****'
        endif
        return
      endif

c...take into account different electric charges for quarks

      if ((IBch.eq.7).or.(IBch.eq.9).or.(IBch.eq.11))
     -    tmpresult=tmpresult*4/9.
      if ((IBch.eq.8).or.(IBch.eq.10).or.(IBch.eq.12))
     -    tmpresult=tmpresult/9.


c... The photon multiplicity is given by the ratio of the squared
c... amplitudes, times a phase space factor:

      tmpresult = (alphem/pi)*m0**2/sqrt(1-mf**2/m0**2)
     -            *(tmpresult/msq2bodyds)


c... Finally, subtract the model-independent part which is contained in the
c... NLO-corrected two-body rate (see 1510.02473), and which should already have
c... been taken into account by Pythia

c      if (IBch.eq.5.or.IBch.eq.6.or.                ! previous version: do this
                                                     ! only for those channels for
c     -    IBch.eq.9.or.IBch.eq.11.or.IBch.eq.12)    ! which Pythia runs exist.
     
      tmpresult = tmpresult - dsIBfsrdxdy(IBch,x,y)

      dsIBffdxdy=tmpresult

      return
      end

