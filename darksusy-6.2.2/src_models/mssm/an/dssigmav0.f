**********************************************************************
*** function dssigmav0 returns the annihilation cross section
*** sigma v at p=0 for neutralino-neutralino annihilation into 2-body
*** final states.
*** 
*** input: pdg1, pdg2 -- PDG codes of final state particles
***
*** Units of returned cross section: cm^3 s^-1
***
*** author: Torsten.Bringmann.fys.uio.no
*** (based on old dssigmav with channel codes)
*** date: 2014-11-14
**********************************************************************

      real*8 function dssigmav0(pdg1,pdg2)
      implicit none
      include 'dsmssm.h'
      include 'dsanyieldmodelcom.h'
      include 'dsandwcom.h'
      include 'dsidtag.h'
      include 'dssvcom.h'
      include 'dsmpconst.h'

c... FIXME: add error handling !?

      integer pdg1,pdg2,p1,p2,partch
      real*8 dsandwdcosnn, tmp
      integer i

      integer dsidnumber
      integer idold
      data idold/-123456789/
      save idold

      dssigmav0=0.d0

      if (abs(pdg1).lt.abs(pdg2)) then
        p1=pdg1
        p2=pdg2
      else
        p1=pdg2
        p2=pdg1        
      endif

       if (p1.eq.25.and.p2.eq.25) then
        partch=3
      elseif (p1.eq.25.and.p2.eq.35) then
        partch=2
      elseif (p1.eq.35.and.p2.eq.35) then
        partch=1
      elseif (p1.eq.36.and.p2.eq.36) then
        partch=4
      elseif (p1.eq.25.and.p2.eq.36) then
        partch=6
      elseif (p1.eq.35.and.p2.eq.36) then
        partch=5
      elseif (abs(p1).eq.37.and.p1.eq.-p2) then
        partch=7
      elseif (p1.eq.23.and.p2.eq.25) then
        partch=9
      elseif (p1.eq.23.and.p2.eq.35) then
        partch=8
      elseif (p1.eq.23.and.p2.eq.36) then
        partch=10
      elseif (abs(p1).eq.24.and.abs(p2).eq.37.and.abs(p1+p2).eq.13) then
        partch=11
      elseif (p1.eq.23.and.p2.eq.23) then
        partch=12 ! Z0 Z0
      elseif (abs(p1).eq.24.and.p2.eq.-p1) then
        partch=13 ! W+ W-
      elseif (abs(p1).eq.12.and.p2.eq.-p1) then
        partch=14 ! nu_e nu_e-bar
      elseif (abs(p1).eq.11.and.p2.eq.-p1) then
        partch=15 ! e- e+
      elseif (abs(p1).eq.14.and.p2.eq.-p1) then
        partch=16 ! nu_mu nu_mu-bar
      elseif (abs(p1).eq.13.and.p2.eq.-p1) then
        partch=17 ! mu- mu+
      elseif (abs(p1).eq.16.and.p2.eq.-p1) then
        partch=18 ! nu_tau nu_tau-bar
      elseif (abs(p1).eq.15.and.p2.eq.-p1) then
        partch=19 ! tau- tau+
      elseif (abs(p1).eq.2.and.p2.eq.-p1) then
        partch=20 ! u u-bar
      elseif (abs(p1).eq.1.and.p2.eq.-p1) then
        partch=21 ! d d-bar
      elseif (abs(p1).eq.4.and.p2.eq.-p1) then
        partch=22 ! c c-bar
      elseif (abs(p1).eq.3.and.p2.eq.-p1) then
        partch=23 ! s s-bar
      elseif (abs(p1).eq.6.and.p2.eq.-p1) then
        partch=24 ! t t-bar
      elseif (abs(p1).eq.5.and.p2.eq.-p1) then
        partch=25 ! b b-bar
      elseif (p1.eq.21.and.p2.eq.21) then
        partch=26 ! gluon gluon (1-loop)
      elseif (p1.eq.22.and.p2.eq.22) then
        partch=28  ! gamma gamma 
      elseif (p1.eq.22.and.p2.eq.23) then
        partch=29  ! gamma Z
      elseif (p1.eq.10000.and.p2.eq.10000) then
         partch=-10000          ! dummy chanel -- not used!!!!
         dssigmav0=0.d0
         return
      else
         write(*,*)'WARNING -- channel not implemented in ',
     &   'dssigmav0: pdg1 = ',
     &   pdg1,' pdg2 = ',pdg2   
        dssigmav0=0d0
        return
      endif



*** For reference, these are the old channel numbers that are used in the code
*** below
***
***   partch   process
***   ------   -------
***        0   All processes, i.e. the full annihilation cross section
***        1   H1 H1
***        2   H1 H2
***        3   H2 H2
***        4   H3 H3
***        5   H1 H3
***        6   H2 H3
***        7   H- H+
***        8   H1 Z
***        9   H2 Z
***       10   H3 Z
***       11   W- H+ and W+ H-
***       12   Z0 Z0
***       13   W+ W-
***       14   nu_e nu_e-bar
***       15   e+ e-
***       16   nu_mu nu_mu-bar
***       17   mu+ mu-
***       18   nu_tau nu_tau-bar
***       19   tau+ tau-
***       20   u u-bar
***       21   d d-bar
***       22   c c-bar
***       23   s s-bar
***       24   t t-bar
***       25   b b-bar
***       26   gluon gluon (1-loop)
***       27   q q gluon (kept for "historic" reasons, even if not 2-body)
***       28   gamma gamma (1-loop)
***       29   Z gamma (1-loop)


      
      if (idold.ne.dsidnumber()) then
         wtot=dsandwdcosnn(0.0d0,0.0d0,kn(1),kn(1))
         mx=mass(kn(1))
c... sigma v = w / (4*E_1^2) = wtot / (2*mx**2) since integration
c... over cos theta gives factor of 2.
         sigmav = gev2cm3s*wtot/(2.d0*mx**2) ! in cm^3/s
         do i=1,29
            sigv(i) = prtial(i)/wtot*sigmav  ! NB: this even calculates 27=qqg
         enddo
         idold=dsidnumber()
      endif
      
      if (partch.gt.0) tmp=sigv(partch)
      
      if (sin(tmp).ne.sin(tmp).or.tmp.lt.0d0) then
        write(*,*) 
     &     'ERROR in dssigmav0: cross section = ',tmp,' for ',pdg1,pdg2
        tmp=0d0
      endif

      dssigmav0=tmp
      
      return

      end
