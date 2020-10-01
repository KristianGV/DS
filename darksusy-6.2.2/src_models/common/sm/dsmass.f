***********************************************************************
*** This function returns the mass of standard model particles, 
*** identified by their PDG code 
*** (http://pdg.lbl.gov/2007/reviews/montecarlorpp.pdf)           
***
***      kPDG  | particle 
***    --------+----------  
***        1   | u quark
***        2   | d quark
***        3   | s quark
***        4   | c quark
***        5   | b quark
***        6   | t quark
***       11   | electron 
***       12   | nu_e
***       13   | muon
***       14   | nu_mu
***       15   | tau
***       16   | nu_tau
***       21   | gluon
***       22   | photon
***       23   | Z
***       24   | W
***       25   | h
***
*** author: torsten.bringmann@fys.uio.no, 2014-05-09
***********************************************************************

      real*8 function dsmass(kPDG)
      implicit none
      integer kPDG, absk

      include 'dsparticles.h'
      include 'dssm.h'

      dsmass=0d0
      absk = abs(kPDG)

      if (absk.eq.1) then
        dsmass=mass(kd)
      elseif (absk.eq.2) then
        dsmass=mass(ku)
      elseif (absk.eq.3) then
        dsmass=mass(ks)
      elseif (absk.eq.4) then
        dsmass=mass(kc)
      elseif (absk.eq.5) then
        dsmass=mass(kb)
      elseif (absk.eq.6) then
        dsmass=mass(kt)
      elseif (absk.eq.11) then
        dsmass=mass(ke)
      elseif (absk.eq.12) then
        dsmass=mass(knue)
      elseif (absk.eq.13) then
        dsmass=mass(kmu)
      elseif (absk.eq.14) then
        dsmass=mass(knumu)
      elseif (absk.eq.15) then
        dsmass=mass(ktau)
      elseif (absk.eq.16) then
        dsmass=mass(knutau)
      elseif (absk.eq.21) then
        dsmass=mass(kgluon)
      elseif (absk.eq.22) then
        dsmass=mass(kgamma)
      elseif (absk.eq.23) then
        dsmass=mass(kz)
      elseif (absk.eq.24) then
        dsmass=mass(kw)
      elseif (absk.eq.25) then
        dsmass=mass(khsm)
      else
        write(*,*) 'ERROR in dsmass: ',kPDG,' is not a valid particle code!'
        stop
      endif

      return
      end
