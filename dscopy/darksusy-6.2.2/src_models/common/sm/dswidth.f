***********************************************************************
*** This function returns the (SM) width of standard model particles, 
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

      real*8 function dswidth(kPDG)
      implicit none
      integer kPDG, absk

      include 'dsparticles.h'
      include 'dssm.h'

      dswidth=0d0
      absk = abs(kPDG)

      if (absk.eq.1) then
        dswidth=width(kd)
      elseif (absk.eq.2) then
        dswidth=width(ku)
      elseif (absk.eq.3) then
        dswidth=width(ks)
      elseif (absk.eq.4) then
        dswidth=width(kc)
      elseif (absk.eq.5) then
        dswidth=width(kb)
      elseif (absk.eq.6) then
        dswidth=width(kt)
      elseif (absk.eq.11) then
        dswidth=width(ke)
      elseif (absk.eq.12) then
        dswidth=width(knue)
      elseif (absk.eq.13) then
        dswidth=width(kmu)
      elseif (absk.eq.14) then
        dswidth=width(knumu)
      elseif (absk.eq.15) then
        dswidth=width(ktau)
      elseif (absk.eq.16) then
        dswidth=width(knutau)
      elseif (absk.eq.21) then
        dswidth=width(kgluon)
      elseif (absk.eq.22) then
        dswidth=width(kgamma)
      elseif (absk.eq.23) then
        dswidth=width(kz)
      elseif (absk.eq.24) then
        dswidth=width(kw)
      elseif (absk.eq.25) then
        dswidth=width(khsm)
      else
        write(*,*) 'ERROR in dswidth: ',kPDG,' is not a valid particle code!'
        stop
      endif

      return
      end
