**********************************************************************
*** Subroutine dsanyieldgammaline returns observables related to 
*** gamma ray lines (or almost lines).
*** Output
***   n =  number of gamma lines
***   nbr[n] = N_gamma * branching fraction to this channel
***              [unitless]
***   eline[n] = energy of line [GeV]
***   wline[n] = width of line [GeV]
***   pdg_second[n] = PDG code of second particle
*** NOTE: nbr, eline and wline have to be defined as real*8 with
*** dimension 10, i.e. real*8 nbr(10),eline(10),wline(10)
*** pdg_second has to be defined as integer with dimension 10.
*** Author: Joakim Edsjo, edsjo@fysik.su.se
*** Date: May, 2014
**********************************************************************

      subroutine dsanyieldgammaline(n,nbr,eline,wline,pdg_second)
      implicit none
      include 'dsmssm.h'
      integer nmax
      parameter(nmax=10)
      integer n
      real*8 nbr(nmax),eline(nmax),wline(nmax)
      integer pdg_second(nmax)
      real*8 dssigmav0,dssigmav0tot,dsmwimp
      real*8 ngaga,ngaz,mwimp,mp2
      real*8 sigv0

      sigv0=dssigmav0tot()

      ngaga=2.d0*dssigmav0(22,22)/sigv0
      ngaz=dssigmav0(22,23)/sigv0

      mwimp=dsmwimp()
c...  Add gamma gamma to array
      n=1
      nbr(1)=ngaga
      eline(1)=mwimp
      wline(1)=0.d0
      pdg_second(1)=22 ! photon

c... Add gamma z if non-zero
      if (ngaz.gt.0.d0) then
         n=2
         mp2=mass(kz)
         nbr(n)=ngaz
         eline(n)=mwimp-m2**2/(4.0d0*mwimp)
         wline(n)=width(kz)
         pdg_second(n)=23 ! Z0
      endif

      return
      end
