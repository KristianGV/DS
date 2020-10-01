*****************************************************************************
***   subroutine dsanyieldsetup prepares the common blocks with annihilation
***   channel branching ratios and Higgs decay widths for annihilation 
***   yield calculations. These common block variables are used by
***   by both the an_yield and the se_yield routines.
***   Note: This routine should be called whenever a new model is generated,
***         hence typically from dsmodelsetup
***   up things for  cascade decays.
***   This routine is the interface between SUSY and the halo
***   and Earth/Sun annihilation routines.
***
*** Author: Joakim Edsjo, edsjo@fysik.su.se
*** Date: 08-01-15
*** Modified: December, 2014
*** modified: torsten bringmann, Dec 2013: removed some common block variables
*****************************************************************************

      subroutine dsanyieldset
      implicit none
      include 'dsmssm.h'
      include 'dsanyieldmodelcom.h'
      include 'dsidtag.h'



      real*8 dssigmav0tot,dsmwimp,br,sv,mx
      integer i,j,kh(4)

c----------------------------------------------- set-up common variables


c...Set up Higgses (these are just juse below in loops)
      kh(1)=kh1
      kh(2)=kh2
      kh(3)=kh3
      kh(4)=khc


c...make sure we computed the branching ratios
      sv=dssigmav0tot()         ! JE CHECK. Is this needed
      mx=dsmwimp()

c...Transfer Higgs widths
c...Neutral Higgses
      do i=1,3
         do j=1,29
            br=hdwidth(j,i)/width(kh(i))
            if (br.gt.ansbrmin) then
               ans0br(j,i)=br
            else
               ans0br(j,i)=0.d0
            endif
         enddo
      enddo

      do i=1,3
         ans0m(i)=mass(kh(i))
      enddo

c...Charged Higgses
      do j=1,15
         br=hdwidth(j,4)/width(kh(4))
         if (br.gt.ansbrmin) then
            anscbr(j)=br
         else
            anscbr(j)=0.d0
         endif
      enddo

      anscm=mass(kh(4))

      end



















