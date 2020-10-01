*******************************************************************************
*** This routine orders all sfermions according to flavour. For sleptons,   ***
*** e.g. this implies the following for the internal particle code ksl_flav ***
*** (while ksl itself is not necessarily ordered):                          ***
***                                                                         ***
***       ksl_flav(1,1) = lightest selectron-like sfermion                  ***
***       ksl_flav(2,1) = lightest smuon-like sfermion                      ***
***       ksl_flav(3,1) = lightest stau-like sfermion                       ***
***       ksl_flav(1,2) = heavier selectron-like sfermion                   ***
***       ksl_flav(2,2) = heavier smuon-like sfermion                       ***
***       ksl_flav(3,2) = heavier stau-like sfermion                        ***
***                                                                         ***
*** Here, the two "most selectron-like" particles are defined as those that ***
*** mix strongest with the light- and right-handed selectron, respectively. ***
*** For squarks, the analogous applies, while for sneutrinos only the first ***
*** 3 entries in the above table exist.                                     ***
***                                                                         ***
***       output: ksnu_flav, ksl_flav, ksqu_flav, ksqd_flav                 ***
***               (stored as common block variables in dsmssm.h)            ***
***                                                                         ***
*** NB: This routine assumes that the mixing and mass common blocks         ***
***     have already been set.                                              ***
***                                                                         ***
*** date:   2016-11-15                                                      ***
*** author: torsten.bringmannl@fys.uio.no,                                  ***
*** modified: tb 2018-12-07, added widths                                   ***
*** modified: tb 2018-12-07, introduced new particle codes ksl_flav etc.    ***
***                          rather than changing any masses or widths      ***
*******************************************************************************

      subroutine dsorder_flavour
      implicit none

      include 'dsmssm.h'
      
      integer i, ktmp, kflav(6)
c      real*8 widthtmp(6), masstmp(6), mixtmpl(6,3), mixtmpr(6,3)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc re-order down-type sleptons as i -> kflav(i)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do i=1,6
        kflav(i)=1
      enddo  

      do i=2,6                                                        ! find the slepton that
        if (abs(sldlmx(i,1)).gt.abs(sldlmx(kflav(1),1))) kflav(1) = i ! looks most like a left-
      enddo                                                           ! handed selectron
      do i=2,6                                               
        if ((i.ne.kflav(1)).and.                           ! find the (remaining!) slepton 
     &      (abs(sldrmx(i,1)).gt.abs(sldrmx(kflav(4),1)))) ! that looks most like a right- 
     &      kflav(4) = i                                   ! handed selectron
      enddo                                                          

      do i=2,6 ! now repeat this for smuons
        if ((i.ne.kflav(1).and.i.ne.kflav(4)).and.                         
     &     (abs(sldlmx(i,2)).gt.abs(sldlmx(kflav(2),2)))) kflav(2) = i 
      enddo                                                           
      do i=2,6                                               
        if ((i.ne.kflav(1).and.i.ne.kflav(4).and.i.ne.kflav(2)).and.                            
     &      (abs(sldrmx(i,2)).gt.abs(sldrmx(kflav(5),2)))) kflav(5) = i                                   
      enddo                                                          

      i=0 ! the remaining ones must be stau-like... 
 10   i = i+1
      if (i.eq.kflav(1).or.i.eq.kflav(2).or.i.eq.kflav(4).or.
     &      i.eq.kflav(5)) goto 10
      kflav(3) = i
 20   i = i+1
      if (i.eq.kflav(1).or.i.eq.kflav(2).or.i.eq.kflav(3).or.
     &      i.eq.kflav(4).or.i.eq.kflav(5)) goto 20 
      kflav(6) = i     

      do i=1,3 ! change mass order if necessary
        if (mass(ksl(kflav(i+3))).lt.mass(ksl(kflav(i)))) then 
          ktmp = kflav(i)
          kflav(i) = kflav(i+3)
          kflav(i+3) = ktmp
        endif
      enddo

c... now write this translation from old to new ordering to common blocks
      do i=1,3
        ksl_flav(i,1) = ksl(kflav(i))
        ksl_flav(i,2) = ksl(kflav(3+i))
      enddo
 
c      write(*,*) 'Sleptons (unordered, ordered): '
c      write(*,*) ksl
c      write(*,*) ksl_flav
 
c... now we have to update the masses and mixing matrices to these new conventions:
c... OBSOLETE -- this only happened in a previous version

c      do i=1,6
c        masstmp(i) = mass(ksl(i)) 
c        widthtmp(i) = width(ksl(i)) 
c        do j=1,3
c          mixtmpl(i,j) = sldlmx(i,j)
c          mixtmpr(i,j) = sldrmx(i,j)
c        enddo
c      enddo
c      do i=1,6         
c        mass(ksl(i)) = masstmp(kflav(i))
c        width(ksl(i)) = widthtmp(kflav(i))
c        do j=1,3
c          sldlmx(i,j) = mixtmpl(kflav(i),j)
c          sldrmx(i,j) = mixtmpr(kflav(i),j)
c        enddo
c      enddo
  
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc re-order sneutrinos as i -> kflav(i)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do i=1,3
        kflav(i)=1
      enddo  

      do i=2,3 ! find the sneutrino that looks most like a snu_e
        if (abs(slulmx(i,1)).gt.abs(slulmx(kflav(1),1))) kflav(1) = i
      enddo  
      do i=2,3 ! find the sneutrino that looks most like a snu_mu
        if ((i.ne.kflav(1)).and.                         
     &     (abs(slulmx(i,2)).gt.abs(slulmx(kflav(2),2)))) kflav(2) = i 
      enddo                                                           
      do i=2,3 ! find the sneutrino that looks most like a snu_tau 
        if (i.ne.kflav(1).and.i.ne.kflav(2)) kflav(3) = i  
      enddo
c... now write this translation from old to new ordering to common blocks
      do i=1,3
        ksnu_flav(i,1) = ksnu(kflav(i))
      enddo

c      write(*,*) 'Sneutrinos (unordered, ordered): '
c      write(*,*) ksnu
c      write(*,*) ksnu_flav

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc re-order down-type squarks as i -> kflav(i)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do i=1,6
        kflav(i)=1
      enddo  

      do i=2,6                                                        ! find the squark that
        if (abs(sqdlmx(i,1)).gt.abs(sqdlmx(kflav(1),1))) kflav(1) = i ! looks most like a left-
      enddo                                                           ! handed down-squark
      do i=2,6                                               
        if ((i.ne.kflav(1)).and.                           ! find the (remaining!) squark 
     &      (abs(sqdrmx(i,1)).gt.abs(sqdrmx(kflav(4),1)))) ! that looks most like a right- 
     &      kflav(4) = i                                   ! handed down-squark
      enddo                                                          

      do i=2,6 ! now repeat this for s-strange
        if ((i.ne.kflav(1).and.i.ne.kflav(4)).and.                         
     &     (abs(sqdlmx(i,2)).gt.abs(sqdlmx(kflav(2),2)))) kflav(2) = i 
      enddo                                                           
      do i=2,6                                               
        if ((i.ne.kflav(1).and.i.ne.kflav(4).and.i.ne.kflav(2)).and.                            
     &      (abs(sqdrmx(i,2)).gt.abs(sqdrmx(kflav(5),2)))) kflav(5) = i                                   
      enddo                                                          

      i=0 ! the remaining ones must be sbottom-like...  
 30   i = i+1
      if (i.eq.kflav(1).or.i.eq.kflav(2).or.i.eq.kflav(4).or.
     &      i.eq.kflav(5)) goto 30
      kflav(3) = i
 40   i = i+1
      if (i.eq.kflav(1).or.i.eq.kflav(2).or.i.eq.kflav(3).or.
     &      i.eq.kflav(4).or.i.eq.kflav(5)) goto 40 
      kflav(6) = i     

      do i=1,3 ! change mass order if necessary
        if (mass(ksqd(kflav(i+3))).lt.mass(ksqd(kflav(i)))) then 
          ktmp = kflav(i)
          kflav(i) = kflav(i+3)
          kflav(i+3) = ktmp
        endif
      enddo

c... now write this translation from old to new ordering to common blocks
      do i=1,3
        ksqd_flav(i,1) = ksqd(kflav(i))
        ksqd_flav(i,2) = ksqd(kflav(3+i))
      enddo
 
c      write(*,*) 'Down-Squarks (unordered, ordered): '
c      write(*,*) ksqd
c      write(*,*) ksqd_flav
  
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc re-order up-type squarks as i -> kflav(i)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do i=1,6
        kflav(i)=1
      enddo  

      do i=2,6                                                        ! find the squark that
        if (abs(squlmx(i,1)).gt.abs(squlmx(kflav(1),1))) kflav(1) = i ! looks most like a left-
      enddo                                                           ! handed up-squark
      do i=2,6                                               
        if ((i.ne.kflav(1)).and.                           ! find the (remaining!) squark 
     &      (abs(squrmx(i,1)).gt.abs(squrmx(kflav(4),1)))) ! that looks most like a right- 
     &      kflav(4) = i                                   ! handed up-squark
      enddo                                                          

      do i=2,6 ! now repeat this for charm-squark
        if ((i.ne.kflav(1).and.i.ne.kflav(4)).and.                         
     &     (abs(squlmx(i,2)).gt.abs(squlmx(kflav(2),2)))) kflav(2) = i 
      enddo                                                           
      do i=2,6                                               
        if ((i.ne.kflav(1).and.i.ne.kflav(4).and.i.ne.kflav(2)).and.                            
     &      (abs(squrmx(i,2)).gt.abs(squrmx(kflav(5),2)))) kflav(5) = i                                   
      enddo                                                          

      i=0 ! the remaining ones must be stau-like...  
 50   i = i+1
      if (i.eq.kflav(1).or.i.eq.kflav(2).or.i.eq.kflav(4).or.
     &      i.eq.kflav(5)) goto 50
      kflav(3) = i
 60   i = i+1
      if (i.eq.kflav(1).or.i.eq.kflav(2).or.i.eq.kflav(3).or.
     &      i.eq.kflav(4).or.i.eq.kflav(5)) goto 60 
      kflav(6) = i     

      do i=1,3 ! change mass order if necessary
        if (mass(ksqu(kflav(i+3))).lt.mass(ksqu(kflav(i)))) then 
          ktmp = kflav(i)
          kflav(i) = kflav(i+3)
          kflav(i+3) = ktmp
        endif
      enddo
 
c... now write this translation from old to new ordering to common blocks
      do i=1,3
        ksqu_flav(i,1) = ksqu(kflav(i))
        ksqu_flav(i,2) = ksqu(kflav(3+i))
      enddo
  
c      write(*,*) 'Up-Squarks (unordered, ordered): '
c      write(*,*) ksqu
c      write(*,*) ksqu_flav

      return
      end





        

