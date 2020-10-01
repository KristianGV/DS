*******************************************************************************
***  subroutine dsgivemodel_decayingDM reads in the parameters to describe  ***
***  a generic decaying DM particle and transfers them to common blocks     ***
***                                                                         ***
***  input:                                                                 ***
***                                                                         ***
***    mDM        - DM mass (in GeV)                                        ***
***    DecRate    - Total decay rate Gamma_tot (in 1/s)                     ***
***    n          - number of (2-body) decay channels                       ***
***    BR         - array of size n containg the branching ratios           ***
***                 Gamma_i/Gamma_tot                                       *** 
***    PDG1, PDG2 - array of size n with PDG code for final state particles ***
***                 (e.g. 5 for bbar, 24 for W^+W^-)                        ***
***                                                                         ***
*** author: Torsten.Bringmann.fys.uio.no                                    ***
*** date 2016-02-06                                                         ***
*******************************************************************************
      subroutine dsgivemodel_decayingDM(mDM,DecRate,n,BR,PDG1,PDG2)
      implicit none
      include 'dsgeneric_decayingDM.h'

      real*8 mDM,DecRate
      integer n, i, j
      integer PDG1(n),PDG2(n)
      real*8 BR(n)
      logical channelfound
c-----------------------------------------------------------------------


      mass(kdm) = mDM
      spin(kdm) = 0 ! for scalar DM 
      kdof(kdm) = 1 ! for scalar DM 

c... write decay rate to common block
      Gammatot = DecRate

c... write decay channels to common block
      do i=1,n
        channelfound = .false.
        j=0
        do 10 while (j.le.numdecch2b-1)
          j=j+1
          if ((PDG1(i).eq.dec_2body(j,1).and.PDG2(i).eq.dec_2body(j,2)).or.
     &        (PDG1(i).eq.dec_2body(j,2).and.PDG2(i).eq.dec_2body(j,1))) then
         
            decBR(j) = BR(i)
            channelfound = .true.
          endif
 10     continue  
        if (.not.channelfound) 
     &       write(*,*) 'dsgivemodel_decayingDM -- the following final states ',
     &                'are not implemented and will be treated as invisible: ',
     &                  PDG1(i),PDG2(i)
      enddo

      end


