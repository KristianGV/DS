*******************************************************************************
*** Subroutine dsddDMCRquenching_set is used to steer which quenching       ***
*** factors should be used for the dsddDMCR routines. The default is no     ***
*** no quenching, i.e. the energy range that dsddDMCRcountrate receives as  ***
*** input refers directly to the nuclear recoil energy.                     ***
*** For other implemented options, see below. Note that most just rely on   ***
*** tabulated queching factors, which makes it straightforward to add new   ***
*** ones (by replacing this function and dsddDMCRquenching).                ***
***                                                                         ***
***                                                                         ***
*** author: Torsten.Bringmann.fys.uio.no                                    ***
*** date 2018-07-04                                                         ***
*******************************************************************************
      subroutine dsddDMCRquenching_set(option)
      implicit none
      include 'dsddcom.h'
      include 'dsio.h'
      character*(*) option
      character*200 filename1,filename2
   
      quenching_set=.false.
      
      if (option.eq.'default'.or.option.eq.'no_quenching') then
         quenchhow = 1
      
      elseif (option.eq.'borexino') then 
         quenchhow = 2 ! read in from table
         call dsdatafile(filename1,'Quench_Borexino.dat')
         call dsdatafile(filename2,'Quench_Borexino_dTdTQ.dat')
         call dsddreadQF(filename1,filename2)
         quenching_set=.true.
      else
         if (prtlevel.ge.1) then
            write(*,*) 'WARNING: dsddDMCRquenching_set has been called'
            write(*,*) 'with an unsupported option.'
            write(*,*) 'Will not change setting from quenchhow = ',quenchhow
         endif
      endif

      return
      end
