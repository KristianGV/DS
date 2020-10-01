      subroutine dsrdset(key,value)
c...set parameters for relic density routines
c...  key - character string
c...  value - character string
c...parameters implemented:
c...  key='dof'
c...    value='help' - show a list of possible values
c...    value='default' - default
c...    value='1' - Gelmini-Gondolo 150MeV
c...    value='2' - Hindmarsch-Philipsen HP-A
c...    value='3' - Hindmarsch-Philipsen HP-B (default)
c...    value='4' - Hindmarsch-Philipsen HP-B2
c...    value='5' - Hindmarsch-Philipsen HP-B3
c...    value='6' - Hindmarsch-Philipsen HP-C
c...    value='7' - Drees-Hajkarim-Schmitz, 1503.03513
c...  key='help'
c...    any value - show a list of possible keys
c...author: paolo gondolo 2007-12-29
c...modified: 2017-05-12 added Drees et al [torsten bringmann]
c...modified: 2018-01-18 directly read in dof tables because they
c...                     are also needed elsewhere [tb]
      implicit none
      include 'dsrdcom.h'
      character*(*) key,value
      character*200 filename


c...print list of keys
      if (key.eq.'help') then
         write (*,*) 'dsrdset: allowed keys are'
         write (*,*) '	dof	-- degrees of freedom in early universe'

cc...default values
ccc [the lines below are a shortcut, but an abuse of the key/value principle...]
c      else if (key.eq.'default') then
c         if (dofcode.ne.1) rdinit=0
c         dofcode=7  ! DS before 6.0: default was dofcode=3

c...degrees of freedom
      else if (key.eq.'dof') then
         if (value.eq.'1') then 
            dofcode=1
            call dsdatafile(filename,'dsdofGG_150.dat')
            call dsrdreaddof(filename)
         else if (value.eq.'2') then
            dofcode=2
            call dsdatafile(filename,'dsdofHP_A.dat')
            call dsrdreaddof(filename)
         else if (value.eq.'3') then
            dofcode=3
            call dsdatafile(filename,'dsdofHP_B.dat')
            call dsrdreaddof(filename)
         else if (value.eq.'4') then
            dofcode=4
            call dsdatafile(filename,'dsdofHP_B2.dat')
            call dsrdreaddof(filename)
         else if (value.eq.'5') then
            dofcode=5
            call dsdatafile(filename,'dsdofHP_B3.dat')
            call dsrdreaddof(filename)
         else if (value.eq.'6') then
            dofcode=6
            call dsdatafile(filename,'dsdofHP_C.dat')
            call dsrdreaddof(filename)
         else if (value.eq.'7'.or.value.eq.'default') then
            dofcode=7
            call dsdatafile(filename,'dsdofDHS.dat')
            call dsrdreaddof(filename)
         else if (value.eq.'help') then
            write (*,*) 'dsrdset: key=dof can have values'
            write (*,*) '	''1''	Gelmini-Gondolo 150MeV'
            write (*,*) '	''2''	Hindmarsch-Philipsen HP-A'
            write (*,*) '	''3''	Hindmarsch-Philipsen HP-B'
            write (*,*) '	''4''	Hindmarsch-Philipsen HP-B2'
            write (*,*) '	''5''	Hindmarsch-Philipsen HP-B3'
            write (*,*) '	''6''	Hindmarsch-Philipsen HP-C'
            write (*,*) '	''7''	Drees-Hajkarim-Schmitz (default)'
         else if (value.eq.'show') then
            write (*,*) 'dsrdset: key=dof has value ',dofcode
         else
            goto 1000
       endif

c...invalid choice
      else
         goto 1000
      endif

      return

 1000 continue
      write (*,*) 'dsrdset: unrecognized key and/or invalid value ',
     &     key,value
      stop
      end
