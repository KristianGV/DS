***********************************************************************
*** The routine dskdset has to be called once before any microhalo (mh)
*** (mh) calculations; it makes all necessary initializations and sets 
*** some default settings.
***
*** author: torsten bringmann (troms@physto.se), 2010-01-23
*** updates: 2013-06-12 (removed model-dependence)
***          2017-04-28 (added two options for quark scattering)
***          2017-05-11 (changed d.o.f. treatment - unified with RD) 
***********************************************************************

      subroutine dskdset(key)
      implicit none
      character*(*) key

      include 'dskdcom.h'
      include 'dsrdcom.h'
      include 'dsio.h'

c      character*200 doffile

      mheps=1D-2   ! this sets the accuracy for the numerical
                   ! routines that determine Tkd

c... previuously, d.o.f were read in here...
c      call dsdatafile(doffile,'dskddof.dat')
c      call dskdreaddof(doffile)  ! reads in tabulated expressions
c                                 ! for eff. degrees of freedom

c...make sure that d.o.f. tables are read in
      if (rdinit.ne.1234) then
         if (prtlevel.gt.0) then 
            call dswrite(0,0,
     &         'dskdset: warning: no previous call to dsrdset')
            call dswrite(0,0,
     &         'dskdset: Choosing default version of d.o.f. tables...')
         endif
         call dsrdset('dof','default')
      endif
      if (dofcode.lt.7) then
         if (prtlevel.gt.0) then 
            call dswrite(0,0,
     &         'dskdset: warning: you called dsrdset with an option for dofcode')
            call dswrite(0,0,
     &         ' that does not provide geff. Hence switching back to default...')
         endif
        call dsrdset('dof','default')
      endif

      
c...how to deal with quark scattering above T_QCD
      if (key.eq.''.or.key.eq.'default'.or.key.eq.'qmax') then   ! procedure adopted in 1205.1914   
         quark_how = 1
      elseif (key.eq.'qmin') then  ! procedure adopted in 0903.0189 
         quark_how = 2
      else
         write (*,*) 'ERROR in dskdset: unrecognized input key: ', key
         stop
      endif

      return
      end





        

