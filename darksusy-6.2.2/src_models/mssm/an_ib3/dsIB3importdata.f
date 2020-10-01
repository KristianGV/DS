**********************************************************************
*** subroutine dsIB3importdata reads data files that contain tabulated
*** yields from qqg final states and stores information in common 
*** blocks.
***
***   input:   qch   -- quark channel (7-up to 12-bottom)
***            spectype -- see beginning of dsIB3yieldtab for definitions
***
*** author: torsten.bringmann@fys.uio.no, 2016-04-02
***********************************************************************

      subroutine dsIB3importdata(qch,spectype)
      implicit none
      integer qch, spectype

      include 'dsIB3com.h'
      include 'dsio.h'

   
      integer mmax, i, j, ienergy
      real*8 tmpdat(nmass), Rtmpdat(nmass)
      character filename*200, mystring*200, shortstr*100, Rpstring*1000

      if (qch.lt.7.or.qch.gt.18   ! 13-18 are 2-body channels (from Pythia 8 runs)
     &    .or.spectype.lt.1.or.spectype.gt.ntype) goto 1100

      fitexists(qch,spectype)=.false.
     
c... get system-specific location of dat/ directory
      call dsdatafile(filename,'')

c... add quark type to filename
      if (qch.eq.7) then
        filename=filename(1:index(filename,' ')-1)//'qqg_tables/uug_ '     
      elseif (qch.eq.8) then
        filename=filename(1:index(filename,' ')-1)//'qqg_tables/ddg_ '       
      elseif (qch.eq.9) then
        filename=filename(1:index(filename,' ')-1)//'qqg_tables/ccg_ '       
      elseif (qch.eq.10) then
        filename=filename(1:index(filename,' ')-1)//'qqg_tables/ssg_ '       
      elseif (qch.eq.11) then
        filename=filename(1:index(filename,' ')-1)//'qqg_tables/ttg_ '       
      elseif (qch.eq.12) then
        filename=filename(1:index(filename,' ')-1)//'qqg_tables/bbg_ '       
      elseif (qch.eq.13) then
        filename=filename(1:index(filename,' ')-1)//'qq_tables/uu_ '     
      elseif (qch.eq.14) then
        filename=filename(1:index(filename,' ')-1)//'qq_tables/dd_ '       
      elseif (qch.eq.15) then
        filename=filename(1:index(filename,' ')-1)//'qq_tables/cc_ '       
      elseif (qch.eq.16) then
        filename=filename(1:index(filename,' ')-1)//'qq_tables/ss_ '       
      elseif (qch.eq.17) then
        filename=filename(1:index(filename,' ')-1)//'qq_tables/tt_ '       
      elseif (qch.eq.18) then
        filename=filename(1:index(filename,' ')-1)//'qq_tables/bb_ '       
      endif
c... then add threebodytype to filename
      if (qch.le.12) then
        if (mod(spectype,4).eq.0) then
          filename=filename(1:index(filename,' ')-1)//'m_sq->Infty_unmixed_ '
        elseif (mod(spectype,4).eq.3) then
          filename=filename(1:index(filename,' ')-1)//'m_sq->Infty_mixed_ '
        elseif (mod(spectype,4).eq.2) then
          filename=filename(1:index(filename,' ')-1)//'VIB_unmixed_ '
        elseif (mod(spectype,4).eq.1) then
          filename=filename(1:index(filename,' ')-1)//'VIB_mixed_ '
        endif
      else
        if (prtlevel.gt.0) write(*,*) 'WARNING from dsIB3importdata: Assessing 2-body yields',
     &             ' rather than 3-body yields!'              
      endif
       
c... and finally the final state type
      if (((spectype-1)/8).eq.0.and.((spectype-1)/4).eq.0) then
        filename=filename(1:index(filename,' ')-1)//'gam_dNdx_int.dat'
      elseif (((spectype-1)/8).eq.0.and.((spectype-1)/4).le.1) then
        filename=filename(1:index(filename,' ')-1)//'gam_dNdx.dat'
      elseif (((spectype-1)/8).eq.1.and.((spectype-9)/4).eq.0) then
        filename=filename(1:index(filename,' ')-1)//'pbar_dNdT_int.dat'
      elseif (((spectype-1)/8).eq.1.and.((spectype-9)/4).le.1) then
        filename=filename(1:index(filename,' ')-1)//'pbar_dNdT.dat'
      else
        goto 1100
      endif

      if (prtlevel.gt.1) write (*,*) 'dsIB3importdata: Reading data from file ',filename
      open (unit=20,file=filename) 

c... Here, we read in the parameters that interpolate between the two extreme
c... 3-body spectra
      if (qch.gt.12) goto 50   ! does not apply to 2-body spectra

      read (20,*,end=1000,err=1000)
      read (20,*,end=1000,err=1000)
      read (20,'(A)',end=1000,err=1000) mystring
      if (mystring(1:40).ne.'# Function for determining full spectrum') goto 990

      read (20,'(A)',end=1000,err=1000) mystring
      if (mystring.eq.'# Y(R) = 1') then ! fit not supplied for this case

        goto 50 ! read in masses

      elseif (mystring.eq.'# Y(R) = R*10**(c1*R**n1 + c2*R**n2 + c3*R**n3)') then

        read (20,*,end=1000,err=1000)
        read (20,'(A)',end=1000,err=1000) mystring
        shortstr=mystring(index(mystring,'c1 =')+4:)
        read(shortstr(1:index(shortstr,',')-1),*,err=1000) ci(1,qch,spectype)  
        shortstr=shortstr(index(shortstr,'c2 =')+4:)
        read(shortstr(1:index(shortstr,',')-1),*,err=1000) ci(2,qch,spectype)   
        shortstr=shortstr(index(shortstr,'c3 =')+4:)
        read(shortstr(1:index(shortstr,',')-1),*,err=1000) ci(3,qch,spectype)
        shortstr=shortstr(index(shortstr,'n1 =')+4:)
        read(shortstr(1:index(shortstr,',')-1),*,err=1000) ni(1,qch,spectype)   
        shortstr=shortstr(index(shortstr,'n2 =')+4:)
        read(shortstr(1:index(shortstr,',')-1),*,err=1000) ni(2,qch,spectype)   
        shortstr=shortstr(index(shortstr,'n3 =')+4:)
        read(shortstr(1:index(shortstr,'.',BACK=.true.)-1),*,err=1000) ni(3,qch,spectype)
       
 20     read (20,'(A)',end=1000,err=1000) mystring
        if (mystring.ne.'# R''') goto 20
        read (20,*,end=1000,err=1000)
        read (20,'(A)',end=1000,err=1000) Rpstring

        fitexists(qch,spectype)=.true.
       
      else ! '# Y(R) = ...' line not recognized

        goto 990

      endif


c... continue to the point where we read in DM masses
 50   read (20,'(A)',end=1000,err=1000) mystring
      if (mystring.ne.'# No. of masses') goto 50


      read (20,*,end=1000,err=1000)
      read (20,*,end=1000,err=1000) mmax
      read (20,*,end=1000,err=1000)
      read (20,'(A)',end=1000,err=1000) mystring

      if (mystring.ne.'# m0_i (GeV)') goto 990
      if (mmax.gt.nmass) then
        if (prtlevel.gt.1) write(*,*) 'WARNING in dsIB3importdata: ',
     &                   'more masses tabulated than allowed by array size!'
        read (20,*,end=1000,err=1000) tmpdat(1:nmass)
        if (prtlevel.gt.1) write(*,*) 'Only importing data up to m_dm = ', tmpdat(nmass)
        mmax=nmass
      else
        read (20,*,end=1000,err=1000) tmpdat(1:mmax)
      endif
      do i = 1, nmass
        if (i.le.mmax) fitmass(qch,spectype,i)=tmpdat(i)
        if (i.gt.mmax) fitmass(qch,spectype,i)=tmpdat(mmax)*1.d2 ! fill up with some large mass value
      enddo
      if (fitexists(qch,spectype)) then ! do the same with the Rprime values
        read (Rpstring,*,end=1000,err=1000) Rtmpdat(1:mmax)  
        do i = 1, nmass
          if (i.le.mmax) Rprime(qch,spectype,i)=Rtmpdat(i)
          if (i.gt.mmax) Rprime(qch,spectype,i)=Rtmpdat(mmax)
        enddo
      endif

      mlo(qch,spectype)=1
      mhi(qch,spectype)=mmax


c... read in spectra
      read (20,*,end=1000,err=1000)
      read (20,*,end=1000,err=1000)
      read (20,*,end=1000,err=1000)
      read (20,'(A)',end=1000,err=1000) mystring


      ienergy = 1
 100  read (20,*,end=200,err=1000) xdat(qch,spectype,ienergy),
     &                              tmpdat(1:mmax)
      do i = 1, nmass
        if (i.le.mmax) dNdEdat(qch,spectype,i,ienergy)=tmpdat(i)
        if (i.gt.mmax) dNdEdat(qch,spectype,i,ienergy)=1.0d-50    !tmpdat(mmax)
c FIXME Maybe add more sophisticated extrapolation -> depends on spectype!
        
        if (dNdEdat(qch,spectype,i,ienergy).eq.0.d0)
     &      dNdEdat(qch,spectype,i,ienergy) = 1.0d-50 ! avoid zero for log-interpolation    
      enddo
      read (20,'(A)',end=200,err=1000) mystring
      ienergy = ienergy + 1
c...  print warning and exit loop for too large ienergy!
      if (ienergy.gt.nenergy) then
        if (prtlevel.gt.1) write(*,*) 'WARNING in dsIB3importdata: more energy bins',
     &                   ' tabulated than allowed by array size!'
        if (prtlevel.gt.1) write(*,*) 'Only importing data up to TGeV = ', 
     &              xdat(qch,spectype,ienergy-1)
        goto 200
      endif
      goto 100
 200  continue  
      ienergy = ienergy - 1
      do i = ienergy + 1, nenergy
        xdat(qch,spectype,i) = xdat(qch,spectype,ienergy)
        do j= 1, nmass
          dNdEdat(qch,spectype,j,i) = 1.0d-50 ! avoid zero for log-interpolation
        enddo
      enddo
      elo(qch,spectype)=1
      ehi(qch,spectype)=ienergy
      
      close (20)
 
      return
 
  990 if (prtlevel.gt.0) write(*,*) 'Unexpected line found when importing data: ', mystring
      return         
 1000 if (prtlevel.gt.0)write(*,*) 'ERROR in dsIB3importdata when reading file ', filename
      write(*,*)
      return
 1100 if (prtlevel.gt.0) write(*,*) 'ERROR in dsIB3importdata -- requested combination ', 
     &           'of quark channel and spectrum type not implemented:',
     &           qch,spectype
      write(*,*)
      return     
      end
