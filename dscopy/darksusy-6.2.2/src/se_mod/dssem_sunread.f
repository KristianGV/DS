      subroutine dssem_sunread

***********************************************************************
*** Reads in data about the solar model used and stores it in a
*** common block (as described in dssem_sun.h).
*** Author: Joakim Edsjo
*** Date: 2003-11-25
*** Modified: 2004-01-28 (calculates potential instead of reading file)
***           2009-11-10 Added possibility to read solar models from
***             Serenelli et al.
***********************************************************************

      implicit none
      include 'dssem_sun.h'
      include 'dssecom.h'
      include 'dsmpconst.h'
      include 'dsio.h'

      real*8 dssem_sunpotint,dssem_suncdensint
      integer i,l,m


      character*200 scr

      real*8 tmp1,tmp2,tmp3,totfr,totfrheavy
      integer ztmp,itmp
      real*8 atmp,aph,apherr,aci,acierr,ab,abe

      real*8 mfrac,amass,spin
      character*4 cel

c...If already initialized, return
      if (sdread) then
        return
      endif

      sdread=.true.

c...First clear data file with mass fractions (in case we have already 
c...read a model previously)
c      write(*,*) 'dssem_sunread: clearing data in memory...'
      do i=1,sdmax
         do m=0,isomax
            do l=1,zmax
               sdmfr(l,m,i)=0.d0
            enddo
         enddo
      enddo

      do m=0,isomax
         do l=1,zmax
            sdaa(l,m)=0.d0
            sdma(l,m)=0.d0
            sdfr(l,m)=0.d0
            sdsp(l,m)=0.d0
         enddo
      enddo

      do l=1,zmax
         sdabund(l)=0.d0
      enddo
c      write(*,*) 'done.'

c----------------------------------------------------------------------
c--------- Bahcall style tables (sunfiletype=1)
c----------------------------------------------------------------------

      if (sunfiletype.eq.1) then ! Bahcall style data files

c...First define the abundances of heavy elements not listed in the
c...Bahcall et al standard solar model file
c...The number fractions are given by n_i = n_H log10(A-12)
c...The abundances listed here are from N. Grevesse and A.J. Sauval,
c...Space Science Reviews 85 (1998) 161.
      sdabund(10)  = 8.08d0  ! Ne
      sdabund(11)  = 6.33d0  ! Na
      sdabund(12)  = 7.58d0  ! Mg
      sdabund(13) = 6.47d0  ! Al
      sdabund(14) = 7.55d0  ! Si
      sdabund(16) = 7.33d0  ! S
      sdabund(18) = 6.40d0  ! Ar
      sdabund(20) = 6.36d0  ! Ca
      sdabund(26) = 7.50d0  ! Fe
      sdabund(28) = 6.25d0  ! Ni
      

c...Mass numbers
      sdaa(1,1)  = 1.0d0  ! H
      sdaa(2,2)  = 4.0d0  ! He4
      sdaa(2,1)  = 3.0d0  ! He3
      sdaa(6,1)  = 12.0d0 ! C12
      sdaa(7,1)  = 14.0d0 ! N14
      sdaa(8,1)  = 16.0d0 ! O16
      sdaa(10,1)  = 20.0d0 ! Ne
      sdaa(11,1)  = 23.0d0 ! Na
      sdaa(12,1)  = 24.0d0 ! Mg
      sdaa(13,1) = 27.0d0 ! Al
      sdaa(14,1) = 28.0d0 ! Si
      sdaa(16,1) = 32.0d0 ! S
      sdaa(18,1) = 40.0d0 ! Ar
      sdaa(20,1) = 40.0d0 ! Ca
      sdaa(26,1) = 56.0d0 ! Fe
      sdaa(28,1) = 59.0d0 ! Ni

c...Masses
      do m=1,isomax
        do l=1,zmax
          sdma(l,m)=sdaa(l,m)*(m_p+m_n)/2.0d0
        enddo
      enddo

c...Now read in data from solar model data file
      write(*,*) 'dssem_sunread: Opening file ',sunfile
      open(unit=13,file=sunfile,
     &  form='formatted',status='old')
      sdn=2  ! will add first line later

c...Read in scratch lines, NOTE: this might have to change with newer files
      do i=1,23
        read(13,'(A)',end=110) scr
      enddo

c...Calculate the total mass fraction (relative to H) of elements
c...heavier than O. This is just for normalization of the fractions
c...as a function of radius below
      totfrheavy=0.0d0
      do i=10,zmax
        totfrheavy=totfrheavy+sdma(i,1)*10**(sdabund(i)-12.0d0)
      enddo

c...Read in table and calculate remaining mass fractions from
c...average abundances sdabund
 100  read(13,*,end=110) sdm(sdn),sdr(sdn),tmp1,sdrho(sdn),tmp2,tmp3,
     &  sdmfr(1,1,sdn),sdmfr(2,2,sdn),sdmfr(2,1,sdn),sdmfr(6,1,sdn),
     &  sdmfr(7,1,sdn),sdmfr(8,1,sdn)
      totfr=sdmfr(1,1,sdn)+sdmfr(2,1,sdn)+sdmfr(2,2,sdn)+sdmfr(6,1,sdn)
     &  +sdmfr(7,1,sdn)+sdmfr(8,1,sdn)
c...Add the heavy elements (>O16)
      do i=10,zmax
        sdmfr(i,1,sdn)=(1.0d0-totfr)*sdma(i,1)*10**(sdabund(i)-12.0d0)
     &    /totfrheavy      
      enddo
      sdn=sdn+1
      if (sdn.gt.sdmax-1) then  ! need one more entry for last line
        goto 110
      endif
      goto 100
 110  continue
      sdn=sdn-1
      close(13)

c...Now add the first line with r/r_sun=0
      sdm(1)=0.0d0
      sdr(1)=0.0d0
      sdrho(1)=sdrho(2)
      do m=0,isomax
         do l=1,zmax
           sdmfr(l,m,1)=sdmfr(l,m,2)
         enddo
      enddo

c...Now add the last line with r/r_sun=1
      sdn=sdn+1
      sdm(sdn)=1.0d0
      sdr(sdn)=1.0d0
      sdrho(sdn)=0.0d0
      do m=0,isomax
         do l=1,zmax
           sdmfr(l,m,sdn)=sdmfr(l,m,sdn-1)
         enddo
      enddo


c...Electron density
c...Now read in data from solar model electron density file
      if (prtlevel.ge.2) write(*,*) 'dssem_sunread: Opening file ',sunnefile
      open(unit=13,file=sunnefile,
     &  form='formatted',status='old')
      sdnne=2  ! will add first line later

c...Read in scratch lines, NOTE: this might have to change with newer files
      do i=1,6
        read(13,'(A)',end=210) scr
      enddo

c...Read in table
 200  read(13,*,end=210) sdrne(sdnne),sdne(sdnne)
      sdnne=sdnne+1
      if (sdnne.gt.sdmax) then  ! need one more entry for last line
        write(*,*) 'ERROR in dssem_sunread: array too small.'
        write(*,*) 'Increase sdmax.'
        goto 210
      endif
      goto 200
 210  continue
      sdnne=sdnne-1
      close(13)

c...Now add the first line with r/r_sun=0
      sdne(1)=sdne(2)
      sdrne(1)=0.0d0

c----------------------------------------------------------------------
c--------- Serenelli style tables (sunfiletype=2)
c----------------------------------------------------------------------

      elseif (sunfiletype.eq.2) then

c...We will here use several sources to obtain the solar abundances and
c...their distribution in the Sun:
c...  - we will use tables from Serenelli for radial profiles for
c...    all elements up to Ni
c...  - for heavier elements, we will use heavy element abundances from
c...    AGSS09 (table 1), scaling them to Fe abundances in
c...    the Serenelli tables
c...  - in Serenelli, isotope abundances are given up to Oxygen, above that
c...    we will use isotope abundance tables in AGSS09 (table 3) 
c...    We will use atomic masses from NIST and spins from webelelments.com
c...    (read in here, but not really used, as we use element specific spin-
c...    dependent scattering tables anyway).

c...First define the abundances of heavy elements not listed in the
c...Serenelli et al files. These are taken from Table 1 in AGSS09
c...(Annu. Rev. Astro. Astrophys. 47 (2009) 481, arXiv:0909.0948)
c...The reported values are photoshperic (or derived from meteorites). To
c...get bulk values, we correct with 0.05 dex for He and 0.04 dex for
c...heavier elements, according to Turcotte & Wimmer-Schweingruber (2002)

c...The number fractions are given by n_i = n_H log10(A-12)
c...Read in these ones from sunabundfile
      if (prtlevel.ge.2) write(*,*) 'dssem_sunread: Opening file ',sunabundfile
      open(unit=13,file=sunabundfile,status='old',form='formatted')

c...Read in scratch lines, NOTE: this might have to change with newer files
      do i=1,4
        read(13,'(A)',end=302) scr
      enddo
     
 301  read(13,*,end=302) ztmp,atmp,aph,apherr,aci,acierr

      if (absrc.eq.'ph') then ! photoshperic
        ab=aph
        abe=apherr
        if (ab.lt.-8.99d0) then ! not given in table
          ab=aci
          abe=acierr
        endif
      else ! meteoritic
        ab=aci
        abe=acierr
        if (ztmp.eq.1.or.ztmp.eq.2.or.ztmp.eq.18
     &     .or.ab.lt.-0.75d0) then ! not reliable
           ab=aph
           abe=apherr
        endif
      endif
      
      if (ztmp.eq.1) then
         sdabund(ztmp)=ab
      elseif (ztmp.eq.2) then
         sdabund(ztmp)=ab+0.05d0
      else
         sdabund(ztmp)=ab+0.04d0
      endif
      sdaberr(ztmp)=abe
c...JE FIXME: Should we get these from here, or aveage below from
c...isotope table instead?
      sdaa(ztmp,0)=atmp ! average mass number for all isotopes
      sdma(ztmp,0)=atmp*(m_p+m_n)/2.d0
      sdaa(ztmp,1)=atmp ! good first approximation, will be refined later JE FIXME
      sdma(ztmp,1)=atmp*(m_p+m_n)/2.d0

      goto 301
 302  close(13)

c...JE FIXME: As an alternative, we could Read in these from AGSS09-table3.txt
c...instead (why?), making sure to assign correct isotope values as those in
c...Serenelli. Where Serenelli
c...uses only one isotope, we should really distribute the mass fractions on
c...the available elements.

c...Read in the isotope abundances, masses and atomic numbers
c...from file
      if (prtlevel.ge.2) write(*,*) 'dssem_sunread: Opening file ',sunisofile
      open(unit=13,file=sunisofile,status='old',form='formatted')

c...Read in scratch lines, NOTE: this might have to change with newer files
      do i=1,12
        read(13,'(A)',end=307) scr
      enddo
     
 306  read(13,*,end=307) ztmp,itmp,atmp,mfrac,amass,spin,cel
c      write(*,*) 'BBB: ',ztmp,itmp,atmp,mfrac,amass,spin,'x',cel,'x'
      sdaa(ztmp,itmp)=atmp
      sdfr(ztmp,itmp)=mfrac/100.d0 ! % to fraction
      sdma(ztmp,itmp)=amass*(m_p+m_n)/2.0d0
      sdsp(ztmp,itmp)=spin
      sdname(ztmp)=cel(1:2)
      goto 306

 307  close(13)
      if (prtlevel.ge.3) write(*,*) 'done.'

c...Mass numbers
c      sdaa(1,1)  = 1.0d0  ! H
c      sdaa(2,2)  = 4.0d0  ! He4
c      sdaa(2,1)  = 3.0d0  ! He3
c      sdaa(6,1)  = 12.0d0 ! C12
c      sdaa(6,2)  = 13.d0  ! C13
c      sdaa(7,1)  = 14.0d0 ! N14
c      sdaa(7,2)  = 15.0d0 ! N15
c      sdaa(8,1)  = 16.0d0 ! O16
c      sdaa(8,2)  = 17.0d0 ! O17
c      sdaa(8,3)  = 18.0d0 ! O18
c      sdaa(10,0) = 20.0d0 ! Ne
c      sdaa(11,0) = 23.0d0 ! Na
c      sdaa(12,0) = 24.0d0 ! Mg
c      sdaa(13,0) = 27.0d0 ! Al
c      sdaa(14,0) = 28.0d0 ! Si
c      sdaa(15,0) = 31.0d0 ! P
c      sdaa(16,0) = 32.0d0 ! S
c      sdaa(17,0) = 35.0d0 ! Cl
c      sdaa(18,0) = 36.0d0 ! Ar
c      sdaa(19,0) = 39.0d0 ! K
c      sdaa(20,0) = 40.0d0 ! Ca
c      sdaa(21,0) = 45.0d0 ! Sc
c      sdaa(22,0) = 48.0d0 ! Ti
c      sdaa(23,0) = 51.0d0 ! V
c      sdaa(24,0) = 52.0d0 ! Cr
c      sdaa(25,0) = 55.0d0 ! Mn
c      sdaa(26,0) = 56.0d0 ! Fe
c      sdaa(27,0) = 59.0d0 ! Co
c      sdaa(28,0) = 59.0d0 ! Ni

c...Masses
c      do m=0,isomax
c        do l=1,zmax
c          sdma(l,m)=sdaa(l,m)*(m_p+m_n)/2.0d0
c        enddo
c      enddo

c...Calculate isotope mix masses, i.e. sdaa(z,0) and sdma(z,0)
c...In principle, we have read in Atomic masses from table 1 in AGSS09, but
c...those are not very accurate, so here we recalculate them from the
c...correct solar abundances of the isotopes.
c...JE FIXME: Clean up the stuff above where table 1 is read
      do l=1,zmax
c         if (sdma(l,0).lt.1.d-5) then ! no isotope average exists
            atmp=0.d0
            amass=0.d0
c            write(*,*) 'CCC1: ',l,sdaa(l,0),sdma(l,0)
            do m=1,isomax
               atmp=atmp+sdaa(l,m)*sdfr(l,m)
               amass=amass+sdma(l,m)*sdfr(l,m)
            enddo
            sdaa(l,0)=atmp
            sdma(l,0)=amass
c            write(*,*) 'CCC2: ',l,atmp,amass
c         endif
      enddo
         

c...Now read in data from solar model data file
      if (prtlevel.ge.2) write(*,*) 'dssem_sunread: Opening file ',sunfile
      open(unit=13,file=sunfile,
     &  form='formatted',status='old')
      sdn=2  ! will add first line later

c...Read in scratch lines, NOTE: this might have to change with newer files
      do i=1,20
        read(13,'(A)',end=312) scr
      enddo

c...Calculate the total mass fraction (relative to H) of elements
c...heavier than O. This is just for normalization of the fractions
c...as a function of radius below
      totfrheavy=0.0d0
      do i=29,zmax
        totfrheavy=totfrheavy+sdma(i,1)*10**(sdabund(i)-12.0d0)
      enddo

c...Read in table and calculate remaining mass fractions from
c...average abundances sdabund
 311  read(13,*,end=312) sdm(sdn),sdr(sdn),tmp1,sdrho(sdn),tmp2,tmp3,
     &  sdmfr(1,1,sdn),sdmfr(2,2,sdn),sdmfr(2,1,sdn),sdmfr(6,1,sdn),
     &  sdmfr(6,2,sdn),sdmfr(7,1,sdn),sdmfr(7,2,sdn),sdmfr(8,1,sdn),
     &  sdmfr(8,2,sdn),sdmfr(8,3,sdn),sdmfr(10,0,sdn),sdmfr(11,0,sdn),
     &  sdmfr(12,0,sdn),sdmfr(13,0,sdn),sdmfr(14,0,sdn),sdmfr(15,0,sdn),
     &  sdmfr(16,0,sdn),sdmfr(17,0,sdn),sdmfr(18,0,sdn),sdmfr(19,0,sdn),
     &  sdmfr(20,0,sdn),sdmfr(21,0,sdn),sdmfr(22,0,sdn),sdmfr(23,0,sdn),
     &  sdmfr(24,0,sdn),sdmfr(25,0,sdn),sdmfr(26,0,sdn),sdmfr(27,0,sdn),
     &  sdmfr(28,0,sdn)
      totfr=0.d0
      do m=1,isomax
         do l=1,zmax
            totfr=totfr+sdmfr(l,m,sdn)
         enddo
      enddo

c...Add the heavy elements (>Ni), scaling to remaining heavy element fraction
c      do i=29,zmax
c        sdmfr(i,1,sdn)=(1.0d0-totfr)*sdma(i,1)*10**(sdabund(i)-12.0d0)
c     &    /totfrheavy      
c      enddo

c...Add heavy elements (>Ni), scaling to Fe
      do i=29,zmax
         sdmfr(i,0,sdn)=sdmfr(26,0,sdn)*
     &     10**(sdabund(i))/10**(sdabund(26))
      enddo

      sdn=sdn+1
      if (sdn.gt.sdmax-1) then  ! need one more entry for last line
        goto 312
      endif

      goto 311
 312  continue
      sdn=sdn-1
      close(13)

c...Now add the first line with r/r_sun=0
      sdm(1)=0.0d0
      sdr(1)=0.0d0
      sdrho(1)=sdrho(2)
      do m=0,isomax
         do l=1,zmax
           sdmfr(l,m,1)=sdmfr(l,m,2)
         enddo
      enddo

c...Now add the last line with r/r_sun=1
      sdn=sdn+1
      sdm(sdn)=1.0d0
      sdr(sdn)=1.0d0
      sdrho(sdn)=0.0d0
      do m=0,isomax
         do l=1,zmax
           sdmfr(l,m,sdn)=sdmfr(l,m,sdn-1)
         enddo
      enddo

c...Add isotope abundances for the stuff heavier than oxygen, i.e.
c...where we only have values for the isotope average (m=0) right now.
c...This will not affect much for spin-independent scattering,
c...but could for spin-dependent. Add it anyway, so that we have a complete
c...isotope abundance description of the Sun. We can always decide later
c...if we want to use this full information or do approximations.

      do l=1,zmax 
         if (sdmfr(l,1,int(sdn/2)).lt.1.d-30
     &      .and.sdmfr(l,0,int(sdn/2)).gt.0.d0) then ! no entry for isotopes
c            write(*,*) 'Adding isotope information for element ',l
            do m=1,isomax
               do i=1,sdn
                  sdmfr(l,m,i)=sdmfr(l,0,i)*sdfr(l,m)
               enddo
            enddo
         endif
       enddo


c...Electron density
c...Now read in data from solar model electron density file
      if (prtlevel.ge.2) write(*,*) 'dssem_sunread: Opening file ',sunnefile
      open(unit=13,file=sunnefile,
     &  form='formatted',status='old')
      sdnne=2  ! will add first line later

c...Read in scratch lines, NOTE: this might have to change with newer files
      do i=1,25
        read(13,'(A)',end=321) scr
      enddo

c...Read in table
 320  read(13,*,end=321) sdrne(sdnne),tmp1,sdne(sdnne)
      sdnne=sdnne+1
      if (sdnne.gt.sdmax) then  ! need one more entry for last line
        write(*,*) 'ERROR in dssem_sunread: array too small.'
        write(*,*) 'Increase sdmax.'
        goto 321
      endif
      goto 320
 321  continue
      sdnne=sdnne-1
      close(13)

c...Now add the first line with r/r_sun=0
      sdne(1)=sdne(2)
      sdrne(1)=0.0d0

      else
         write(*,*) 'DS ERROR in dssem_sunread: invalid sunfiletype: ',
     &     sunfiletype
         stop
      endif

c----------------------------------------------------------------------
c======= General setups for the loaded files =======
c----------------------------------------------------------------------

c...Sun parameters
      m_sun = 1.9889d30   ! solar mass in kg
      r_sun = 6.9598d8    ! solar radius in m

c...Selection of which elements to include is now moved to
c...dssenu_selectelements

c=======================================================================
c=======================================================================
c=======================================================================

c...Now calculate and tabulate the potential inside the Sun

      if (prtlevel.ge.3) write(*,*) 
     &  'dssem_sunread: tabulating potential inside the Sun...'
      do i=1,sdn
        sdphi(i)=dssem_sunpotint(sdr(i)*r_sun)
      enddo
      if (prtlevel.ge.3) write(*,*) '  ...done'

      if (prtlevel.ge.3) write(*,*) 
     &  'dssem_sunread: tabulating column density inside the Sun...'
      do i=1,sdn
        sdcdens(i,0)=dssem_suncdensint(sdr(i)*r_sun,'N')
        sdcdens(i,1)=dssem_suncdensint(sdr(i)*r_sun,'p')
        sdcdens(i,2)=dssem_suncdensint(sdr(i)*r_sun,'n')
      enddo
      if (prtlevel.ge.3) write(*,*) '  ...done'
      cd_sun(0)=sdcdens(sdn,0)  ! total (p+n) column density g/cm^2
      cd_sun(1)=sdcdens(sdn,1)  ! total proton column density g/cm^2
      cd_sun(2)=sdcdens(sdn,2)  ! total neutron column density g/cm^2

c      write(*,*) 'Ip = ',cd_sun(1)/r_sun/100.0d0
c      write(*,*) 'In = ',cd_sun(2)/r_sun/100.0d0

      return

      end


      

