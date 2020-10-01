**********************************************************************
*** program generic_wimp_oh2 loops over WIMP mass and finds the
*** annihilation cross sections (for the chosen annihilation channel)
*** that give a good relic density (within the Planck 2015, 3 sigma
*** limits)
*** Author: Joakim Edsjo, edsjo@fysik.su.se
*** Date: July 1, 2016
**********************************************************************

      program generic_wimp_oh2
      implicit none
	
      integer i,n,j
      real*8 svlow(4),svmid(4),svhigh(4),mwimp
      real*8 oh2low,oh2mid,oh2high
      real*8 dsmass
      real*8 mmin,mmax
      real*8 f,fth
      logical nearthreshold
      logical selfconj
      integer pdgann(4) ! use three pdg codes
      real*8 findsv
      parameter(oh2low=0.1165d0,
     &     oh2mid=0.1193d0,
     &     oh2high=0.1221d0)    ! Planck, arXiv:1502.01589,
                                ! tab 4, TT,TE,EE+lowP+lensing , 0.1193±0.0014

      real*8 pp,dsanwx,dsrdthav,tmp,xx,dsrdwintp ! JE TMP
      external dsanwx,dsrdwintp ! JE TMP
      character*80 filename

c...Initialize DarkSUSY
      call dsinit
      
c...Setup
      mmin=0.1d0                ! starting mass
c      mmin=82.1382d0 ! JE TMP
c      mmin=88.9091d0
      mmax=1.d4                 ! ending mass
      f=1.1d0     ! mass increase between steps
      fth=1.02d0  ! mass increase between steps when close to threshold
      
c...Filename for output file
      selfconj=.true.           ! self-conjugated WIMP
      pdgann(1)=12 ! nu_e nu_e-bar
c      pdgann(2)=5 ! b b-bar ! not a suitable choice as we don't have free
                   ! quarks at this temperature
      pdgann(2)=15 ! tau- tau+
      pdgann(3)=6 ! t t-bar
      pdgann(4)=24 ! W- W+
      filename='generic_wimp_oh2-planck-sigmav.dat' 

      open(unit=71,file=filename,status='unknown',form='formatted')
      write(71,'(A)')
     &     '# Thermal sigma v for some different final states'
      write(71,'(A)')
     &     '# The final states are:'
      do i=1,4
         write(71,'(A,I1,A,I2)') '# Case ',i,': PDG = ',pdgann(i)
      enddo
      
      write(71,'(10(A))') '#  i     Mass [GeV]      ',
     &   'sigvmid_1     sigvperr_1     sigvmerr_1      ',
     &   'sigvmid_2     sigvperr_2     sigvmerr_2      ',
     &   'sigvmid_3     sigvperr_3     sigvmerr_3      ',
     &   'sigvmid_4     sigvperr_4     sigvmerr_4,',
     &  ' sigv in units of cm^3/s^1'

c...Start scanning over masses, we do this by setting mwimp to the starting
c...mass and then increase it by a factor f for the next point. When
c...close to a threshold, the increase is instead fth to sampel this region
c...better. The scaning ends when mwimp has reached mmax

      mwimp=mmin
      i=0
      
 10   i=i+1                     ! loop starts here
      do j=1,4
c         write(*,*) j  ! JE TMP
         svmid(j)=findsv(mwimp,oh2mid,pdgann(j),selfconj,0)         
         svlow(j)=findsv(mwimp,oh2high,pdgann(j),selfconj,0)
         svhigh(j)=findsv(mwimp,oh2low,pdgann(j),selfconj,0)
      enddo

c      if (i.eq.2) then ! JE TMP this if statement
c         do j=0,1000
c            pp=0.d0+dble(j)/dble(1000)*500.d0
c            tmp=dsanwx(pp)
c            write(47,*) pp,tmp
c         enddo
c         stop
c      endif

c      if (i.eq.12) then          ! JE TMP this if statement
c         do j=0,1000
c            xx=dble(j)/dble(1000.d0)*30.d0+2.d0
c            tmp=dsrdthav(xx,dsrdwintp)
c            write(48,*) xx,tmp
c         enddo
c         stop
c      endif
      
            
            
         
      write(71,100) i,mwimp,(svmid(j),svhigh(j)-svmid(j),
     &   svlow(j)-svmid(j),j=1,4)
      write(*,100)  i,mwimp,svmid(1),svmid(2),svmid(3),svmid(4)


c...check if we are near a threshold
      nearthreshold=.false.
      do j=1,4
c         if (mwimp.ge.0.80d0*dsmass(pdgann(j))
c     &         .and.mwimp.le.1.1d0*dsmass(pdgann(j))) then
c            nearthreshold=.true.
c         endif
         if (mwimp.ge.0.65d0*dsmass(pdgann(j))
     &         .and.mwimp.le.1.2d0*dsmass(pdgann(j))) then
            nearthreshold=.true.
         endif
      enddo

c...increase mwimp for next step      
      if (nearthreshold) then
         mwimp=mwimp*fth
      else
         mwimp=mwimp*f
      endif

c...Special treatment of last point to make sure we finish at mmax      
      if (mwimp.lt.mmax.and.mwimp*f.gt.mmax) mwimp=mmax

      if (mwimp.le.mmax) goto 10 ! not yet reached the final mass

c...Scanning over masses done      
         
      close(71)

 100  format (I4,100(1x,E14.6))
      write(*,*) 'Done.'
      end


***************
*** Subroutines
***************

      real*8 function findsv(mwimp,oh2goal,pdgann,selfconj,printlevel)
c...Subroutine to find the signa v that gives relic density oh2goal
c...It makes a binary search in sigmav until oh2 with the required
c...relative error has been found
c...This binary search assumes that omega h^2 is monotonically decreasing
c...with sigma v. As there can be some numerical undertainties in the
c...relic density estimate that are larger than the uncertainty we require
c...we stop after nmax (see below) omega evaluations even if required
c...accuracy has not been found.
c...Inputs
c...   mwimp    = WIMP mass in GeV (real*8)
c...   oh2goal  = what relic density Omega h^2 to search for (real*8)
c...   pdgann   = PDG code of final state particle (integer)
c...   selfconj = if WIMP is self-conjugate or not (logical)
c...   printlevl: 0 = no warnings are printed
c...              1 = warnings are printed       

      implicit none
      real*8 mwimp,oh2goal,oh2,oh2b,oh2a,xf,si,sv
      integer pdgann
      integer ierr,iwar,nfc,printlevel
      logical selfconj
      real*8 dsrdomega
      integer npoints,nmax
      parameter(nmax=100)

      real*8 svupper,svlower,svb,sva,oh2relerr
      parameter(oh2relerr=1.d-4, ! relative error on oh2 when exiting
     &     svupper=1e-18,       ! upper sigma v searched for
     &     svlower=1e-29)       ! lower sigma v searched for

      data svb/svupper/
      data sva/svlower/
      save sva,svb

      si=1.d-5 ! scattering cross section in pb, irrelelevant for us here

      if (svb.ne.svupper) then
         svb=svb*2.d0
         sva=sva/2.d0
      else         
         svb=svupper
         sva=svlower
      endif

      npoints=0
         
 5    call dsgivemodel_generic_wimp(mwimp,selfconj,
     &     svb,pdgann,SI)
      oh2b=dsrdomega(0,0,xf,ierr,iwar,nfc)

      call dsgivemodel_generic_wimp(mwimp,selfconj,
     &     sva,pdgann,SI)
      oh2a=dsrdomega(0,0,xf,ierr,iwar,nfc)
      npoints=npoints+2
      
      if (((oh2b.gt.oh2goal).or.(oh2a.lt.oh2goal))) then
         if (npoints.eq.2) then ! first guess not so good
            svb=svupper
            sva=svlower
            goto 5
         endif
         if (printlevel.gt.0) then
            write(*,*)
     &         'Warning in findsv: oh2 out of bounds, returning 0'
            write(*,*) 'for svb = ',svb,' we get oh2 = ',oh2b
            write(*,*) 'for sva = ',sva,' we get oh2 = ',oh2a
         endif
         findsv=1.d0 ! returns one to make plotting easier
         return
      endif


 10   sv=10**((log10(sva)+log10(svb))/2.d0) ! binary searh in log sigma v
      call dsgivemodel_generic_wimp(mwimp,selfconj,
     &     sv,pdgann,SI)
      oh2=dsrdomega(0,0,xf,ierr,iwar,nfc)
      npoints=npoints+1

      if ((abs((oh2-oh2goal)/oh2goal).lt.oh2relerr).or.npoints.ge.nmax) then
         findsv=sv
         return
      endif

c...Divide interval in half, decide which half      
      if (oh2.gt.oh2goal) then  ! too low sigma v
         sva=sv
         oh2a=oh2
      else
         svb=sv
         oh2b=oh2
      endif
      goto 10

      return
      end
      
         
      
      
