      subroutine dsflavour(brbsg,delta0,brbmumu,brbmumuuntag,
     &   brbtaunu,brbdtaunu,amu,rmu23,iflag) 
c_______________________________________________________________________
c  Routine that calls SuperIso and calculates various flavour
c  observables. More are easily added here.
c  The interface goes via a SLHA file. For SuperIso to be able to read it
c  we need to set it to no-CP-violation and MFV.
c  output:
c    brbsg:   Br(b -> s gamma)
c    brbmumu: Br(B_s -> mu+ mu-)
c    iflag: flag that if non-zero indicates a problem
c  author: Joakim Edsjo, edsjo@fysik.su.se
c  date: 2013-05-01
c=======================================================================
      implicit none
      
      include 'dsidtag.h'
      include 'dsio.h'
      real*8 brbsg,brbmumu,brbmumuuntag,delta0
      real*8 brbtaunu,brbdtaunu,amu,rmu23
      character*80 slhafile
      character*8 cpid
      integer iflag,test,pid

c---------- functions
      integer test_slha,getpid
      real*8 bsgamma_calculator
      real*8 delta0_calculator
      real*8 bsmumu_calculator
      real*8 bsmumu_untag_calculator
      real*8 btaunu_calculator
      real*8 bdtaunu_calculator
      real*8 muon_gm2_calculator
      real*8 rmu23_calculator

c--------------------      
  
      pid=getpid() ! process id, to generate unique file name
      write(cpid,'(I8)') pid
c...Generate slha filename and save model to file
      slhafile=idtag//'-'//cpid//'-tmp.slha'
      call dscharadd(slhafile,char(0))

      if (prtlevel.gt.2) 
     &   write(*,*) 'Will write SLHA file to file: ',slhafile

      call dsslhawrite(slhafile,2)

c...Now ask SuperIso to calculate flavour observables. This is more or less
c... a copy of what is done in SuperIso's slha.c

      test=test_slha(slhafile)
c      write(*,*) 'test_slha: ',test ! JE TMP
      if (test.gt.0) then
         brbsg=bsgamma_calculator(slhafile)
         delta0=delta0_calculator(slhafile)
         brbmumu=bsmumu_calculator(slhafile)
         brbmumuuntag=bsmumu_untag_calculator(slhafile)
         brbtaunu=btaunu_calculator(slhafile)
         brbdtaunu=bdtaunu_calculator(slhafile)
         amu=muon_gm2_calculator(slhafile)
         rmu23=rmu23_calculator(slhafile)

c         brbdmumu=bdmumu_calculator(slhafile)
      else
         brbsg=0.d0
         delta0=0.d0
         brbmumu=0.d0
         brbmumuuntag=0.d0
         brbtaunu=0.d0
         brbdtaunu=0.d0
         amu=0.d0
         rmu23=0.d0
         iflag=test
      endif

c      write(*,*) ' Br(b->s gamma):   ',brbsg ! JE TMP
c      write(*,*) ' Delta-0:          ',delta0 ! JE TMP
c      write(*,*) ' Br(B_s->mu+ mu-): ',brbmumu ! JE TMP

c...Now delete the file
      open(unit=77,file=slhafile,status="old")
      close(77,status="delete")

      return

      end
