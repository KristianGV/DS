      subroutine dsreadnuclides
      implicit none
      include 'dsnuclides.h'
      include 'dsmpconst.h'
      include 'dsio.h'
      character*200 filename,scratch*2
      character*300 msg
      integer itmp,i
      real*8 tmp,etmp
      call dsdatafile(filename,'nuclides.dat')
      if (prtlevel.ge.2) then
        write (msg,'(a,a)') 'dsreadnuclides: reading ',filename
        call dswrite(0,0,msg)
      endif  
      open (unit=13,file=filename,status='unknown',
     &     form='formatted')
      read(13,'(a)',err=2000) scratch
      do i=1,nnucld
         read(13,*,err=2000) nucldsym(i),nucldz(i),nuclda(i),nucldn(i),
     &        tmp,itmp,nucldnatab(i),nucldj(i),nucldmu(i),etmp
         nucldam(i)=tmp*atomicmassunit
         nucldm(i)=(nucldz(i)*m_p_amu+nucldn(i)*m_n_amu)*atomicmassunit
         if (etmp.ne.9999.d0) then
            nucldm(i)=nucldm(i)-etmp*nuclda(i)*1.d-3
         endif
         if (itmp.eq.0) nucldstable(i)=.true.
         if (itmp.eq.1) nucldstable(i)=.false.
      enddo
      close(13)
      inucld=1
ccc begin dbg
ccc      do i=1,nnucld
ccc         write(*,*) nucldsym(i),nucldz(i),nuclda(i),nucldn(i),nucldm(i),
ccc     &        nucldstable(i),nucldnatab(i),nucldj(i),nucldmu(i)
ccc      enddo
ccc end dbg
      return
 2000 continue
      close(13)
      write (*,*) 'dsreadnuclides: error while reading ',filename
      write (*,*) 'DarkSUSY will stop'
      stop
      end
