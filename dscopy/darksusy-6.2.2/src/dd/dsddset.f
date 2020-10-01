      recursive subroutine dsddset(key,value)
c...
c...  type : commonly used
c...  desc : Set parameters for scattering cross sections
c...  c - character string specifying choice to be made
c...author: paolo gondolo 2000-07-07
c...modified by Gintaras Duda 2007-06-27 for new FF options
c...modified by Paolo Gondolo 2008-02-18
c...modified by Torsten Bringmann 2014-12-10: removed model-specific part
c...modified by Paolo Gondolo 2016-02-06: completely removed mssm part
c...modified by Paolo Gondolo 2016-11-19: added dsddsf option with values sisd, eft (make better later)
c...modified by Paolo Gondolo 2016-11-20: cleaned up the options
      implicit none
      include 'dsddcom.h'
      character*(*) key,value
      integer ierr

      ierr=1

      if (key.eq.'sf') then
        call dsddset('sf_m',value)
        call dsddset('sf_sigma',value)
        call dsddset('sf_delta',value)
        call dsddset('sf_phi',value)
        call dsddset('sf_deltasigma',value)
        call dsddset('sf_phim',value)
        ierr=0

      else if (key.eq.'sf_m') then
        if (value.eq.'default') then
          call dsddsetsf(ddsf(1),'best')
          ierr=0
        else if (
     &         value.eq.'q=0'
     &         .or.value.eq.'best'
     &         .or.value.eq.'FB'
     &         .or.value.eq.'SOG'
     &         .or.value.eq.'Fermi'
     &         .or.value.eq.'L-S'
     &         .or.value.eq.'haxton'
     &         .or.value.eq.'ds4.1'
     &         .or.value.eq.'gauss'
     &         .or.value.eq.'gould'
     &         ) then
          call dsddsetsf(ddsf(1),value)
          ierr=0
        endif

      else if (key.eq.'sf_sigma') then
        if (value.eq.'default') then
          call dsddsetsf(ddsf(2),'best')
          ierr=0
        else if (
     &         value.eq.'q=0'
     &         .or.value.eq.'best'
     &         .or.value.eq.'ISM'
     &         .or.value.eq.'haxton'
     &         .or.value.eq.'OddG'
     &         .or.value.eq.'SPSM'
     &         .or.value.eq.'SIMS'
     &         .or.value.eq.'ISMR'
     &         ) then
          call dsddsetsf(ddsf(2),value)
          ierr=0
        endif

      else if (key.eq.'sf_delta') then
        if (value.eq.'default') then
          call dsddsetsf(ddsf(3),'best')
          ierr=0
        else if (
     &         value.eq.'q=0'
     &         .or.value.eq.'best'
     &         .or.value.eq.'haxton'
     &         ) then
          call dsddsetsf(ddsf(3),value)
          ierr=0
        endif

      else if (key.eq.'sf_phi') then
        if (value.eq.'default') then
          call dsddsetsf(ddsf(4),'best')
          ierr=0
        else if (
     &         value.eq.'q=0'
     &         .or.value.eq.'best'
     &         .or.value.eq.'haxton'
     &         ) then
          call dsddsetsf(ddsf(4),value)
          ierr=0
        endif

      else if (key.eq.'sf_deltasigma') then
        if (value.eq.'default') then
          call dsddsetsf(ddsf(5),'best')
          ierr=0
        else if (
     &         value.eq.'q=0'
     &         .or.value.eq.'best'
     &         .or.value.eq.'haxton'
     &       ) then
          call dsddsetsf(ddsf(5),value)
          ierr=0
        endif

      else if (key.eq.'sf_phim') then
        if (value.eq.'default') then
          call dsddsetsf(ddsf(6),'best')
          ierr=0
        else if (
     &         value.eq.'q=0'
     &         .or.value.eq.'best'
     &         .or.value.eq.'haxton'
     &         ) then
          call dsddsetsf(ddsf(6),value)
          ierr=0
        endif

      else if (key.eq.'me') then
        call dsddset('me_s',value)
        call dsddset('me_v',value)
        call dsddset('me_vm',value)
        call dsddset('me_t',value)
        call dsddset('me_a',value)
        call dsddset('me_am',value)
        call dsddset('me_p',value)
        call dsddset('me_t2e',value)
        call dsddset('me_t2o',value)
        ierr=0

      else if (key.eq.'me_s') then
        if (value.eq.'default') then
          ddme(1) = 'gls91'
          ierr=0
        else if (
     &         value.eq.'gls91'
     &         ) then
          ddme(1) = value
          ierr=0
        endif

      else if (key.eq.'me_v') then
        if (value.eq.'default') then
          ddme(2) = ''
          ierr=0
        else if (
     &         .false.
     &         ) then
          ddme(2) = value
          ierr=0
        endif

      else if (key.eq.'me_vm') then
        if (value.eq.'default') then
          ddme(3) = ''
          ierr=0
        else if (
     &         .false.
     &         ) then
          ddme(3) = value
          ierr=0
        endif

      else if (key.eq.'me_t') then
        if (value.eq.'default') then
          ddme(4) = ''
          ierr=0
        else if (
     &         .false.
     &         ) then
          ddme(4) = value
          ierr=0
        endif

      else if (key.eq.'me_a') then
        if (value.eq.'default') then
          ddme(5) = 'smc'
          ierr=0
        else if (
     &         value.eq.'smc'
     &         .or.value.eq.'emc'
     &         ) then
          ddme(5) = value
          ierr=0
        endif

      else if (key.eq.'me_am') then
        if (value.eq.'default') then
          ddme(6) = ''
          ierr=0
        else if (
     &         .false.
     &         ) then
          ddme(6) = value
          ierr=0
        endif

      else if (key.eq.'me_p') then
        if (value.eq.'default') then
          ddme(7) = ''
          ierr=0
        else if (
     &         .false.
     &         ) then
          ddme(7) = value
          ierr=0
        endif

      else if (key.eq.'me_t2e') then
        if (value.eq.'default') then
          ddme(8) = ''
          ierr=0
        else if (
     &         .false.
     &         ) then
          ddme(8) = value
          ierr=0
        endif

      else if (key.eq.'me_t2o') then
        if (value.eq.'default') then
          ddme(9) = ''
          ierr=0
        else if (
     &         .false.
     &         ) then
          ddme(9) = value
          ierr=0
        endif

      endif

      if (ierr.ne.0) then
        write (*,*) 'dsddset: unrecognized option ''',key,''' ''',value,''''
        return
      endif

! obsolescent
! the rest of this subrotoutine sets the values of the matrix elements
! and will eventually be moved to its own subroutine
! it will also be switched to include all matrix elements for all currents

      if (ddme(1).eq.'gls91') then
c...gaisser,leutwyler,sainio,plb253(91)252
c...drees,nojiri,prd48(93)3483
c...bergstrom,gondolo,astrop.ph.
        ftp(7)=0.023
        ftp(8)=0.034
        ftp(9)=0.0595
        ftp(10)=0.14
        ftp(11)=0.0595
        ftp(12)=0.0595
        ftn(7)=0.019
        ftn(8)=0.041
        ftn(9)=0.0592
        ftn(10)=0.14
        ftn(11)=0.0592
        ftn(12)=0.0592
      endif
            
      if (ddme(5).eq.'smc') then
c...bergstrom,gondolo,astrop.ph.
c...d.adams et al cern-ppe/94-57
        delu=0.74
        deld=-0.40
        dels=-0.12
      else if (ddme(5).eq.'emc') then
c...drees,nojiri,prd48(93)3483
        delu=0.77
        deld=-0.49
        dels=-0.15
      endif

c      do i=1,ddnsf
c        write (*,'(a,i1,a,a,a)') 'PG> ddsf(',i,')=/',ddsf(i),'/'
c      enddo
c      do i=1,ddnme
c        write (*,'(a,i1,a,a,a)') 'PG> ddme(',i,')=/',ddme(i),'/'
c      enddo
      
      return
      end


      subroutine dsddsetsf(sf,value)
      character*(*) sf,value
      if (value.eq.'q=0') then
        sf(1:1) = '0'
      else
        sf = sf(1:1)//value
      endif
      end
