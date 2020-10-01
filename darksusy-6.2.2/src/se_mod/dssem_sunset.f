      subroutine dssem_sunset(c)
c...determine which solar model to use
c...  c - character string specifying choice to be made
c...author: joakim edsjo, 2009-11-10
      implicit none

      include 'dssem_sun.h'

      integer fl,l,m
      character*(*) c

c...Reset read variable as we need to reread solar table evertime the solar
c...model is changed
      sdread=.false.

      if (c.eq.'bs05op') then
c...     This is from the Standard solar model, BS05(OP)
c...     The mass fractions for heavier elements are from N. Grevesse and
c...     A.J. Sauval, Space Science Reviews 85 (1998) 161 normalized such that
c...     their total mass fractions matches that of the heavier elements in 
c...     the BS05op model.
         sunid='bs05op'
         sunfiletype=1
         call dsdatafile(sunfile,'bs05op.dat')
c...     The electron density is from the Standard solar model, BS05(OP)
         call dsdatafile(sunnefile,'nele_bs05op.dat')

      elseif (c.eq.'bs05agsop') then
c...     This is form the alternative Standard solar model, BS05(AGS,OP)
c...     with new heavy element measurements (fits worse with helioseismology
c...     so we don't use it as a defualt).
c...     The mass fractions for heavier elements are from N. Grevesse and
c...     A.J. Sauval, Space Science Reviews 85 (1998) 161 normalized such that
c...     their total mass fractions matches that of the heavier elements in 
c...     the BS05agsop model.
         sunid='bs05agsop'
         sunfiletype=1
         call dsdatafile(sunfile,'bs05_agsop.dat')
c...     The electron density is from the Standard solar model, BS05(OP)
         call dsdatafile(sunnefile,'nele_bs05op.dat')

      elseif (c.eq.'AGSS09'.or.c.eq.'agss09'.or.c.eq.'default') then
         sunid='agss09'
         sunfiletype=2
         call dsdatafile(sunfile,
     &      'Serenelli-model_agss09.dat')
         call dsdatafile(sunnefile,
     &      'Serenelli-flux_distrib_agss09.dat')
         call dsdatafile(sunabundfile,
     &      'AGSS09-table1.txt')
         call dsdatafile(sunisofile,
     &      'solar-isotopes.dat')
         absrc='ci' ! CI meteorites

      elseif (c.eq.'AGSS09ph'.or.c.eq.'agss09ph') then
         sunid='agss09ph'
         sunfiletype=2
         call dsdatafile(sunfile,
     &      'Serenelli-model_agss09ph.dat')
         call dsdatafile(sunnefile,
     &      'Serenelli-flux_distrib_agss09ph.dat')
         call dsdatafile(sunabundfile,
     &      'AGSS09-table1.txt')
         call dsdatafile(sunisofile,
     &      'solar-isotopes.dat')
         absrc='ph' ! photoshperic

      elseif (c.eq.'AGS05'.or.c.eq.'ags05') then
         sunid='ags05'
         sunfiletype=2
         call dsdatafile(sunfile,
     &      'Serenelli-model_ags05.dat')
         call dsdatafile(sunnefile,
     &      'Serenelli-flux_distrib_ags05.dat')
         call dsdatafile(sunabundfile,
     &      'AGSS09-table1.txt')
         call dsdatafile(sunisofile,
     &      'solar-isotopes.dat')
         absrc='ph' ! photoshperic

      elseif (c.eq.'GS98'.or.c.eq.'gs98') then
         sunid='gs98'
         sunfiletype=2
         call dsdatafile(sunfile,
     &      'Serenelli-model_gs98.dat')
         call dsdatafile(sunnefile,
     &      'Serenelli-flux_distrib_gs98.dat')
         call dsdatafile(sunabundfile,
     &      'AGSS09-table1.txt')
         call dsdatafile(sunisofile,
     &      'solar-isotopes.dat')
         absrc='ph' ! photoshperic
      
      else
         write(*,*) 'Invalid solar model: ',c
         write(*,*) 'Check dssem_sunset for valid options.'
         stop
      endif

c...delete possible spaces in file name
      fl=200
      do l=1,fl
 40     if (sunfile(l:l).eq.' ') then
          fl=fl-1
          do m=l,fl
            sunfile(m:m)=sunfile(m+1:m+1)
          enddo
          if (fl.eq.l) goto 50
          goto 40
        endif
      enddo
 50   continue

c...delete possible spaces in file name
      fl=200
      do l=1,fl
 60     if (sunnefile(l:l).eq.' ') then
          fl=fl-1
          do m=l,fl
            sunnefile(m:m)=sunnefile(m+1:m+1)
          enddo
          if (fl.eq.l) goto 70
          goto 60
        endif
      enddo
 70   continue

      if (sunfiletype.eq.2) then

c...delete possible spaces in file name
      fl=200
      do l=1,fl
 80     if (sunabundfile(l:l).eq.' ') then
          fl=fl-1
          do m=l,fl
            sunabundfile(m:m)=sunabundfile(m+1:m+1)
          enddo
          if (fl.eq.l) goto 81
          goto 80
        endif
      enddo
 81   continue

c...delete possible spaces in file name
      fl=200
      do l=1,fl
 90      if (sunisofile(l:l).eq.' ') then
          fl=fl-1
          do m=l,fl
            sunisofile(m:m)=sunisofile(m+1:m+1)
          enddo
          if (fl.eq.l) goto 91
          goto 90
        endif
      enddo
 91   continue

      endif


      return
      end

         

         

         
