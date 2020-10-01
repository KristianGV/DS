      subroutine dssenu_selectelements

***********************************************************************
*** After having read a solar model file, we here select which
*** elements to include in capture rate calculation. This is determined
*** by the common block variable sesunacc.
*** Author: Joakim Edsjo
*** Date: 2015-06-15 (based on dssem_sunread code)
***********************************************************************

      implicit none
      include 'dssem_sun.h'
      include 'dssecom.h'
      include 'dsmpconst.h'
      include 'dsio.h'

      logical inciso
      integer l,m,m_1,m_2

      integer sesunaccsave
      data sesunaccsave/-1/
      save sesunaccsave

      if (sesunacc.eq.sesunaccsave) return ! no need to recalculate tables

      sesunaccsave=sesunacc

c----------------------------------------------------------------------
c======= General setups for the loaded files =======
c----------------------------------------------------------------------

c...Define list of elements and isotopes to include
c----Spin-independent
      sdelsi=0

      do l=1,zmax
c...    determine if isotopes should be included
        inciso=.false.
        if (sdmfr(l,1,int(sdn/2)).gt.0.d0) inciso=.true. ! iso information=yes
        if (inciso) then
           m_1=1
           m_2=isomax
        else
           m_1=0
           m_2=0
        endif
c...The above is the accurate approach (sesunacc=1). Now modify for less
c...accurate options
c...For sesunacc>=2, include isotopic averages above Oxygen
        if (sesunacc.ge.2.and.l.ge.10) then 
           m_1=0
           m_2=0
        endif

c        write(*,*) 'AAA: ',l,sdmfr(l,1,int(sdn/2)),sdmfr(l,0,int(sdn/2))
c     &    ,m_1,m_2

        do m=m_1,m_2
           if (sdmfr(l,m,int(sdn/2)).gt.0.d0) then
c              write(*,*) 'AA:  adding isotope ',m
c...For fast calculation (sesunacc=2,3), only include elements with 
c...appreciable abundances
c...Old selection
c              if (sesunacc.eq.2.and.
c     &          (10**(sdabund(l)-12)*sdaa(l,m)**4).lt.1.d-2) goto 980
c              if (sesunacc.eq.3.and.
c     &          (10**(sdabund(l)-12)*sdaa(l,m)**4).lt.1.d-1) goto 980
c...New better selection
              if (sesunacc.eq.2.and.
     &          (sdmfr(l,m,int(sdn/2))*sdaa(l,m)**2).lt.1.d-2) goto 980
              if (sesunacc.eq.3.and.
     &          (sdmfr(l,m,int(sdn/2))*sdaa(l,m)**2).lt.0.3d0) goto 980
              sdelsi=sdelsi+1
              sdzsi(sdelsi)=l
              sdisosi(sdelsi)=m
              if (prtlevel.ge.2) 
     &          write(*,*) 'Adding SI Z=',l,' i=',m,' A=',sdaa(l,m)
c              write(*,*) 'AA: Adding (Z,I) = ',l,m,' ',
c     &           sdname(l),'-',int(sdaa(l,m)+0.5),
c     &           10**(sdabund(l)-12)*sdaa(l,m)**4
              if (sdelsi.ge.sdelmax) then
                 write(*,*) 'DS ERROR in dssnu_selectelements: ',
     &             'sdelmax too small for SI'
                 stop
              endif
 980          continue
            endif
        enddo
      enddo
      write(*,'(A,I3,A)') 'dssenu_selectelements: ',sdelsi,
     &   ' elements will be included in the Sun SI capture calculation.'

c---Spin-dependent elements
      sdelsd=0
      do l=1,zmax
         do m=1,isomax
            if (sdsp(l,m).ne.0.d0) then
               if (sdmfr(l,m,int(sdn/2)).gt.0.d0) then
c...For fast calculation (sesunacc=2,3), only include elements with 
c...appreciable abundances
c...Old code                  
c                  if (sesunacc.eq.2.and.
c     &            (10**(sdabund(l)-12)*sdaa(l,m)**2).lt.1.d-4) goto 990
c                  if (sesunacc.eq.3.and.
c     &            (10**(sdabund(l)-12)*sdaa(l,m)**2).lt.1.d-3) goto 990
c...New better selection
              if (sesunacc.eq.2.and.
     &          (sdmfr(l,m,int(sdn/2)).lt.1.d-6)) goto 990
              if (sesunacc.eq.3.and.
     &          (sdmfr(l,m,int(sdn/2)).lt.5.d-5)) goto 990
                  sdelsd=sdelsd+1
                  sdzsd(sdelsd)=l
                  sdisosd(sdelsd)=m
                  if (prtlevel.ge.2) 
     &              write(*,*) 'Adding SD Z=',l,' i=',m,' A=',sdaa(l,m)
               endif
               if (sdelsd.ge.sdelmax) then
                  write(*,*) 'DS ERROR in dssenu_selectelements: ',
     &              'sdelmax too small for SD'
                  stop
               endif
 990           continue
            endif
         enddo
      enddo

      write(*,'(A,I3,A)') 'dssenu_selectelements: ',sdelsd,
     &     ' elements will be included in the Sun SD capture calculation.'

      return

      end


      

