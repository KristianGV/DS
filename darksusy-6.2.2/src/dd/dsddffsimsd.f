      subroutine dsddffsimsd(q,a,z,l2jjpp,l2jjnn,l2jjpn,ierr)
c_______________________________________________________________________
c  Spin-dependent structure functions for direct detection.
c  Simplified estimates where no more accurate model descriptions exist,
c  like for isotopes with two odd neucleons. This routine is based on the
c  the routines dsddoddgff, which radii from Ellis and Flores where
c  applicable for similiar isotopes, but spin part estimated rather crudely.
c  With or without gaussian form factor.
c  input:
c    q : real*8  : momentum transfer in GeV ( q=sqrt(M*E/2/mu^2) )
c    a : integer : mass number
c    z : integer : atomic number
c  output:
c    l2jjnn, l2jjpp, l2jjpn : the factors \lambda^2 J (J+1) made functions
c       of the momentum transfer q, for pp, nn, and pn; they are
c       related to Spp, Snn, and Spn through
c           S = \lambda^2 J (J+1) (2J+1)/\pi
c       In turn, Spp, Snn, and Spn are related to the spin-dependent 
c       structure functions S_{00}(q), S_{01}(q), S_{11}(q) defined in
c       Engel (PLB264,114,1991) via
c           Spp = S00+S11+S01
c           Snn = S00+S11-S01
c           Spn = 2*(S00-S11)
c  author of original routine: paolo gondolo (paolo@physics.utah.edu) 2008
c  author of this modified version: Joakim Edsjo (edsjo@fysik.su.se) 2010
c=======================================================================
      implicit none
      include 'dsnuclides.h'
      include 'dsio.h'
      include 'dsmpconst.h'
      real*8 q
      integer a,z,ierr
      integer n,dsnucldindx,i
      real*8 l2jjpp,l2jjnn,l2jjpn,l2jjef
      real*8 y,j,ff,r,rsprc,mun,mup,gn,gp,mu,lam
      real*8 avsn,avsp,lamp,lamn
c     j = total nuclear spin
      
      ierr=0
      l2jjpp=0.d0
      l2jjnn=0.d0
      l2jjpn=0.d0

      n=a-z

      if (mod(n,2).eq.0.and.mod(z,2).eq.0) return ! no spin

      i=dsnucldindx(a,z)
      if (i.eq.0) goto 1000

      j=nucldj(i)
      if (j.eq.0.d0) return
      if (j.eq.9999.d0) goto 1000

      
      mu=nucldmu(i)
      mun=nucldmu(1)
      mup=nucldmu(2)
      gn=2.d0*mun               ! g-factor neutron
      gp=2.d0*mup               ! g-factor proton

      
      
      if (mod(n,2).eq.1.and.mod(z,2).eq.1) then ! odd odd
        avsn=j/2.d0             ! very approximate
        avsp=j/2.d0             ! very approximate
        if (a.eq.2.and.z.eq.1) then ! H2, EF NPB 307 (1988) 883
c...The following assumes <S_p> = <S_n> 
          l2jjef = 0.445d0
          avsp=sqrt(l2jjef/(4.d0*(j+1.d0)/j))
          avsn=avsp
        elseif (a.eq.6.and.z.eq.3) then ! B10, EF NPB 307 (1988) 883
c...The following assumes <S_p> = <S_n> 
          l2jjef = 0.36d0
          avsp=sqrt(l2jjef/((j+1.d0)/j))
          avsn=avsp
        elseif (a.eq.10.and.z.eq.5) then ! B10, EF NPB 307 (1988) 883
c...The following assumes <S_p> = <S_n> 
          l2jjef = 1.d0/3.0d0
          avsp=sqrt(l2jjef/((j+1.d0)/j))
          avsn=avsp
        elseif (a.eq.14.and.z.eq.7) then ! N-14, EF NPB 307 (1988) 883
c...The following assumes <S_p> = <S_n> 
          l2jjef = 0.03d0       ! EF88
          avsp=sqrt(l2jjef/((j+1.d0)/j))
          avsn=avsp
        endif

        lamn=avsn/j
        lamp=avsp/j
        l2jjpp=lamp*lamp*j*(j+1.d0)
        l2jjnn=lamn*lamn*j*(j+1.d0)
        l2jjpn=2.d0*lamp*lamn*j*(j+1.d0)


      elseif (mod(n,2).eq.1.and.mod(z,2).eq.0) then ! odd neutron
        if (mu.eq.9999.d0) goto 1000
        lam=mu/gn/j
        l2jjpp=0.d0
        l2jjnn=lam*lam*j*(j+1.d0)
        l2jjpn=0.d0

      else                      ! odd proton
        if (mu.eq.9999.d0) goto 1000
        lam=(mu-j)/(gp-1.d0)/j
        l2jjpp=lam*lam*j*(j+1.d0)
        l2jjnn=0.d0
        l2jjpn=0.d0

      endif

      if (q.gt.0.d0) then
! table from Ellis and Flores PLB293,259(1991)
        if (a.eq.7.and.z.eq.3) then ! Li7
          rsprc=1.09d0
        else if (a.eq.11.and.z.eq.4) then ! B11
          rsprc=1.05d0
        else if (a.eq.14.and.z.eq.7) then ! N14 ! approximation
          rsprc=1.04d0
        else if (a.eq.15.and.z.eq.7) then ! N15
          rsprc=1.04d0
        else if (a.eq.19.and.z.eq.9) then ! F19
          rsprc=0.960d0
        else if (a.eq.23.and.z.eq.11) then ! Na23
          rsprc=1.08d0
        else if (a.eq.27.and.z.eq.13) then ! Al27
          rsprc=1.07d0
        else if (a.eq.31.and.z.eq.15) then ! P31
          rsprc=0.913d0
        else if (a.eq.35.and.z.eq.17) then ! Cl35
          rsprc=1.07d0
        else if (a.eq.51.and.z.eq.23) then ! V51
          rsprc=1.10d0
        else if (a.eq.55.and.z.eq.25) then ! Mn55
          rsprc=1.09d0
        else if (a.eq.59.and.z.eq.27) then ! Co59
          rsprc=1.08d0
        else if (a.eq.69.and.z.eq.31) then ! Ga69
          rsprc=0.918d0
        else if (a.eq.71.and.z.eq.31) then ! Ga71
          rsprc=0.918d0
        else if (a.eq.75.and.z.eq.33) then ! As75
          rsprc=0.928d0
        else if (a.eq.79.and.z.eq.35) then ! Br79
          rsprc=0.900d0
        else if (a.eq.81.and.z.eq.35) then ! Br81
          rsprc=0.900d0
        else if (a.eq.93.and.z.eq.41) then ! Nb93
          rsprc=1.11d0
        else if (a.eq.107.and.z.eq.47) then ! Ag107
          rsprc=0.894d0
        else if (a.eq.109.and.z.eq.47) then ! Ag109
          rsprc=0.894d0
        else if (a.eq.121.and.z.eq.51) then ! Sb121
          rsprc=0.922d0
        else if (a.eq.123.and.z.eq.51) then ! Sb123
          rsprc=1.09d0
        else if (a.eq.127.and.z.eq.53) then ! I127
          rsprc=0.928d0
        else if (a.eq.133.and.z.eq.55) then ! Cs133
          rsprc=1.13d0
        else if (a.eq.139.and.z.eq.57) then ! La139
          rsprc=1.12d0
        else if (a.eq.191.and.z.eq.77) then ! Ir191
          rsprc=0.903d0
        else if (a.eq.193.and.z.eq.77) then ! Ir193
          rsprc=0.903d0
        else if (a.eq.203.and.z.eq.81) then ! Tl203
          rsprc=0.848d0
        else if (a.eq.205.and.z.eq.81) then ! Tl205
          rsprc=0.848d0
        else if (a.eq.209.and.z.eq.83) then ! Bi209
          rsprc=1.10d0
        else if (a.eq.3.and.z.eq.2) then ! He3
          rsprc=1.00d0
        else if (a.eq.9.and.z.eq.4) then ! Be9
          rsprc=1.07d0
        else if (a.eq.17.and.z.eq.8) then ! O17
          rsprc=1.12d0
        else if (a.eq.29.and.z.eq.14) then ! Si29
          rsprc=0.908d0
        else if (a.eq.47.and.z.eq.22) then ! Ti47
          rsprc=1.10d0
        else if (a.eq.49.and.z.eq.22) then ! Ti49
          rsprc=1.10d0
        else if (a.eq.67.and.z.eq.30) then ! Zn67
          rsprc=1.08d0
        else if (a.eq.73.and.z.eq.32) then ! Ge73
          rsprc=1.13d0
        else if (a.eq.91.and.z.eq.40) then ! Zr91
          rsprc=0.947d0
        else if (a.eq.99.and.z.eq.44) then ! Ru99
          rsprc=0.952d0
        else if (a.eq.101.and.z.eq.44) then ! Ru101
          rsprc=0.952d0
        else if (a.eq.111.and.z.eq.48) then ! Cd111
          rsprc=0.873d0
        else if (a.eq.113.and.z.eq.48) then ! Cd113
          rsprc=0.867d0
        else if (a.eq.115.and.z.eq.50) then ! Sn115
          rsprc=0.863d0
        else if (a.eq.117.and.z.eq.50) then ! Sn117
          rsprc=0.863d0
        else if (a.eq.129.and.z.eq.54) then ! Xe129
          rsprc=0.857d0
        else if (a.eq.131.and.z.eq.54) then ! Xe131
          rsprc=0.918d0
        else if (a.eq.155.and.z.eq.64) then ! Gd155
          rsprc=0.884d0
        else if (a.eq.157.and.z.eq.64) then ! Gd157
          rsprc=0.884d0
        else if (a.eq.183.and.z.eq.74) then ! W183
          rsprc=0.873d0
        else if (a.eq.199.and.z.eq.80) then ! Hg199
          rsprc=0.870d0
        else if (a.eq.201.and.z.eq.80) then ! Hg201
          rsprc=0.870d0
        else if (a.eq.207.and.z.eq.82) then ! Pb207
          rsprc=0.865d0
        else
c            if (prtlevel.gt.0) write (*,*)
c     &           'dsddspsmff: OddG gaussian form factor unavailable ',
c     &           'for A=',a,' Z=',z
c            l2jjpp=-1.d0
c            l2jjnn=-1.d0
c            l2jjpn=-1.d0
c            ierr=1
c            return
          if (prtlevel.gt.0) then
            write (*,*)
     &           'dsddspsmff: SimpleOddOdd spin to charge radius',
     &           ' relation not found ',
     &           'for A=',a,' Z=',z
            write (*,*) 
     &           '   will assume spin radius is the same as charge radius'
          endif
          rsprc=1.d0            ! assume spin radius is equal to charge radius
        endif
        r=0.91d0*exp(log(dble(a))/3.d0)+0.3d0
        y=q*r*rsprc*fermiGeV
        ff=exp(-y*y/3.d0)
        l2jjpp=l2jjpp*ff
        l2jjnn=l2jjnn*ff
        l2jjpn=l2jjpn*ff
      endif

      return

 1000 continue
      if (prtlevel.gt.0) write (*,*)
     &     'dsddsimsdff: SimpleOddOdd SD form factor unavailable for A=',
     &     a,' Z=',z
      l2jjpp=-1.d0
      l2jjnn=-1.d0
      l2jjpn=-1.d0
      ierr=1
      return
      end
