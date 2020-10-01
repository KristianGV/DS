*     -*- mode: fortran -*-
      integer ddng,ddnsf,ddnme
      parameter (ddng=27)
      parameter (ddnsf=6)
      parameter (ddnme=9)
! sf : 1_m, 2_sigma, 3_delta, 4_phi, 5_deltasigma, 6_mphi
      character*16 ddsf(6)
! me : 1_s, 2_v, 3_vm, 4_t, 5_a, 6_am, 7_p, 8_t2e, 9_t2o
      character*16 ddme(9)
      real*8 fme(9,2,4) ! me ;  p,n ; u,d,s,g
      common /ddcom/
     &     fme,
     &     ddsf,
     &     ddme
      save /ddcom/

! added by TB -- quenching for Borexino etc
      real*8 lntqdat(300),lnTrecdat(300),lndTdTdat(300)
      integer nqdat,khi,klo,quenchhow
      logical quenching_set
      common /ddquench/ lntqdat,lnTrecdat,lndTdTdat,
     &                  nqdat,khi,klo,quenchhow,quenching_set
      save /ddquench/

! added by TB -- astro parameters, so far only used by dsddCR routines
      real*8 rholocal, vlocal, Deff
      integer attenuation_how
      common /DDCRastro/ rholocal, vlocal, Deff, attenuation_how
      save /DDCRastro/
c... Elements considered: /O,Si,Mg,Fe,Ca,P,Na,S,Ni,Al,Cr/
      integer Nelements
      parameter (Nelements=11)
      integer AN(Nelements), ZN(Nelements)
      real*8 mNaU(Nelements)
      data AN   /16,28,24,56,40,30,23,32,59,27,52/
      data ZN   /8,14,12,26,20,15,11,16,28,13,24/ 
      data mNaU /16.0d0,28.1d0,24.3d0,55.85d0,40.0d0,30.0d0,23.0d0,
     &           32.0d0,59.0d0,27.0d0,52.0d0/        ! mass in au

! obsolescent
      real*8 ftp(7:12),ftn(7:12),delu,deld,dels
      common /ddcomlegacy/
     &     ftp,ftn,delu,deld,dels
      save /ddcomlegacy/
