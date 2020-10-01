c...The include below includes the standard templates for user-definable
c...functions. The templates are in src/templates/. If you want to use your
c...own functions, copy one of these templates and the templates-standard.f
c...and modify templates-standard.f to load your new modified user function.
      include 'templates_standard.f'

      program dstemp
c
c     This program is just for temporary testing.
c
c     In addition to the DarkSUSY programs and data files, this test
c     program needs the input file dstest.mod, provided with the
c     DarkSUSY distribution. The test program outputs one file for
c     checking purposes (in an actual application, the user should
c     decide what information to output and in which format):
c     dstest1.tmp is output for the first model. You should compare it
c     with dstest.tmp provided with the DarkSUSY distribution.
c     You will also find dstest.output with the distribution. This file is
c     the output (to the terminal) when running the test program. You
c     should compare this with your output as well. All numbers should
c     to within the numerical accuracy. Some numbers, like the Z gamma
c     cross sections may wary a few % when the cross sections are low
c     due to numerical differences.

c     
c-----This line is 72 columns long--------------------------------------
c
      implicit none
      real*8 oh2a,oh2b,xf,dsrdomega                      ! relic density
      real*8 sigsip,sigsin,sigsdp,sigsdn                 ! direct detection
      integer i,nnuc,a(10),z(10),stoich(10)              ! direct detection
      real*8 e,si(10),sd(10),siff(10),sdff(10),t,rsi,rsd ! direct detection
      real*8 phiep,dsepdiff                              ! positron flux
      real*8 fluxgacdiff,fluxgac,fluxgaga,fluxgaz        ! gamma-ray flux
      real*8 jpsi,cospsi0,delta,fdelta,dshmjave,dshmj    ! gamma-ray flux
      real*8 jfunc,epshx,eps1hx                          ! gamma-ray flux
      external jfunc                                     ! gamma-ray flux
      integer nminhx,nmaxhx,niterhx,ierrhx               ! gamma-ray flux
      real*8 nsigvgaga,nsigvgaz,nsigvgacont,nsigvgacdiff ! gamma-rays
      integer istat                                      ! gamma-ray flux, etc
      real*8 egam,dshrgacontdiff,egath,dshrgacont        ! gamma-ray flux
      real*8 tpbess(3),pb_a,pb_b,pb_c,dshrpbardiff       ! pbar flux
      real*8 eth,thmax,rateea,ratesu,energy,theta        ! neutrino telescopes
      integer rtype,ptype                                ! neutrino telescopes
      real*8 sigpred, bgpred, lnLike, pval, refLike, dof ! IceCube likelihood
      real*8 ICtheoryError,phi_cut, totobsdbl            ! IceCube likelihood
      integer totobs, likechoice                         ! IceCube likelihood
      logical uselogNorm, pvalFromRef                    ! IceCube likelihood
      logical BGLikePrecompute                           ! IceCube likelihood
      character (len=256) eventf, edispf, BGf, efareaf   ! IceCube likelihood
      real*8 gm2amp,dsgm2muon                            ! g-2 amplitude
      integer unphys,excl,hwarning,acceptable,iend,ierr,iwar,nfc
      real*8 dshrmuhalo,dshrmudiff,nmusigmav,phimuhalo,dsabsq,tp,phidb,
     &  dshrdbardiff,dsepspecm
      real*8 tkd,dskdtkd,dskdmcut,mcut
      real*8 l2jjpp,l2jjnn,l2jjpn
      real*8 enu, dsepdphidphaaxi
      real*8 phiin, rholoc, dsrhosph_def, dspbdphidthaaxi, dsdbdphidthaaxi
      character message*80
      logical first
      data first/.true./

      logical UNITOK, UNITOP
      
      
c TB: set up 'real' testing of DS
      integer testlevel, testunit
      data testlevel/1/           ! 1 -- [default] standard output to screen,
                                  !      warnings if results are different from savfile 
                                  ! 2 -- minimal output, ~only warning
                                  ! 9 -- rewrite savfile 
                                  !      [NB: always change savshort to something
                                  !           else than default for this option!!!]
      character*200 savfile, savshort
      data savshort/'DStest_sav.dat'/    ! default: 'DStest_sav.dat'
      data testunit/98/
      real*8 comp


c
c     Here we include the file dssusy.h which contains common block
c     definitions for particle masses, susy parameters, and results.
c     We also include dsidtag.h which gives the possibility to tag
c     every model with a 12 character id tag (variable idtag). This
c     id tag is printed in case of errors for that model.
c
      include 'dsmssm.h'
      include 'dsio.h'
      include 'dsidtag.h'

      

c
c     This call initializes the DarkSUSY package. This call should
c     be the first call in any program using DarkSUSY. This call initializes
c     some global variables and calls various other modules to initialize
c     them with their default values. Check out src/ini/dsinit.f if you
c     want to see what it does.
c
      call dsinit

      if (moduletag.ne.'MSSM') then
         write(*,*) 'Sorry, dstest is only designed for the MSSM',
     &     ' module and'
        write(*,*) 'you compiled it with module ',moduletag
        write(*,*)
        write(*,*) 'Set DS_MODULE=ds_mssm in /examples/test/makefile' //
     &              'and try again...'
        stop
      endif

      write (*,*) 'DarkSUSY initialized'

      call dsgive_model(500.d0,1000.d0,400.d0,10.d0,3000.d0,0.d0,0.d0)

      call dsmodelsetup(unphys,hwarning)


      write(*,*) '  Neutralino mass = ',mass(kn(1))
      write(*,*) '  Gaugino fraction = ',
     &    dsabsq(neunmx(1,1))+dsabsq(neunmx(1,2))
      write(*,*) '  H1 mass =  ',mass(kh1),width(kh1)
      write(*,*) '  H2 mass =  ',mass(kh2),width(kh2)
      write(*,*) '  H3 mass =  ',mass(kh3),width(kh3)
      write(*,*) '  H+- mass = ',mass(khc),width(khc)
      

      write(*,*) 'Calculating omega h^2 with coannihilations,',
     &     ' please be patient...'
      oh2b=dsrdomega(1,1,xf,ierr,iwar,nfc)
      write(*,*) '  with coannihilations Oh2 = ',oh2b,ierr,iwar
      
      end
