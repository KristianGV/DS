c Dilogarithm from the CHAPLIN package 1106.5739
      
c---  Li2

      double complex  function cli2(z)
      implicit none
      double complex ris, z, bsli2_inside,bsli2_outside, wcli2
      double complex zlocal
      double precision zabs, pi, zeta2, border, tiny, arg

      pi=3.1415926535897932385D0
      zeta2=pi**2/6d0

      border = 0.3d0 
      tiny = 1d-14
      zabs = abs(z)
      zlocal=z

      if (zabs.gt.1d0+tiny) then
         ris=-wcli2(1d0/z)-zeta2-0.5d0*log(-z)**2
      elseif (zabs.le.border) then 
         ris=bsli2_inside(z)
      else
         if (zabs.gt.1d0) then
            arg=atan2(dimag(zlocal),dreal(zlocal))
            zlocal=dcmplx(cos(arg),sin(arg))
         endif
         ris=bsli2_outside(zlocal)
      endif

      cli2=ris
      return
      end
      
c---  recursion
      
      double complex  function wcli2(z)
      implicit none
      double complex z, cli2
      wcli2 =  cli2(z)
      return
      end

C---- expansion of dilogarithm in y = - log(1-z) with Bernoulli numbers  

      double  complex function bsli2_inside(z)
      implicit none
      integer i, Nmax
      double complex ris, z, zb 
      double precision bern(11)

c     bern(i+1) = BernoulliB(2i)/(2i)!
      data bern /1.D0,0.8333333333333333D-1,-0.1388888888888889D-2
     &,0.3306878306878307D-4,-0.8267195767195767D-6,0.208767569878681D-7
     &,-0.5284190138687493D-9,0.1338253653068468D-10
     &,-0.3389680296322583D-12,0.8586062056277845D-14
     &,-0.2174868698558062D-15/
      parameter (Nmax=11)       ! this is half the order we want (coz odd bernoulli numbers are zero except BernoulliB(1)=-0.5d0)
      
      zb = dcmplx(1d0,0d0)-z
      zb = -log(zb)
      ris = -zb**2/4d0          !accounting for BernoulliB(1) = -0.5d0
      do i=1,Nmax
         ris = ris + zb**(2*i-1)*bern(i)/(2*i-1)
      enddo
      bsli2_inside = ris
      return 
      end

C---- expansion of the dilogarithm in log(z) with Zeta values  
C-------- used for border < |z| < 1
      
      double  complex function bsli2_outside(z)
      implicit none
      integer i, Nmax
      double complex ris, z, zb
      double precision zeta(29),zeta0,zeta2
      
c     zeta(i) = Zeta(2-2i-1)/(2i+1)! i.e. Zeta(-1)/6, Zeta(-3)/120, Zeta(-5)/7!....
      data zeta /-0.01388888888888889d0,0.00006944444444444444d0
     &,-7.873519778281683d-7,1.148221634332745d-8,-1.897886998897100d-10
     &,3.387301370953521d-12,-6.372636443183180d-14,1.246205991295067d-
     &15,-2.510544460899955d-17,5.178258806090624d-19,-1.088735736830085
     &d-20,2.325744114302087d-22,-5.035195213147390d-24,1.10264992943812
     &2d-25,-2.438658550900734d-27,5.440142678856252d-29,-1.222834013121
     &735d-30,2.767263468967951d-32,-6.300090591832014d-34,1.44208683884
     &1848d-35,-3.317093999159543d-37,7.663913557920658d-39,-1.777871473
     &383066d-40,4.139605898234137d-42,-9.671557036081102d-44,2.26671870
     &1676613d-45,-5.327956311328254d-47,1.255724838956433d-48,-2.967000
     &542247094d-50/

      parameter (Nmax=29) ! this is half the order we want (coz even zetaval2 are zero except for 0,2)
      parameter (zeta0 = 1.644934066848226d0)
      parameter (zeta2 = -0.2500000000000000d0)

      zb = log(z)
      ris = dcmplx(zeta0, 0d0) + zb*(1d0 -log(-zb)) 
     &     + zb**2*zeta2
      do i=1,Nmax 
         ris = ris + zb**(2*i+1)*zeta(i)
      enddo
      
      bsli2_outside=ris 
      return 
      end

