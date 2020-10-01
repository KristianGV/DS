
c...
c...auxiliary function for inner integrand
c...Compared to dssenu_csint2.f, this routine integrates numerically over
c...the form factors
c...input: velocity relaitve to sun in km/s
c...output: integrand in cm^-4
c...We here follow the analysis in Gould, ApJ 321 (1987) 571 and more
c...specifically the more general expressions in appendix A, but further
c...generalized to allow for cut-off velocities at e.g. Jupiter and to 
c...allow for numerial integration over form factors
      real*8 function dssenu_csintff2(u,foveru)
      implicit none
      include 'dssem_sun.h'
      include 'dssecom.h'
      include 'dsmpconst.h'

      real*8 u,mu,mx,sigsi,sigsd,muplus,muminus,ma
      real*8 r,ni,v,dssem_sunvesc,vp,
     &  dssem_sundenscomp,foveru
      real*8 dssenu_csintff3,omvw,ww
      real*8 lo,hi,res
      integer vtype
      external foveru,dssenu_csintff3
      common/seint/mx,ma,sigsi,sigsd,r,vtype

      mu=mx/ma
      muplus=(mu+1.d0)/2.d0
      muminus=(mu-1.d0)/2.d0
      v=dssem_sunvesc(r/100.0d0)   ! escape velocity at r
c...now correct escape velocity so that we demand a slightly
c...lower escape velocity after to make sure it does not reach Jupiter
      vp=sqrt(v**2-veout**2)

      ni=dssem_sundenscomp(r/100.0d0,zx,ix) ! target number density in sun

      ww=sqrt(u**2+v**2)
      wwx=ww ! transfer to common needed by dssenu_csintff3

c...Now integrate the cross section with form-factor over momentum
c...transfer. We will use y=(delta-E)/(mx*w*2/2) as integration variable 
c...What we are integrating now is essentiall (A5) in Gould, but with 
c...sigma in the integral and with general form factors. 
c...We are adding ni*w later (or actaully w**2 as that is what we need).
      lo=(u**2+veout**2)/ww**2
      hi=mu/muplus**2
      call dshiprecint3(dssenu_csintff3,lo,hi,res)
      omvw=ni*(ww/c_light)**2*res ! this is (A5) in Gould *w

c...In the expressions abve, we have modified the expressions in Gould
c...(or if you prefer Lundberg & Edsjo, 2004) to allow for scatterings 
c...not only to velocities less than the escape velocity, but
c...instead require a more conservative value from requiring that
c...the WIMP does not reach out further than to a distance where the
c...escape velocity is veout. (veout=0 gives us back the old 
c...Gould expressions). This modification shows up in the lower limit
c...for the integrand.

      dssenu_csintff2=foveru(u)*omvw ! gould (A6) + (2.8) + modifications

c...we have two factors of c missing, add them to get units cm^-4
      dssenu_csintff2=dssenu_csintff2*(3.0d10)**2
c...dssenu_csintff2 is now equal to f(u)/u * w Omega_v^-(w), i.e.
c...f(u)/u * (A6) in Gould, except that we have calculated it for more
c...general form factors and allowing for a different escape velocity (i.e.
c...not reaching out to infinity).

c...Note, in principle we should have a theta function here to allow
c...only such scatterings that go to a velocity lower than the escape
c...velocity, theta(mu/muplus^2 - u^2/(u^2+v^2)),  however, this
c...is taken care of in the limits for the u integration and is thus not
c...needed here.
c      if (u**2/(u**2+v**2).gt.mu/muplus**2) then
c        dssenu_csintff2=0.0d0
c      endif

      return
      end
