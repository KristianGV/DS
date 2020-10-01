***********************************************************************
*** note. this routine uses the full expressions for the capture
*** rate in the earth from gould, apj 521 (1987) 571.
*** this routine replaces dssenu_capearth which use the approximations
*** given in the jkg review.
***********************************************************************


       real*8 function dssenu_capearth2(mx,rho,sigsi)
c----------------------------------------------------------------------
c         capture rate in the earth
c       uses the full routines instead of jkg (as in dssenu_capearth).
c *** full: use formulas by gould as reported in jkg
c
c       mx: WIMP mass
c       sigsi: spin independent cross section in units of cm^2
c       vbar: 3D WIMP velocity dispersion in the halo
c       vstar: Sun's velocity through the halo
c       lars bergstrom 1998-09-15
c----------------------------------------------------------------------
       implicit none
       include 'dshmcom.h' ! needed for v_sun and vd_3d
       real*8 mx,rho,sigsi,dssenu_capearthfull
       real*8 v_star,v_bar
       v_star=v_sun
       v_bar=vd_3d
       dssenu_capearth2=dssenu_capearthfull(mx,sigsi,v_star,v_bar,rho)
       return
       end





