* Resum4H23.F
* resum the 4H coupling taking into account four scales:
*   low scale = MT
*   medium scale 1 = MCha
*   medium scale 2 = MGlu
*   high scale = MSUSY
* this file is part of FeynHiggs
* last modified 11 Feb 17 th

#include "externals.h"
#include "types.h"

#define loglevelXt 2


* Main function calculating the correction to Mh due to resummation

	subroutine Resum4H23(res)
	implicit none
	RealType res

#include "Resum4H.h"
#define __SUBROUTINE__ "Resum4H23"

	integer fail, i
	RealType x(nSM)
	RealType yOS, msMS, logarg
	RealType g1uC, g1uC2, g1dC, g1dC2
	RealType g2uC, g2uC2, g2dC, g2dC2
	RealType htC, htC2, lC
	RealType gMs, gMs2, gyMs, gyMs2
	RealType htMT2sub, sublog1L, sublog2L
	RealType subnonlog1L, subnonlog2La, subnonlog2Lb
	RealType dlam_asatOS, Mh2
	RealType dup1, dup2

	external threshdiff

	RealType ggyMT2
	parameter (ggyMT2 = 1/24D0*gMT2*gyMT2*
     &    (gyMT2 - 25*gMT2)/(gMT2 + gyMT2)**2)

	yOS = xOS + mueOS*(1/TB + TB)

	dlam_asatMS = 0

	if( loglevelXt .eq. 2 ) then
	  xMS = ggyMT2 + .75D0*htMT2*xOS1
	  xMS = xOS*(1 + k1L*(4*g3MT2 - xMS)*llog)
	  mueMS = mueOS
	  if( loglevel .gt. 2 ) then
	    call loopfun4H(lfmueMS, mueMS)
	    call TLthresholdasat(dlam_asatMS, xMS*MSUSYOS/MGl)
	  endif
	else
	  call loopfun4H(lfmueOS, mueOS)
	  xMS = 32/3D0*g3MT2 + (6 - pi*sqrt3)*htMT2*yOS/TB2
	  logarg = abs(xOS)*Mf(tT2,3)/MSUSYOS
	  if( logarg .gt. 1D-12 ) xMS = xMS + xOS*(
     &      4/3D0*g3MT2*(5 + 3*llog - xOS) -
     &      ggyMT2*llog -
     &      .75D0*htMT2/SB2*(
     &        2*lfmueOS(4) - mueOS**2*lfmueOS(2) - 3 +
     &        (llog - 35/6D0)*SB2 +
     &        (4/3D0*log2 + 2*log(logarg))*SB2*xOS2 -
     &        (4/9D0*pi*sqrt3 - 2)*CB2*yOS**2 )
     &      )
	  msMS = 1 + k1L*(-16/3D0*g3MT2 +
     &      .75D0*htMT2/SB2*(lfmueOS(5) + CB2 +
     &        2*SB2*xOS2 - (pi/sqrt2 - 2)*CB2*yOS**2))
	  xMS = (xOS + k1L*xMS)/msMS
	  tSUSYMS = tSUSYOS + log(msMS**2)
	  mueMS = mueOS/msMS
	  call loopfun4H(lfmueMS, mueMS)
	  call TLthresholdasat(dlam_asatMS, xMS*msMS*MSUSYOS/MGl)
	endif
	xMS2 = xMS**2

* evaluate loopfunctions necessary for threshold corrections
	call loopfunSf(lfSf, MSS(3,3)/MSS(4,3))
	call loopfunIno(lfM12, Re(M_1)/Re(MUEc), Re(M_2)/Re(MUEc))

	do i = 1, 5
	  tSS(i) = 2/3D0*log(abs(MSS(i,1)*MSS(i,2)*MSS(i,3)))
	enddo
	tSQ = 2*log(abs(MSS(3,3)))
	tSU = 2*log(abs(MSS(4,3)))
	tA0 = log(MA02)
	tMUE = log(Sq(MUEc))
	tM_2 = log(Sq(M_2))

* initial values from tree-level relations
	x(ig) = gMT
	x(igy) = gyMT
	x(ig3) = g3MT
	x(iht) = htMT
	x(itb) = TB

* solve for x using Newton procedure
	call FindMinimum(fail, threshdiff, x, nSM - 1)

#include "FOSubtract.h"

	if( looplevel .eq. 0 ) then
	  Mh2 = x(ilam) + subnonlog2Lb 
	  if( loglevel .lt. 3 ) Mh2 = Mh2 + subnonlog2La
	  Mh2 = vev**2*Mh2
	  Mh2 = Mh2 - k1L/vev**2*(
     &      6*MT2*(Mh2 - 4*MT2)*Re(B0(Mh2,MT2,MT2)) +
     &      .5D0*(Mh2*(Mh2 - 4*MZ2) + 12*MZ2**2)*Re(B0(Mh2,MZ2,MZ2)) +
     &      (Mh2*(Mh2 - 4*MW2) + 12*MW2**2)*Re(B0(Mh2,MW2,MW2)) +
     &      Mh2*( 6*Re(A0(MT2)) +
     &        (7*MW2 - 4*MZ2)*Re(A0(MZ2))/(MZ2 - MW2) -
     &        3*MW2*Re(A0(Mh2))/(Mh2 - MW2) +
     &        (3*(Mh2/(Mh2 - MW2) - MW2/(MZ2 - MW2)) - 11)*Re(A0(MW2)) +
     &        4.5D0*Mh2*Re(B0(Mh2,Mh2,Mh2)) +
     &        .5D0*(Mh2 + MZ2 + 2*MW2 - 6*MT2) ) -
     &      4*(MZ2**2 + 2*MW2**2) )
	  res = Mh2
	else if( looplevel .eq. 1 ) then
	  res = sublog1L + subnonlog1L - subnonlog2Lb
	  if( loglevel .lt. 3 ) res = res - subnonlog2La
	  res = x(ilam) - res
	else
	  res = sublog1L + sublog2L + subnonlog1L
	  if( loglevel .ge. 3 ) res = res + subnonlog2La
	  res = x(ilam) - res
	endif

	if( debuglevel .gt. 4 ) then
	  DHIGGS "final values:"
	  DHIGGS "lam(MT)      =", x(ilam)	ENDL
	  DHIGGS "subnonlog1L  =", subnonlog1L	ENDL
	  DHIGGS "sublog1L     =", sublog1L	ENDL
	  DHIGGS "subnonlog2La =", subnonlog2La	ENDL
	  DHIGGS "subnonlog2Lb =", subnonlog2Lb	ENDL
	  DHIGGS "sublog2L     =", sublog2L	ENDL
	endif
	end


************************************************************************
* threshdiff is the function which is minimized by FindMinimum.
* It takes n=5 input parameters (the guesses for lambda(MT) and
* the gaugino-gaugino-higgs couplings g[12][ud](MCha)) and uses
* them to run all couplings up to Q = MSUSY.
* At Q = MSUSY it checks how well the threshold corrections are
* fulfilled and saves the discrepancies in diff(1..5).
* |diff|^2 is given back to FindMinimum.

	subroutine threshdiff(n, x, dx, dxmin)
	implicit none
	integer n
	RealType x(nSM), dx(nSM), dxmin

#include "Resum4H.h"

	RealType y(nSplit)
	RealType ySB, ySB2, yCB, yCB2, yC2B, yS2B
	RealType yg2, ygy2, yht2
	RealType yg1d2, yg2d2, yg1u2, yg2u2
	RealType dx2, yMS

	external betaSM2L, betaSplit2L

* read in guesses
	y(ig) = x(ig)
	yg2 = y(ig)**2
	y(igy) = x(igy)
	ygy2 = y(igy)**2
	y(ig3) = x(ig3)
	y(iht) = x(iht)
	yht2 = y(iht)**2
	y(itb) = x(itb)

	yCB2 = 1/(1 + y(itb)**2)
	yCB = sqrt(yCB2)
	ySB = y(itb)*yCB
	ySB2 = ySB**2
	yC2B = yCB2*(1 - y(itb))*(1 + y(itb))
	yS2B = 2*ySB*yCB

* threshold corrections at Q = MSUSY

	y(ig1d) = 1/48.D0*(yCB*(48 - k1L*(12*xMS2*yht2 + 
     &           9*yg2*(2 - (7 - 6*tA0 + 6*tSUSYMS)*ySB2) + 
     &           ygy2*(4*(33 - 2*tA0 + 42*tSUSYMS) - 
     &              3*(7 - 6*tA0 + 6*tSUSYMS)*ySB2 - 
     &              8*(3*tSS(1) + 6*tSS(2) + tSS(3) + 8*tSS(4) + 
     &                 2*tSS(5)))))*y(igy))

	y(ig2d) = 1/48.D0*(yCB*(48 - k1L*(6*(ygy2 + 2*xMS2*yht2) - 
     &           3*(7 - 6*tA0 + 6*tSUSYMS)*ygy2*ySB2 + 
     &           yg2*(32 - 8*tA0 + 104*tSUSYMS + 
     &              (33 - 42*tA0 + 42*tSUSYMS)*ySB2 - 24*tSS(1) - 
     &              72*tSS(3))))*y(ig))

	y(ig1u) = 1/48.D0*((48*ySB**2 + k1L*
     &          (12*yht2*(9 + 7*tSQ - 13*tSU + 6*tSUSYMS - 
     &               xMS2*ySB**2) - 
     &            ySB**2*(9*(2 - (7 - 6*tA0 + 6*tSUSYMS)*yCB2)*
     &                yg2 + 
     &               ygy2*(132 - tA0*(8 - 18*yCB2) + 
     &                  6*tSUSYMS*(28 - 3*yCB2) - 21*yCB2 - 
     &                  8*(3*tSS(1) + 6*tSS(2) + tSS(3) + 
     &                     8*tSS(4) + 2*tSS(5))))))*y(igy))/ySB

	y(ig2u) = 1/48.D0*((48*ySB**2 + k1L*
     &          (12*yht2*(9 - 9*tSQ + 3*tSU + 6*tSUSYMS - 
     &               xMS2*ySB**2) - 
     &            ySB**2*(3*(2 - (7 - 6*tA0 + 6*tSUSYMS)*yCB2)*
     &                ygy2 + 
     &               yg2*(32 - 8*tA0 + 104*tSUSYMS + 
     &                  (33 - 42*tA0 + 42*tSUSYMS)*yCB2 - 
     &                  24*tSS(1) - 72*tSS(3)))))*y(ig))/ySB

* lambda(MSUSY): tree-level     
	y(ilam) = 1/4.D0*(yC2B**2*(yg2 + ygy2))

* lambda(MSUSY): 1L DRbar -> MSbar
	y(ilam) = y(ilam) - 
     &    k1L/12.D0*((9 - 2*yC2B**2)*yg2**2 + 6*yg2*ygy2 + 3*ygy2**2)

* lambda(MSUSY): 1L heavy Higgs
	y(ilam) = y(ilam) - 
     &    k1L*(3/4.D0*(yC2B**2*(yg2 + ygy2)**2*yS2B**2) - 
     &       1/192.D0*((tA0 - tSUSYMS)*
     &          (53*yg2**2 + 42*yg2*ygy2 + 29*ygy2**2 - 
     &            4*(ygy2**2 + yg2*(7*yg2 + 6*ygy2))*(yC2B - yS2B)*
     &             (yC2B + yS2B) - 
     &            9*(yg2 + ygy2)**2*
     &             (yS2B**4 + yC2B**2*(yC2B**2 - 6*yS2B**2)))))

* lambda(MSUSY): 1L sfermions
	y(ilam) = y(ilam) - 
     &    k1L/12.D0*(yht2*(6*
     &           (yC2B*(3*tSUSYMS*yg2 - tSQ*(3*yg2 - ygy2) + 
     &                (-4*tSU + 3*tSUSYMS)*ygy2) - 
     &             6*(tSQ + tSU - 2*tSUSYMS)*yht2) - 
     &          3*xMS2*(yht2*(24*lfSf(1) - 2*xMS2*lfSf(2)) + 
     &             yC2B*(3*(ygy2*lfSf(3) + yg2*lfSf(4)) - 
     &                yC2B*(yg2 + ygy2)*lfSf(5)))) + 
     &       yC2B**2*(4*tSUSYMS*(3*yg2**2 + 5*ygy2**2) - 
     &          3*yg2**2*(tSS(1) + 3*tSS(3)) - 
     &          ygy2**2*(3*tSS(1) + 6*tSS(2) + tSS(3) + 8*tSS(4) + 
     &             2*tSS(5))))

* 2L correction to lambda(MSUSY)     
	if( loglevel .gt. 2 ) then
* asat threshold
	  y(ilam) = y(ilam) + 4*k2L*yht2**2*y(ig3)**2*dlam_asatMS
* atat threshold
	  yMS = xMS + mueMS*(1/y(itb) + y(itb))
	  y(ilam) = y(ilam) - k2L/12.D0*(yht2**2*
     &        (3*yht2*(6*xMS2**3*ySB2 + 
     &             12*(3 + 16*Kfac)*(yCB2*yMS*(4*xMS + yMS)) + 
     &             xMS2**2*
     &              (yCB2*(75 + (19 + 96*Kfac)*yMS**2) - 
     &                6*(11 - lfmueMS(1) + lfmueMS(2))) - 
     &             4*xMS2*(yCB2*
     &                 (45 + 12*xMS*yMS + 14*yMS**2 - 
     &                   Kfac*(72 - 48*xMS*yMS - 72*yMS**2)) - 
     &                3*(7 - 4*lfmueMS(1) + 4*lfmueMS(2))) - 
     &             6*(1 - (13 - 120*Kfac - 2*pi**2)*yCB2 - 
     &                16*lfmueMS(1) + 
     &                mueMS**2*
     &                 (2*(6 - 2*lfmueMS(1) - 
     &                      xMS2*
     &                       (6 + 6*lfmueMS(1) - 
     &                       xMS2*(1 + lfmueMS(1)))) + 
     &                   (6 - (12 - xMS2)*xMS2)*lfmueMS(2)) + 
     &                8*lfmueMS(3)))))/ySB2
     &       
	endif

* run from Q = MSUSY to Q = MCha with Split RGEs

	call ODESolve(nSplit, betaSplit2L, y,
     &    tSUSYMS, tCha, odeeps, odeh1, 0)

* threshold conditions @ MCha

	y(ig) = y(ig)*(1 - k1L/3.D0*(3*tCha - tMUE - 2*tM_2)*y(ig))
	yg2 = y(ig)**2

	y(igy) = y(igy)*(1 - k1L/3.D0*(tCha - tMUE)*y(igy))
	ygy2 = y(igy)**2

	yg1u2 = y(ig1u)**2
	yg1d2 = y(ig1d)**2
	yg2u2 = y(ig2u)**2
	yg2d2 = y(ig2d)**2

	y(iht) = (1 - 
     &      k1L/12.D0*((3*tCha - 3*tMUE)*
     &          (yg1d2 + yg1u2 + 3*(yg2d2 + yg2u2)) - 
     &         (yg1d2 + yg1u2)*lfM12(6,1) - 
     &         3*(yg2d2 + yg2u2)*lfM12(6,2) - 
     &         2*lfM12(5,1)*y(ig1d)*y(ig1u) - 
     &         6*lfM12(5,2)*y(ig2d)*y(ig2u)))*y(iht)

	y(ilam) = y(ilam) - 
     &    k1L/12.D0*(7*(yg1d2**2 + yg1u2**2)*lfM12(1,1) + 
     &       27*(yg2d2**2 + yg2u2**2)*lfM12(2,2) + 
     &       14*(yg1d2*yg2d2 + yg1u2*yg2u2)*lfM12(2,3) + 
     &       18*yg1d2*yg1u2*lfM12(3,1) + 
     &       2*(yg1u2*yg2d2 + yg1d2*yg2u2)*lfM12(3,3) + 
     &       42*yg2d2*yg2u2*lfM12(4,2) + 
     &       16*lfM12(4,3)*(y(ig1u)*y(ig2d) + y(ig1d)*y(ig2u))*
     &        (y(ig1d)*y(ig2d) + y(ig1u)*y(ig2u)) + 
     &       8*y(ig1d)*y(ig1u)*
     &        (4*lfM12(1,3)*y(ig2d)*y(ig2u) + 
     &          lfM12(5,1)*(2*(yg1d2 + yg1u2) - y(ilam))) + 
     &       24*lfM12(5,2)*y(ig2d)*y(ig2u)*
     &        (2*(yg2d2 + yg2u2) - y(ilam)) - 
     &       12*((1/3.D0*((yg1d2 + yg1u2)*lfM12(6,1)) + 
     &             (yg2d2 + yg2u2)*lfM12(6,2))*y(ilam) + 
     &          (tCha - tMUE)*
     &           (yg1d2**2 + yg1u2**2 + 
     &             2*(yg1u2 + yg2d2)*(yg1d2 + yg2u2) + 
     &             5*(yg2d2**2 + yg2u2**2) + 
     &             4*y(ig1d)*y(ig1u)*y(ig2d)*y(ig2u) - 
     &             (yg1d2 + yg1u2 + 3*(yg2d2 + yg2u2))*y(ilam))))

* run from Q = MCha to Q = MT with SM RGEs

	call ODESolve(nSM, betaSM2L, y,
     &    tCha, tTop, odeeps, odeh1, 0)

* return lambda(MT)
* (must do this explicitly since minimization does not extend over lambda)

	x(ilam) = y(ilam)

* calc difference to known values of SM couplings @MT        

	dx(ig3) = y(ig3) - g3MT
	dx2 = dx(ig3)**2

	dx(iht) = y(iht) - htMT
	dx2 = dx2 + dx(ig3)**2

	dx(ig) = y(ig) - gMT
	dx2 = dx2 + dx(ig)**2

	dx(igy) = y(igy) - gyMT
	dx2 = dx2 + dx(igy)**2

	dx(itb) = y(itb) - TB
	dx2 = dx2 + dx(itb)**2

	dxmin = dx2
	end
