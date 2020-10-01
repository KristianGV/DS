
* FOSubtract.h
* the fixed-order subtraction terms
* this file is part of FeynHiggs
* generated by conv/convSubtr.m 14-Sep-2015 10:20


        g1uC = gyMT*SB*
     &    (1 - cL/80.D0*(15*(2 - 7*CB2)*gMT2 - 
     &         (4*(45*htMT2 + gyMT2*(-44*SB2 + 7/4.D0*S2B**2)))/SB2))

#ifdef DETAILED_DEBUG
	DHIGGS "g1uC =", g1uC ENDL
#endif

        g1dC = CB*gyMT*
     &    (1 - cL/80.D0*(gMT2*(30 - 105*SB2) + gyMT2*(176 - 28*SB2)))

#ifdef DETAILED_DEBUG
	DHIGGS "g1dC =", g1dC ENDL
#endif

        g2uC = gMT*SB*(1 + 
     &      cL/240.D0*(540*htMT2 + gyMT2*(21*S2B**2 - 24*SB2) - 
     &          gMT2*(160*SB2 + 165/4.D0*S2B**2))/SB2)

#ifdef DETAILED_DEBUG
	DHIGGS "g2uC =", g2uC ENDL
#endif

        g2dC = CB*gMT*(1 - 
     &      cL/240.D0*(gyMT2*(24 - 84*SB2) + gMT2*(160 + 165*SB2)))

#ifdef DETAILED_DEBUG
	DHIGGS "g2dC =", g2dC ENDL
#endif

        sublog2L = cL**2/2.D0*
     &    (htMT2**2*vev**2*(tSUSY - tTop)*
     &      (g3MT**2*(64 - 96*xOS) - 
     &        3*htMT2*(2 + xOS)*(10 - (8 - xOS)*xOS) - 
     &        (96*g3MT**2 - 18*htMT2)*(tSUSY - tTop)))

#ifdef DETAILED_DEBUG
	DHIGGS "sublog2L =", sublog2L ENDL
#endif

        sublog1L = 1/(24.D0*Pi)*
     &    (Alfa1L*(2*(10*C2B**2*MW2*MZ2 - 
     &            4*MW2**2*(8 - 5*S2B**2) - MZ2**2*(5 + S2B**2))*
     &          (tCha - tTop) + 
     &         (36*MTy2**2 + (18*C2B*MTy2 - 84*C2B**2*MW2)*MZ2 + 
     &            MW2**2*(68 - 62*S2B**2) + 
     &            3*MZ2**2*(14 - S2B**2*(10 + 3*S2B**2)))*
     &          (tSUSY - tTop)))/(MW2*SW2)

#ifdef DETAILED_DEBUG
	DHIGGS "sublog1L =", sublog1L ENDL
#endif

	g1udC = g1dC + g1uC

	g2udC = g2dC + g2uC

	htC = 1/12.D0*(htMT*(12 - cL*(g1udC**2 + 3*g2udC**2)))

#ifdef DETAILED_DEBUG
	DHIGGS "htC =", htC ENDL
#endif

        lC = 1/12.D0*((gMT2 + gyMT2)*
     &       (3*C2B**2*(1 + cL*((gMT2 + gyMT2)*S2B**2)) + 
     &         9*cL*(C2B*htC**2*xOS)) - 
     &      cL*(6*gMT2*gyMT2 + 3*gyMT2**2 + 
     &         gMT2**2*(7 + 2*S2B**2) - 6*htC**4*(12 - xOS)*xOS))

#ifdef DETAILED_DEBUG
	DHIGGS "lC =", lC ENDL
#endif

        subnonlog = 1/12.D0*
     &    (vev**2*(3*C2B**2*(gMT2 + gyMT2) - 
     &        cL*(7*(g1dC**4 + g1uC**4) + 
     &           16*g1dC*g1uC*(g1dC**2 + g1uC**2 + g2udC**2) + 
     &           2*(g1uC**2*g2udC*(6*g2uC + g2udC) + 
     &              g1dC**2*(9*g1uC**2 + 6*g2dC*g2udC + g2udC**2))-
     &             4*g1udC**2*lC + 
     &           gMT2**2*(7 + (2 + 9*C2B**2)*S2B**2) + 
     &           3*(g2udC**2*
     &               (g2dC*(9*g2dC - 2*g2uC) + 9*g2uC**2 - 4*lC) + 
     &              gyMT2*(2*gMT2 + gyMT2)*(1 + 3*C2B**2*S2B**2))-
     &             ((9 - 3*C2B)*C2B*(gMT2 + gyMT2)*htC**2 + 
     &              htC**4*(72 - 6*xOS))*xOS)))

#ifdef DETAILED_DEBUG
	DHIGGS "subnonlog =", subnonlog ENDL
#endif
