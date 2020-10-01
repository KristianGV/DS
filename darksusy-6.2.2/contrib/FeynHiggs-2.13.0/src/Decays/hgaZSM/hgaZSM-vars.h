#if 0
* hgaZSM-vars.h
* variable declarations
* generated by FormCalc 8.4 on 27-Feb-2015 17:19
* this file is part of FeynHiggs
#endif

#ifndef VARS_H
#define VARS_H

#include "externals.h"
#include "types.h"
#include "debug.h"

#else

#include "Decay.h"
#include "hgaZSM-renconst.h"

	ComplexType Sub11, Sub10, Sub14, Sub17, Sub16, Sub20, Sub1(3)
	ComplexType Sub4(3), Sub7(3), Sub2(3), Sub5(3), Sub8(3)
	ComplexType pave1, pave2, pave3, pave4(Ncc), pave5(3)
	ComplexType pave6(Ncc,3), pave7(3), pave8(Ncc,3), pave9(3)
	ComplexType pave10(Ncc,3)
	common /hgaZSM_varXs/ Sub11, Sub10, Sub14, Sub17, Sub16
	common /hgaZSM_varXs/ Sub20, Sub1, Sub4, Sub7, Sub2, Sub5
	common /hgaZSM_varXs/ Sub8, pave1, pave2, pave3, pave4
	common /hgaZSM_varXs/ pave5, pave6, pave7, pave8, pave9
	common /hgaZSM_varXs/ pave10

	ComplexType Pair1, Pair2, Pair3, Eps1, Abb6, Abb10, Abb9, Abb11
	ComplexType Abb8, Abb5, Abb7, Abb1, Abb2, Abb4, Abb3, Sub23
	ComplexType Sub12, Sub19, Sub26, Sub21, Sub24, Sub18, Sub15
	ComplexType Sub13, Sub25, Sub22, Sub3(3)
	ComplexType Sub6(3), Sub9(3)
	common /hgaZSM_varXh/ Pair1, Pair2, Pair3, Eps1, Abb6, Abb10
	common /hgaZSM_varXh/ Abb9, Abb11, Abb8, Abb5, Abb7, Abb1
	common /hgaZSM_varXh/ Abb2, Abb4, Abb3, Sub23, Sub12, Sub19
	common /hgaZSM_varXh/ Sub26, Sub21, Sub24, Sub18, Sub15
	common /hgaZSM_varXh/ Sub13, Sub25, Sub22, Sub3, Sub6, Sub9

	integer seq(2), Hel(3)
	common /hgaZSM_helind/ seq, Hel

	integer Gen4
	common /hgaZSM_indices/ Gen4

	ComplexType Cloop(1)
	common /hgaZSM_formfactors/ Cloop

#endif