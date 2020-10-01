      
c####################################################################
********************************************************************
c####################################################################
      REAL*8 FUNCTION omega(FF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

*----------------------------------------------------
*OMEH2=CISLO*PLANCK*VOLNO*TEMP*(TGAMMA**3)/SJINT/CRIT

	CISLO=1.66
	VOLNO=9.

*PLANCK=0.818933748
* MP     = 2.1767E-08      kg  
* 1  eV  = 1.782661731E-36 kg
* 1 GeV  = 1.782661731E-27 kg
* MP     = 2.1767/1.782661731 *1.E+19 GeV
	PLANCK=1.782661731/2.1767    
*1.E-19 left

*TEMP=(T_x/T_gamma)^3=2./(2.+7./8.*(6.*2.+5.*3.*2.))
	TEMP=2./(2.+7./8.*(6.*2.+5.*3.*2.))

 

* k boltz= 1.3806503E-23 J/K
*        = 8.617342E-05 ev/K
*        = 8.617342E-14 GeV/K
* T_gamma^3 = (2.726*0.8617342)^3  ! 10^(-39) left
	TG3=(2.726*0.8617342)**3     
* 10^(-39) left

* ro_crit= 8.0992E-47 GeV^(4)
	CRIT = 8.0992                
* 10^47 left

*CS: GeV^(-2) =0.3893796623  *1.E+09 pb
*    1pb      =1./0.3893796623*1.E-09 GeV^2
*    1fb      =1./0.3893796623*1.E-12 GeV^2
*

	CSCONV = 3.893796623               
* 10^11 left
 
*  -19-39+47+11  =0 

        omega=10000
        IF(FF.ne.0.) then
	  omega=CISLO*PLANCK*VOLNO*TEMP*TG3/CRIT*CSCONV/1000./FF
	endif
      
      return
      
*------------------------------------------------------
      END
