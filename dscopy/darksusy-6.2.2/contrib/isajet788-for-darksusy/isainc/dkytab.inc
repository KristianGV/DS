C          LOOK must be dimensioned to the maximum value of INDEX.
      INTEGER   MXLOOK
      PARAMETER (MXLOOK=500)
      INTEGER   MXDKY
      PARAMETER (MXDKY=3000)
      COMMON/DKYTAB/LOOK(MXLOOK),CBR(MXDKY),MODE(5,MXDKY),MELEM(MXDKY)
      SAVE /DKYTAB/
      INTEGER   LOOK,MODE,MELEM
      REAL      CBR
#ifdef LEVEL2_X
      LEVEL 2,/DKYTAB/
#endif
