      INTEGER   LIMPOM
      PARAMETER (LIMPOM=20)
      COMMON/MBGEN/POMWT(LIMPOM),POMGEN(LIMPOM),MNPOM,MXPOM,PDIFFR,
     $NPOM,XBARY(2),DXBARY(2),XPOM(LIMPOM,2)
      SAVE /MBGEN/
      INTEGER   MNPOM,MXPOM,NPOM
      REAL      POMWT,POMGEN,PDIFFR,XBARY,DXBARY,XPOM