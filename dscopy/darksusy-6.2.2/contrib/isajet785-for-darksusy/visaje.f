CDECK  ID>, VISAJE.
      CHARACTER*40 FUNCTION VISAJE()
      COMMON/IDRUN/IDVER,IDG(2),IEVT,IEVGEN
      SAVE /IDRUN/
      INTEGER   IDVER,IDG,IEVT,IEVGEN
      VISAJE = ' ISAJET     V7.85   04-NOV-2015 13:42:47'
      IDVER = 785
      RETURN
      END
