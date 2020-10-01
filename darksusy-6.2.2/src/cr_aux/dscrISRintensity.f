*******************************************************************************
*** Function dscrISRintensity provides the differential intensity of the    ***
*** local interstellar cosmic rays (LIS). In the current implementation,    ***
*** use the parameterzations by Della Torre et al. (1701.02363) and         ***
*** Boschini et al. (ApJ 840:115, 2017).                                    ***
***                                                                         ***
***  Input:                                                                 ***
***    r  - rigidity of CR partisle [GV]                                    ***
***    CRtype - 1 for protons                                               ***
***             2 for helium                                                ***
***                                                                         ***
***                                                                         ***
***  Units of output: 1/ [m**2 s sr GV]                                     ***
***                                                                         ***
*** author: Torsten.Bringmann.fys.uio.no                                    ***
*** date 2018-06-22                                                         ***
*******************************************************************************
      real*8 function dscrISRintensity(r, CRtype)
      implicit none
      include 'dsio.h'

      real*8 r, res
      integer CRtype
      
      res=0.0d0
      
      if (CRtype.eq.1) then ! protons
        if (r.gt.1) then
           res = res + 10500. ! original value 10800; 10500 connecte the two
     &                        ! regimes more smoothly 
     &           + 8590./r - 4.23d6/(3190.+r) + 2.74d5/(17.4+r) 
     &           - 39400./(0.464+r)
        elseif (r.gt.0.2) then
           res = res + 94.1 - 831*r + 16700*r**3 -10200*r**4
        else ! extrapolate by hand
           res = res + 8.61d4*r**4.7        
        endif

      elseif (CRtype.eq.2) then ! helium; updated 2019-02-01

        if (r.gt.1.9) then
           res = res + 3120. - 5530./r + 3370./(1.29+r) + 1.34d5/(88.5+r) 
     &           - 1.17d6/(861.+r) +0.03*r
        elseif (r.gt.0.2) then
           res = res + 1.14 - 118*r**2 + 578*r**3 - 87*r**5        
        else ! extrapolate by hand
           res = res + 0.0d0     
        endif

      else
        if (prtlevel.gt.1) then
          write(*,*) 'WARNING: dscrISRintensity called with' 
          write(*,*) 'unrecognized parameter CRtype = ',CRtype 
          write(*,*) 'I will return zero...'
        endif
      endif 

      res=res/r**2.7

c... added TB  2019-03-15     
      if (r.gt.4.0d6) then
        if (prtlevel.gt.1) then
          write(*,*) 'WARNIGN in dscrISRflux: entering CR knee region'
          write(*,*) 'Model this with softening of CR spectrum by Dn=0.6'
          write(*,*) 'CRtype, R = ', CRtype, r
        endif
        res=res/(r/4.0d6)**0.6
        if (r.gt.5.0d10) then ! ~ GZK
           res=res/(r/5.0d10)**4 ! 0.0d0 leads to numerical problems
                                 ! but any strong cut will do
        endif
      endif
      
      dscrISRintensity=res          
      return
      end


