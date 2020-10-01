*******************************************************************************
*** Function dssisigtm provides the momentum-transfer cross section per DM  ***
*** mass, in conventional units of cm**2/g.                                 ***
***                                                                         ***
***  type : interface                                                       ***
***  desc : momentum-transfer cross section                                 ***
***                                                                         ***
***  Input:                                                                 ***
***    vkms - relative velocity of scattering DM particles (in km/s)        ***
***                                                                         ***
***  This interfacte function is required by dssisigtmav in src/.           ***
***                                                                         ***
*** author: Torsten.Bringmann.fys.uio.no                                    ***
*** date 2015-05-18                                                         ***
*******************************************************************************
      real*8 function dssisigtm(vkms)
      implicit none
      include 'dsvdSIDM.h'
      include 'dsmpconst.h'

      real*8 vkms

      real*8 alphaDM, res, vrel, beta, R
      real*8 dssisigmatborn, dssisigmatclassical, dssisigmatres

      dssisigtm = 0.0d0
      
      alphaDM = gDM**2/4./pi
      vrel = vkms/c_light
      
      if (alphaDM*mass(kDM)/mass(kmed).lt.0.3) then ! Born regime

         beta = 2*alphaDM*mass(kmed)/mass(kDM)/vrel**2
         R = vrel*mass(kDM)/mass(kmed) 
         res = dssisigmatborn(mass(kmed),beta,R)

      elseif (vrel*mass(kDM)/mass(kmed).gt.1.0) then ! classical regime

         beta = 2*alphaDM*mass(kmed)/mass(kDM)/vrel**2
         if (DMtype.eq.1.and.Mediatortype.eq.2) then ! DM fermion, vector mediator:
                                                     ! average repulsive and attractive part 
            res = (dssisigmatclassical(mass(kmed),beta,0)
     &            + dssisigmatclassical(mass(kmed),beta,1))/2.    

         elseif (DMtype.eq.1.and.Mediatortype.eq.3) then ! DM fermion, scalar mediator:
                                                         ! only attrative potential
            res = dssisigmatclassical(mass(kmed),beta,0)

         else
           write(*,*) 'ERROR in dssisigtm: called with unsupported '
           write(*,*) 'DM/DR/mediator spin combination: ',  
     &                DMtype, Mediatortype, DRtype   
           write(*,*) 'program stopping...'
           stop
         endif

      else ! use approximation for Hulthén potential

         if (DMtype.eq.1.and.Mediatortype.eq.2) then ! DM fermion, vector mediator:
                                                     ! average repulsive and attractive part 
            res = (dssisigmatres(mass(kDM),mass(kmed),vrel,alphaDM,0)
     &           + dssisigmatres(mass(kDM),mass(kmed),vrel,alphaDM,1))/2.   
              
         elseif (DMtype.eq.1.and.Mediatortype.eq.3) then ! DM fermion, scalar mediator:
                                                         ! only attrative potential
            res = dssisigmatres(mass(kDM),mass(kmed),vrel,alphaDM,0)

         else
           write(*,*) 'ERROR in dssisigtm: called with unsupported '
           write(*,*) 'DM/DR/mediator spin combination: ',  
     &                DMtype, Mediatortype, DRtype   
           write(*,*) 'program stopping...'
           stop
         endif

      endif
          
c... divide by DM mass and convert to cm**2/g
      dssisigtm = res/mass(kDM) * gev3cm2g         
      
      return
      end

