c_______________________________________________________________________
c  Function dsfisigma
c     Calculates cross section as function of s.
c
c
c
c
c  author: Kristian Gjestad Vangsnes (kristgva@uio.no)   2020-09-07
c=======================================================================
      real*8 function sigma(s)
      implicit none

      include 'dsficom.h'
      include 'dsmpconst.h'

      real*8 s,p,v,W,dsanwx,dssigmavpartial,dssigmav
      real*8, external :: v_ij
      p=sqrt(s/4-mdm**2)

c       vmoeller = 2.0d0*p*sqrt(s)/(s-2.0d0*dsmwimp()**2)


      
      if(isnan(dssigmavpartial(ichannel_22,p)/v)) then
            sigma=0
      else
c           dssigmavpartial is given in cm^3 s^-1            
            sigma=dssigmavpartial(ichannel_22,p)/v_ij(mdm,mdm,s)/gev2cm3s
            ! sigma=dssigmavpartial(ichannel_22,p)/v_ij(m1_22,m2_22,s)
      end if
      return
      end

      real*8 function v_ij(mi,mj,s)
      implicit none
      
      real*8 mi,mj,s

      v_ij=sqrt(s-(mi+mj)**2)*sqrt(s-(mi-mj)**2)/(s-mi**2-mj**2)

      return
      end