*******************************************************************
*** subroutine dsib2ffZsVIB add to the common block array amp(i,j) 
*** the helicity amplitudes of the s channel VIB diagrams. 
*** The Z is emitted from the mediators (H3, Z), becoming H_{i}, i= 1, 2. H_{i}.
*** The sum is performed only over ih.
*** 
*** Author: Francesca Calore, 2013-04-10
*******************************************************************
      subroutine dsib2ffZsVIB(f, ih) 
      implicit none
      include 'dsmssm.h'
      include 'dsib2com.h'

      integer f, ih

      integer kh(2),i 
      real*8 Mz, Mh3, Mw, Mhk(3)

      complex*16 deltaH3, deltaZ
           
      kh(1)=kh1
      kh(2)=kh2

c...  masses/mx mediators
      Mz = mass(kz)/mx
      Mh3 = mass(kh3)/mx
      do i= 1,2
         Mhk(i) = mass(kh(i))/mx
      end do
      Mw = mass(kw)/mx


c... denominators
c... H3 
      deltaH3 = (dcmplx(4*E1J**2 - Mhk(ih)**2,Mhk(ih)*width(kh(ih))/mx)*
     -     dcmplx(4*E1J*(E1J + EvJ) - Mh3**2 + MB**2,Mh3*width(kh3)/mx))
c... Z
      deltaZ = (dcmplx(4*E1J**2 - Mhk(ih)**2,Mhk(ih)*width(kh(ih))/mx)*
     -     dcmplx(4*E1J*(E1J + EvJ) - Mz**2 + MB**2,Mz*width(kz)/mx))


c... helicity amplitudes
      amp(0, 2) = amp(0, 2) + (
     -     (((0,16)*E1J*kJ*kv*gl(kz,kh3,kh(ih))*
     -     (Conjg(gr(kh3,kn(1),kn(1))) - gr(kh3,kn(1),kn(1)))*
     -     gr(kh(ih),f,f))/MB)/deltaH3) ! H3 contribution
     -     + ((                 ! Z contribution
     -     ((0,16)*E1J**2*(E1J + EvJ)*kJ*kv*Mw*
     -     (Conjg(gl(kz,kn(1),kn(1))) 
     -     + gl(kz,kn(1),kn(1)))*gl(kh(ih),kz,kz)*gr(kh(ih),f,f))
     -     /MB**3)/deltaZ)

      return
      end


